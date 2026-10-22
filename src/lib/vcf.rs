//! VCF text written field by field: a header of meta-information lines and sample names, and records
//! whose values are formatted as htslib formats them, floats included.
use std::fmt::Write;

use ahash::{AHashMap, AHashSet};

/// The `##fileformat` line every header starts with.
const FILE_FORMAT: &str = "##fileformat=VCFv4.2";

/// The columns every header line of column names starts with.
const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";

/// The bits of the float written as `.`, a NaN no arithmetic produces.
const MISSING_FLOAT_BITS: u32 = 0x7F80_0001;

/// An integer or float value with a sentinel that is written as `.`.
pub trait Numeric: Copy {
    /// The value written as `.`.
    fn missing() -> Self;

    /// Whether this is the value written as `.`.
    fn is_missing(&self) -> bool;
}

impl Numeric for i32 {
    fn missing() -> Self {
        i32::MIN
    }

    fn is_missing(&self) -> bool {
        *self == i32::MIN
    }
}

impl Numeric for f32 {
    fn missing() -> Self {
        f32::from_bits(MISSING_FLOAT_BITS)
    }

    fn is_missing(&self) -> bool {
        self.to_bits() == MISSING_FLOAT_BITS
    }
}

/// One allele of an unphased genotype.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenotypeAllele {
    /// An allele by its index, 0 for REF.
    Unphased(i32),
    /// An allele not called, written as `.`.
    UnphasedMissing,
}

/// A VCF header: meta-information lines in the order they were pushed, after `##fileformat`, and
/// the sample names. A line with the key and ID of a structured line pushed before, or repeating an
/// unstructured line, is dropped, as htslib drops it.
#[derive(Clone, Debug, Default)]
pub struct Header {
    lines: Vec<String>,
    samples: Vec<String>,
    pushed: AHashSet<String>,
    contig_ids: AHashMap<String, u32>,
    info: AHashSet<String>,
    format: AHashSet<String>,
    filter: AHashSet<String>,
}

impl Header {
    /// Adds a sample column.
    pub fn push_sample(&mut self, name: &str) -> &mut Self {
        self.samples.push(name.to_string());
        self
    }

    /// Adds a meta-information line, such as `##source=vartovcf` or `##INFO=<ID=DP,...>`.
    pub fn push_record(&mut self, line: &str) -> &mut Self {
        let structured = structured_id(line);
        let identity = match structured {
            Some((key, id)) => format!("##{key}=<ID={id}"),
            None => line.to_string(),
        };
        if !self.pushed.insert(identity) {
            return self;
        }
        match structured {
            Some(("contig", id)) => {
                let rid = u32::try_from(self.contig_ids.len()).expect("Over 2^32 contigs!");
                self.contig_ids.insert(id.to_string(), rid);
            }
            Some(("INFO", id)) => {
                self.info.insert(id.to_string());
            }
            Some(("FORMAT", id)) => {
                self.format.insert(id.to_string());
            }
            Some(("FILTER", id)) => {
                self.filter.insert(id.to_string());
            }
            _ => {}
        }
        self.lines.push(line.to_string());
        self
    }

    /// The index of a contig among the `##contig` lines, or an error naming a contig they lack.
    pub fn contig_id(&self, contig: &str) -> Result<u32, String> {
        self.contig_ids
            .get(contig)
            .copied()
            .ok_or_else(|| format!("Contig '{contig}' has no ##contig line in the VCF header!"))
    }

    /// The header as VCF text, ending with the line of column names.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut text = String::new();
        for line in std::iter::once(FILE_FORMAT).chain(self.lines.iter().map(String::as_str)) {
            text.push_str(line);
            text.push('\n');
        }
        text.push_str(COLUMNS);
        if !self.samples.is_empty() {
            text.push_str("\tFORMAT");
            for sample in &self.samples {
                text.push('\t');
                text.push_str(sample);
            }
        }
        text.push('\n');
        text.into_bytes()
    }

    /// Refuses a field the header does not declare.
    fn check(&self, kind: &str, id: &str) -> Result<(), String> {
        let declared = match kind {
            "INFO" => &self.info,
            "FORMAT" => &self.format,
            _ => &self.filter,
        };
        if declared.contains(id) {
            Ok(())
        } else {
            Err(format!(
                "The {kind} field '{id}' is not declared in the VCF header!"
            ))
        }
    }
}

/// The key and ID of a structured meta-information line, `##key=<ID=id,...>`.
fn structured_id(line: &str) -> Option<(&str, &str)> {
    let (key, value) = line.strip_prefix("##")?.split_once('=')?;
    let value = value.strip_prefix("<ID=")?;
    let end = value.find([',', '>']).unwrap_or(value.len());
    Some((key, &value[..end]))
}

/// One VCF record, built field by field and written as one line; ID and QUAL are always `.`.
pub struct Record<'h> {
    header: &'h Header,
    rid: u32,
    site: String,
    filters: String,
    info: String,
    format: String,
    samples: Vec<String>,
}

impl<'h> Record<'h> {
    /// Starts a record of a contig the header declares, at a 1-based position, with its REF and ALT.
    pub fn new(
        header: &'h Header,
        contig: &str,
        position: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Self, String> {
        let rid = header.contig_id(contig)?;
        let mut site =
            String::with_capacity(contig.len() + ref_allele.len() + alt_allele.len() + 16);
        write!(
            site,
            "{contig}\t{position}\t.\t{ref_allele}\t{alt_allele}\t."
        )
        .unwrap();
        Ok(Record {
            header,
            rid,
            site,
            filters: String::new(),
            info: String::new(),
            format: String::new(),
            samples: vec![String::new(); header.samples.len()],
        })
    }

    /// The index of the record's contig among the header's `##contig` lines.
    pub fn rid(&self) -> u32 {
        self.rid
    }

    /// Adds a FILTER label. `PASS` replaces every label added before it and gives way to any label
    /// added after it, and a label added twice is written once.
    pub fn push_filter(&mut self, label: &str) -> Result<(), String> {
        self.header.check("FILTER", label)?;
        if label == "PASS" || self.filters == "PASS" {
            self.filters.clear();
        }
        if self.filters.split(';').any(|added| added == label) {
            return Ok(());
        }
        if !self.filters.is_empty() {
            self.filters.push(';');
        }
        self.filters.push_str(label);
        Ok(())
    }

    /// Adds a String INFO field; an empty value leaves the field out.
    pub fn push_info_string(&mut self, key: &str, value: &str) -> Result<(), String> {
        self.header.check("INFO", key)?;
        if !value.is_empty() {
            self.start_info(key);
            self.info.push_str(value);
        }
        Ok(())
    }

    /// Adds an Integer INFO field; no values leave the field out.
    pub fn push_info_integer(&mut self, key: &str, values: &[i32]) -> Result<(), String> {
        self.push_info(key, values, push_integer)
    }

    /// Adds a Float INFO field; no values leave the field out.
    pub fn push_info_float(&mut self, key: &str, values: &[f32]) -> Result<(), String> {
        self.push_info(key, values, push_float)
    }

    fn push_info<V: Numeric>(
        &mut self,
        key: &str,
        values: &[V],
        push: fn(&mut String, V),
    ) -> Result<(), String> {
        self.header.check("INFO", key)?;
        if !values.is_empty() {
            self.start_info(key);
            push_values(&mut self.info, values, push);
        }
        Ok(())
    }

    fn start_info(&mut self, key: &str) {
        if !self.info.is_empty() {
            self.info.push(';');
        }
        self.info.push_str(key);
        self.info.push('=');
    }

    /// Adds the GT FORMAT field from the alleles of each sample in turn; no alleles leave it out.
    pub fn push_genotypes(&mut self, alleles: &[GenotypeAllele]) -> Result<(), String> {
        self.push_format("GT", alleles, |sample, alleles| {
            for (index, allele) in alleles.iter().enumerate() {
                if index > 0 {
                    sample.push('/');
                }
                match allele {
                    GenotypeAllele::Unphased(allele) => push_integer(sample, *allele),
                    GenotypeAllele::UnphasedMissing => sample.push('.'),
                }
            }
        })
    }

    /// Adds an Integer FORMAT field from the values of each sample in turn; no values leave it out.
    pub fn push_format_integer(&mut self, key: &str, values: &[i32]) -> Result<(), String> {
        self.push_format(key, values, |sample, values| {
            push_values(sample, values, push_integer)
        })
    }

    /// Adds a Float FORMAT field from the values of each sample in turn; no values leave it out.
    pub fn push_format_float(&mut self, key: &str, values: &[f32]) -> Result<(), String> {
        self.push_format(key, values, |sample, values| {
            push_values(sample, values, push_float)
        })
    }

    /// Adds a FORMAT field, with an equal share of the values for each sample.
    fn push_format<V>(
        &mut self,
        key: &str,
        values: &[V],
        push: impl Fn(&mut String, &[V]),
    ) -> Result<(), String> {
        self.header.check("FORMAT", key)?;
        if values.is_empty() {
            return Ok(());
        }
        let samples = self.samples.len();
        if samples == 0 || !values.len().is_multiple_of(samples) {
            return Err(format!(
                "Expected a multiple of {samples} values for the FORMAT field '{key}', found {}!",
                values.len()
            ));
        }
        let first = self.format.is_empty();
        if !first {
            self.format.push(':');
        }
        self.format.push_str(key);
        let width = values.len() / samples;
        for (sample, values) in self.samples.iter_mut().zip(values.chunks_exact(width)) {
            if !first {
                sample.push(':');
            }
            push(sample, values);
        }
        Ok(())
    }

    /// The record as one line of VCF text, ending with a newline; without a FORMAT field it has no
    /// FORMAT or sample columns, as htslib writes it.
    pub fn to_vcf_string(&self) -> String {
        let mut line = String::with_capacity(
            self.site.len()
                + self.filters.len()
                + self.info.len()
                + self.format.len()
                + self.samples.iter().map(String::len).sum::<usize>()
                + 8,
        );
        line.push_str(&self.site);
        line.push('\t');
        line.push_str(if self.filters.is_empty() {
            "."
        } else {
            &self.filters
        });
        line.push('\t');
        line.push_str(if self.info.is_empty() {
            "."
        } else {
            &self.info
        });
        if !self.format.is_empty() {
            line.push('\t');
            line.push_str(&self.format);
            for sample in &self.samples {
                line.push('\t');
                line.push_str(sample);
            }
        }
        line.push('\n');
        line
    }
}

/// Appends values separated by commas, `.` for a missing one.
fn push_values<V: Numeric>(out: &mut String, values: &[V], push: fn(&mut String, V)) {
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        if value.is_missing() {
            out.push('.');
        } else {
            push(out, *value);
        }
    }
}

/// Appends an integer in decimal.
fn push_integer(out: &mut String, value: i32) {
    write!(out, "{value}").unwrap();
}

/// Appends a float, widened to a double, as htslib's `kputd` writes it: zero as `0` or `-0`, a
/// magnitude from 0.0001 to 999999 with 6 significant digits, rounded half to even and without
/// trailing zeros, and any other value as C's `%g` writes it.
fn push_float(out: &mut String, value: f32) {
    let mut value = f64::from(value);
    if value == 0.0 {
        out.push_str(if value.is_sign_negative() { "-0" } else { "0" });
        return;
    }
    if value < 0.0 {
        out.push('-');
        value = -value;
    }
    if !(0.0001..=999999.0).contains(&value) {
        push_general(out, value);
        return;
    }

    // As in `kputd`: the 6 digits, or 7 after rounding up, are written right to left, ending `shift`
    // places before the end of the buffer, and the point is then put among them.
    let (scale, shift) = match value {
        v if v < 0.001 => (1e9, 1),
        v if v < 0.01 => (1e8, 2),
        v if v < 0.1 => (1e7, 3),
        v if v < 1.0 => (1e6, 4),
        v if v < 10.0 => (1e5, 5),
        v if v < 100.0 => (1e4, 6),
        v if v < 1000.0 => (1e3, 7),
        v if v < 10000.0 => (1e2, 8),
        v if v < 100000.0 => (1e1, 9),
        _ => (1.0, 10),
    };
    let mut digits = (value * scale).round_ties_even() as u32;
    let mut buffer = [0u8; 21];
    let mut start = 20 - shift;
    for _ in 0..3 {
        start -= 2;
        let pair = (digits % 100) as u8;
        buffer[start] = b'0' + pair / 10;
        buffer[start + 1] = b'0' + pair % 10;
        digits /= 100;
    }
    if digits >= 1 {
        start -= 1;
        buffer[start] = b'0' + digits as u8;
    }

    let mut length = 20 - start;
    let mut end;
    if length <= 10 {
        end = start + 5;
        while length < 10 {
            start -= 1;
            buffer[start] = b'0';
            length += 1;
        }
        start -= 2;
        buffer[start] = b'0';
        buffer[start + 1] = b'.';
    } else {
        start -= 1;
        end = start + 6;
        let mut point = start;
        while length > 10 {
            buffer[point] = buffer[point + 1];
            point += 1;
            length -= 1;
        }
        buffer[point] = b'.';
    }
    while buffer[end] == b'0' && end > start {
        end -= 1;
    }
    if buffer[end] != 0 && buffer[end] != b'.' {
        end += 1;
    }
    out.push_str(std::str::from_utf8(&buffer[start..end]).unwrap());
}

/// Appends a positive value as C's `%g` writes it: 6 significant digits, rounded half to even,
/// without trailing zeros, in exponent notation when its exponent is below -4 or above 5; a NaN is
/// `nan`, or `-nan` with its sign bit set, as glibc writes it.
fn push_general(out: &mut String, value: f64) {
    if value.is_nan() {
        out.push_str(if value.is_sign_negative() {
            "-nan"
        } else {
            "nan"
        });
        return;
    }
    if value.is_infinite() {
        out.push_str("inf");
        return;
    }
    let scientific = format!("{value:.5e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    if (-4..6).contains(&exponent) {
        let fixed = format!("{value:.*}", (5 - exponent) as usize);
        out.push_str(trim_fraction(&fixed));
    } else {
        out.push_str(trim_fraction(mantissa));
        let sign = if exponent < 0 { '-' } else { '+' };
        write!(out, "e{sign}{:02}", exponent.abs()).unwrap();
    }
}

/// Removes the trailing zeros of a number's fraction, and its point when no fraction is left.
fn trim_fraction(number: &str) -> &str {
    if number.contains('.') {
        number.trim_end_matches('0').trim_end_matches('.')
    } else {
        number
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    /// Floats as htslib 1.19's `kputd` writes them, by their bits: chosen edge cases, then values
    /// drawn from every bit pattern but NaN.
    const HTSLIB_FLOATS: &[(u32, &str)] = &[
        (0x00000000, "0"),
        (0x80000000, "-0"),
        (0x3f800000, "1"),
        (0xbf800000, "-1"),
        (0x3f000000, "0.5"),
        (0x3f2d0785, "0.675896"),
        (0x2fa36bc4, "2.97261e-10"),
        (0x38d1b717, "0.0001"),
        (0xb8d1b717, "-0.0001"),
        (0x3a83126f, "0.001"),
        (0x3c23d70a, "0.01"),
        (0x3dcccccd, "0.1"),
        (0x41200000, "10"),
        (0x42c80000, "100"),
        (0x447a0000, "1000"),
        (0x461c4000, "10000"),
        (0x47c35000, "100000"),
        (0x497423f0, "999999"),
        (0x497423f6, "999999"),
        (0x497423f8, "1e+06"),
        (0x49742400, "1e+06"),
        (0x4996b438, "1.23457e+06"),
        (0x4b800000, "1.67772e+07"),
        (0x7f7fffff, "3.40282e+38"),
        (0xff7fffff, "-3.40282e+38"),
        (0x00800000, "1.17549e-38"),
        (0x34000000, "1.19209e-07"),
        (0x00000001, "1.4013e-45"),
        (0x80000001, "-1.4013e-45"),
        (0x38d1b716, "0.0001"),
        (0x38d1b68e, "9.9999e-05"),
        (0x38d1b719, "0.0001"),
        (0x3a83126e, "0.001"),
        (0x3c23d70a, "0.01"),
        (0x3dcccccc, "0.1"),
        (0x3f7fffff, "1"),
        (0x411fffff, "10"),
        (0x42c7ffff, "100"),
        (0x4479fffe, "1000"),
        (0x461c3fff, "10000"),
        (0x47c34fff, "100000"),
        (0x3901742d, "0.000123457"),
        (0x3aa1d139, "0.00123457"),
        (0x3c4a4587, "0.0123457"),
        (0x3dfcd6e9, "0.123457"),
        (0x3f9e0651, "1.23457"),
        (0x414587e6, "12.3457"),
        (0x42f6e9df, "123.457"),
        (0x449a522b, "1234.57"),
        (0x4640e6b6, "12345.7"),
        (0x47f12064, "123457"),
        (0x3eaaaaab, "0.333333"),
        (0x3f2aaaab, "0.666667"),
        (0x40490fdb, "3.14159"),
        (0x402df854, "2.71828"),
        (0x3d4ccccd, "0.05"),
        (0x3e4ccccd, "0.2"),
        (0x3e800000, "0.25"),
        (0x3eb33333, "0.35"),
        (0x3f4ccccd, "0.8"),
        (0x3727c5ac, "1e-05"),
        (0x3851b717, "5e-05"),
        (0x377ba882, "1.5e-05"),
        (0x34bbe7a2, "3.5e-07"),
        (0x0da24260, "1e-30"),
        (0x7149f2ca, "1e+30"),
        (0x000116c2, "9.99995e-41"),
        (0x39d0eeb4, "0.000398507"),
        (0x39aeae65, "0.000333178"),
        (0x3b1a8ee2, "0.00235837"),
        (0x3643b1b3, "2.91607e-06"),
        (0x3ceb61ae, "0.0287331"),
        (0x3a3663aa, "0.00069576"),
        (0x3f7fffef, "0.999999"),
        (0x3f7ffb48, "0.999928"),
        (0x3f1a0ac8, "0.601727"),
        (0x42720000, "60.5"),
        (0x42053333, "33.3"),
        (0x399d4952, "0.0003"),
        (0x3951b717, "0.0002"),
        (0x45000000, "2048"),
        (0x47800040, "65536.5"),
        (0x40900000, "4.5"),
        (0x3e000000, "0.125"),
        (0x7f800000, "inf"),
        (0xff800000, "-inf"),
        (0x7fc00000, "nan"),
        (0xe124b63a, "-1.899e+20"),
        (0x8b9a74ab, "-5.94942e-32"),
        (0x64e1b3ac, "3.33077e+22"),
        (0x00174626, "2.13738e-39"),
        (0xf2adbbaf, "-6.88228e+30"),
        (0xfed75123, "-1.43103e+38"),
        (0x8a94501a, "-1.4282e-32"),
        (0x12751a71, "7.7341e-28"),
        (0x96573f6c, "-1.73876e-25"),
        (0x46ea7191, "30008.8"),
        (0x13d2ea5d, "5.32426e-27"),
        (0x9db4cf31, "-4.78598e-21"),
        (0x8e0f4e18, "-1.76637e-30"),
        (0x9e43c23e, "-1.03634e-20"),
        (0x268a56bc, "9.59918e-16"),
        (0xe7e1f2d2, "-2.13402e+24"),
        (0xeec01fef, "-2.97299e+28"),
        (0x4a8ca751, "4.60894e+06"),
        (0x12bbe422, "1.18576e-27"),
        (0xa9cdf49d, "-9.14626e-14"),
        (0xfc95b972, "-6.21931e+36"),
        (0x3cc0494f, "0.0234725"),
        (0x88dfc4db, "-1.34676e-33"),
        (0x78d703d9, "3.48882e+34"),
        (0x8d219e6f, "-4.98027e-31"),
        (0x63680239, "4.2798e+21"),
        (0x06cd666e, "7.72629e-35"),
        (0xea1e9eae, "-4.79399e+25"),
        (0x00a30b2b, "1.49732e-38"),
        (0x590d22c8, "2.48289e+15"),
        (0x57dfd022, "4.9217e+14"),
        (0x16a31f2f, "2.63537e-25"),
        (0xdd9e740c, "-1.42722e+18"),
        (0x70e04de3, "5.5535e+29"),
        (0x52de38ed, "4.77219e+11"),
        (0x2db9938c, "2.10976e-11"),
        (0xe6cb9168, "-4.80661e+23"),
        (0x083db87b, "5.7092e-34"),
        (0x59627ba2, "3.98433e+15"),
        (0xd4d02589, "-7.15186e+12"),
        (0xdc4cda99, "-2.30645e+17"),
        (0xa4e4fbd6, "-9.93059e-17"),
        (0x485ae539, "224149"),
        (0x8b4427a7, "-3.77781e-32"),
        (0xf9a8cf9f, "-1.09565e+35"),
        (0xeb30a9f2, "-2.13573e+26"),
        (0x3fdc4855, "1.72096"),
        (0x6c00d4fe, "6.22993e+26"),
        (0xa57ad991, "-2.17578e-16"),
        (0x37585015, "1.28932e-05"),
        (0x960739b8, "-1.09234e-25"),
        (0x57302520, "1.93673e+14"),
        (0x211591aa, "5.06759e-19"),
        (0xf7339f7a, "-3.64319e+33"),
        (0x1f4f3f94, "4.38865e-20"),
        (0xef05ba8a, "-4.1387e+28"),
        (0x52ce02a0, "4.42404e+11"),
        (0xc1d3364d, "-26.4015"),
        (0x44427dc0, "777.965"),
        (0x74b57f9d, "1.15038e+32"),
        (0xb390f5fe, "-6.75027e-08"),
        (0x08c30e49, "1.17395e-33"),
        (0x4849434c, "206093"),
        (0x643e98dc, "1.40636e+22"),
        (0x538d2a8e, "1.21261e+12"),
        (0x2d4eade0, "1.17484e-11"),
        (0xe6a8e2b9, "-3.9877e+23"),
        (0xa5084706, "-1.18202e-16"),
        (0x10f2efb2, "9.58215e-29"),
        (0xed95af30, "-5.79063e+27"),
        (0x5603e229, "3.62518e+13"),
        (0x629c364a, "1.4408e+21"),
        (0x6ef58860, "3.79943e+28"),
        (0x20c5141c, "3.33864e-19"),
        (0xca9c72df, "-5.12651e+06"),
        (0xdc31a73c, "-2.0002e+17"),
        (0xf21c39b7, "-3.09436e+30"),
        (0xd0768762, "-1.65443e+10"),
        (0x13c222cf, "4.90068e-27"),
        (0xa4e6c942, "-1.00088e-16"),
        (0xc4184305, "-609.047"),
        (0x43682219, "232.133"),
        (0xa24f100c, "-2.80622e-18"),
        (0x4998b54b, "1.25099e+06"),
        (0xb90ea0b3, "-0.00013602"),
        (0xce0631df, "-5.62854e+08"),
        (0x0f876de1, "1.33544e-29"),
        (0xa55ca37c, "-1.91373e-16"),
        (0x17544745, "6.85908e-25"),
        (0x6829bbfb, "3.20619e+24"),
        (0xb5887e50, "-1.01696e-06"),
        (0xf2064d51, "-2.66013e+30"),
        (0x4e226067, "6.81057e+08"),
        (0x47feaf70, "130399"),
        (0xd00c2978, "-9.40611e+09"),
        (0xf1437ec9, "-9.68046e+29"),
        (0x4dd82104, "4.53255e+08"),
        (0x76e83af8, "2.3551e+33"),
        (0x47574643, "55110.3"),
        (0x5c71400c, "2.71624e+17"),
        (0xfa6fbcb4, "-3.11197e+35"),
        (0xb2de7348, "-2.58966e-08"),
        (0xea5eef73, "-6.73781e+25"),
        (0xc1a201cb, "-20.2509"),
        (0xb2ff01c6, "-2.96867e-08"),
        (0x0a3afc05, "9.00297e-33"),
        (0xe2f4add8, "-2.25677e+21"),
        (0x9ebd599f, "-2.00482e-20"),
        (0x845ac858, "-2.57178e-36"),
        (0x776578f0, "4.65425e+33"),
        (0xd7198d6d, "-1.68833e+14"),
        (0x303f98d7, "6.97026e-10"),
        (0xa78631e5, "-3.72466e-15"),
    ];

    fn float(value: f32) -> String {
        let mut out = String::new();
        push_float(&mut out, value);
        out
    }

    #[test]
    fn test_push_float_writes_floats_as_htslib() {
        for &(bits, expected) in HTSLIB_FLOATS {
            assert_eq!(float(f32::from_bits(bits)), expected, "{bits:#010x}");
        }
    }

    #[rstest]
    #[case(1000005.0, "1e+06")]
    #[case(-1000005.0, "-1e+06")]
    #[case(1000105.0, "1.0001e+06")]
    #[case(1000025.0, "1.00002e+06")]
    #[case(10000050.0, "1e+07")]
    #[case(f32::from_bits(0x7fc0_0000), "nan")]
    #[case(f32::from_bits(0xffc0_0000), "-nan")]
    fn test_push_float_writes_ties_and_nans_as_htslib_on_glibc(
        #[case] value: f32,
        #[case] expected: &str,
    ) {
        assert_eq!(float(value), expected);
    }

    #[test]
    fn test_missing_values_are_sentinels() {
        assert!(i32::missing().is_missing());
        assert!(f32::missing().is_missing());
        assert!(!f32::NAN.is_missing());
        assert!(!0.is_missing());
        assert!(!0.0f32.is_missing());
    }

    /// A header of two samples declaring the fields the records below use.
    fn header() -> Header {
        let mut header = Header::default();
        header.push_sample("S1").push_sample("S2");
        for line in [
            "##source=probe",
            r#"##INFO=<ID=S,Number=1,Type=String,Description="s">"#,
            r#"##INFO=<ID=I,Number=A,Type=Integer,Description="i">"#,
            r#"##INFO=<ID=I,Number=1,Type=Float,Description="dup">"#,
            r#"##INFO=<ID=F,Number=.,Type=Float,Description="f">"#,
            r#"##FILTER=<ID=PASS,Description="pass">"#,
            r#"##FILTER=<ID=A,Description="a">"#,
            r#"##FILTER=<ID=B,Description="b">"#,
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">"#,
            r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">"#,
            r#"##FORMAT=<ID=X,Number=1,Type=Float,Description="x">"#,
            r#"##FORMAT=<ID=E,Number=1,Type=Integer,Description="e">"#,
            "##contig=<ID=c1,length=10>",
            "##contig=<ID=c2,length=20>",
            "##contig=<ID=c1,length=30>",
            "##source=probe",
            "##reference=x",
        ] {
            header.push_record(line);
        }
        header
    }

    #[test]
    fn test_header_keeps_the_order_lines_were_pushed_in_and_drops_repeats() {
        let expected = [
            "##fileformat=VCFv4.2",
            "##source=probe",
            r#"##INFO=<ID=S,Number=1,Type=String,Description="s">"#,
            r#"##INFO=<ID=I,Number=A,Type=Integer,Description="i">"#,
            r#"##INFO=<ID=F,Number=.,Type=Float,Description="f">"#,
            r#"##FILTER=<ID=PASS,Description="pass">"#,
            r#"##FILTER=<ID=A,Description="a">"#,
            r#"##FILTER=<ID=B,Description="b">"#,
            r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">"#,
            r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">"#,
            r#"##FORMAT=<ID=X,Number=1,Type=Float,Description="x">"#,
            r#"##FORMAT=<ID=E,Number=1,Type=Integer,Description="e">"#,
            "##contig=<ID=c1,length=10>",
            "##contig=<ID=c2,length=20>",
            "##reference=x",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        ];
        let header = header();
        assert_eq!(
            String::from_utf8(header.to_bytes()).unwrap(),
            format!("{}\n", expected.join("\n"))
        );
        assert_eq!(header.contig_id("c1"), Ok(0));
        assert_eq!(header.contig_id("c2"), Ok(1));
        assert_eq!(
            header.contig_id("c9"),
            Err("Contig 'c9' has no ##contig line in the VCF header!".to_string())
        );
    }

    type Fill = fn(&mut Record) -> Result<(), String>;

    #[rstest]
    #[case::pass_gives_way(|r: &mut Record| { r.push_filter("PASS")?; r.push_filter("A") }, "A\t.")]
    #[case::pass_replaces(|r: &mut Record| { r.push_filter("A")?; r.push_filter("PASS") }, "PASS\t.")]
    #[case::label_once(|r: &mut Record| { r.push_filter("A")?; r.push_filter("B")?; r.push_filter("A") }, "A;B\t.")]
    #[case::empty_info(|r: &mut Record| { r.push_info_string("S", "")?; r.push_info_integer("I", &[]) }, ".\t.")]
    #[case::missing_info(
        |r: &mut Record| {
            r.push_info_integer("I", &[i32::missing(), 3])?;
            r.push_info_float("F", &[f32::missing()])?;
            r.push_info_string("S", "x")
        },
        ".\tI=.,3;F=.;S=x"
    )]
    #[case::empty_format(
        |r: &mut Record| { r.push_format_integer("AD", &[])?; r.push_format_integer("E", &[i32::missing(), -5]) },
        ".\t.\tE\t.\t-5"
    )]
    #[case::per_sample(
        |r: &mut Record| {
            r.push_genotypes(&[
                GenotypeAllele::Unphased(0),
                GenotypeAllele::Unphased(1),
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing,
            ])?;
            r.push_format_integer("AD", &[i32::missing(), 4, 7, i32::missing()])?;
            r.push_format_float("X", &[0.675896, f32::missing()])
        },
        ".\t.\tGT:AD:X\t0/1:.,4:0.675896\t./.:7,.:."
    )]
    #[case::nothing(|_: &mut Record| Ok(()), ".\t.")]
    fn test_record_writes_fields_as_htslib(#[case] fill: Fill, #[case] expected: &str) {
        let header = header();
        let mut record = Record::new(&header, "c2", 5, "A", "C").unwrap();
        fill(&mut record).unwrap();
        assert_eq!(record.rid(), 1);
        assert_eq!(
            record.to_vcf_string(),
            format!("c2\t5\t.\tA\tC\t.\t{expected}\n")
        );
    }

    #[test]
    fn test_record_refuses_what_the_header_does_not_declare() {
        let header = header();
        assert_eq!(
            Record::new(&header, "c9", 5, "A", "C").err(),
            Some("Contig 'c9' has no ##contig line in the VCF header!".to_string())
        );
        let mut record = Record::new(&header, "c1", 5, "A", "C").unwrap();
        assert_eq!(
            record.push_info_integer("Z", &[1]),
            Err("The INFO field 'Z' is not declared in the VCF header!".to_string())
        );
        assert_eq!(
            record.push_format_integer("Z", &[1, 2]),
            Err("The FORMAT field 'Z' is not declared in the VCF header!".to_string())
        );
        assert_eq!(
            record.push_filter("Z"),
            Err("The FILTER field 'Z' is not declared in the VCF header!".to_string())
        );
        assert_eq!(
            record.push_format_integer("E", &[1, 2, 3]),
            Err("Expected a multiple of 2 values for the FORMAT field 'E', found 3!".to_string())
        );
    }
}
