//! VCF records built from vartovcf's values, whose missing ones are sentinels, as noodles-vcf
//! records, each field checked against the header, and written with noodles-vcf.
use std::error::Error;

use noodles::core::Position;
use noodles::vcf::variant::io::Write as _;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
use noodles::vcf::variant::record_buf::samples::Keys;
use noodles::vcf::variant::record_buf::samples::sample::Value as FormatValue;
use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
use noodles::vcf::variant::record_buf::samples::sample::value::{Array as FormatArray, Genotype};
use noodles::vcf::variant::record_buf::{AlternateBases, Samples};
use noodles::vcf::{self, Header, variant::RecordBuf};

/// The bits of the float that stands for a missing value, a NaN no arithmetic produces.
const MISSING_FLOAT_BITS: u32 = 0x7F80_0001;

/// An integer or float value with a sentinel that stands for a missing value.
pub trait Numeric: Copy {
    /// The value that stands for a missing one.
    fn missing() -> Self;

    /// Whether this is the value that stands for a missing one.
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

/// The header as VCF text, ending with the line of column names.
pub fn header_bytes(header: &Header) -> Result<Vec<u8>, String> {
    let mut writer = vcf::io::Writer::new(Vec::new());
    writer
        .write_header(header)
        .map_err(|error| format!("Could not write the VCF header: {error}!"))?;
    Ok(writer.into_inner())
}

/// One VCF record, built field by field and written with noodles-vcf; ID and QUAL are always
/// missing, and a field the header does not declare is refused.
pub struct Record<'h> {
    header: &'h Header,
    rid: u32,
    record: RecordBuf,
    keys: Keys,
    samples: Vec<Vec<Option<FormatValue>>>,
}

impl<'h> Record<'h> {
    /// Starts a record of a contig the header declares, at a 1-based position, with its REF and an
    /// ALT that is `.` when there is none. A REF base other than A, C, G, T or N is refused, as VCF
    /// has no place for it and noodles-vcf would write an ambiguity code as one of its bases.
    pub fn new(
        header: &'h Header,
        contig: &str,
        position: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Result<Self, String> {
        let rid = header
            .contigs()
            .get_index_of(contig)
            .and_then(|rid| u32::try_from(rid).ok())
            .ok_or_else(|| format!("Contig '{contig}' has no ##contig line in the VCF header!"))?;
        let start = usize::try_from(position)
            .ok()
            .and_then(|position| Position::try_from(position).ok())
            .ok_or_else(|| {
                format!("Expected a position of at least 1, found {contig}:{position}!")
            })?;
        if let Some(base) = ref_allele
            .bytes()
            .find(|base| !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
        {
            return Err(format!(
                "Expected REF bases A, C, G, T or N at {contig}:{position}, found '{}' in {ref_allele}!",
                char::from(base)
            ));
        }
        let alternate_bases = match alt_allele {
            "." => Vec::new(),
            alt_allele => vec![alt_allele.to_string()],
        };
        let record = RecordBuf::builder()
            .set_reference_sequence_name(contig)
            .set_variant_start(start)
            .set_reference_bases(ref_allele)
            .set_alternate_bases(AlternateBases::from(alternate_bases))
            .build();
        Ok(Record {
            header,
            rid,
            record,
            keys: Keys::default(),
            samples: vec![Vec::new(); header.sample_names().len()],
        })
    }

    /// The index of the record's contig among the header's `##contig` lines.
    pub fn rid(&self) -> u32 {
        self.rid
    }

    /// Adds a FILTER label. `PASS` replaces every label added before it and gives way to any label
    /// added after it, and a label added twice is written once.
    pub fn push_filter(&mut self, label: &str) -> Result<(), String> {
        declared("FILTER", self.header.filters().contains_key(label), label)?;
        let filters = self.record.filters_mut().as_mut();
        if label == "PASS" || (filters.len() == 1 && filters.contains("PASS")) {
            filters.clear();
        }
        filters.insert(label.to_string());
        Ok(())
    }

    /// Adds a String INFO field; an empty value leaves the field out.
    pub fn push_info_string(&mut self, key: &str, value: &str) -> Result<(), String> {
        self.check_info(key)?;
        if !value.is_empty() {
            let value = InfoValue::String(value.to_string());
            self.record.info_mut().insert(key.to_string(), Some(value));
        }
        Ok(())
    }

    /// Adds an Integer INFO field; no values leave the field out.
    pub fn push_info_integer(&mut self, key: &str, values: &[i32]) -> Result<(), String> {
        self.push_info(key, values, InfoValue::Integer, InfoArray::Integer)
    }

    /// Adds a Float INFO field; no values leave the field out.
    pub fn push_info_float(&mut self, key: &str, values: &[f32]) -> Result<(), String> {
        self.push_info(key, values, InfoValue::Float, InfoArray::Float)
    }

    /// Adds an INFO field of one value, which may be missing, or of several.
    fn push_info<V: Numeric>(
        &mut self,
        key: &str,
        values: &[V],
        scalar: fn(V) -> InfoValue,
        array: fn(Vec<Option<V>>) -> InfoArray,
    ) -> Result<(), String> {
        self.check_info(key)?;
        let value = match values {
            [] => return Ok(()),
            [value] => present(*value).map(scalar),
            values => Some(InfoValue::Array(array(options(values)))),
        };
        self.record.info_mut().insert(key.to_string(), value);
        Ok(())
    }

    /// Adds the GT FORMAT field from the alleles of each sample in turn; no alleles leave it out.
    pub fn push_genotypes(&mut self, alleles: &[GenotypeAllele]) -> Result<(), String> {
        self.push_format("GT", alleles, |alleles| {
            let genotype: Genotype = alleles
                .iter()
                .map(|allele| {
                    let position = match allele {
                        GenotypeAllele::Unphased(index) => usize::try_from(*index).ok(),
                        GenotypeAllele::UnphasedMissing => None,
                    };
                    Allele::new(position, Phasing::Unphased)
                })
                .collect();
            Some(FormatValue::Genotype(genotype))
        })
    }

    /// Adds an Integer FORMAT field from the values of each sample in turn; no values leave it out.
    pub fn push_format_integer(&mut self, key: &str, values: &[i32]) -> Result<(), String> {
        self.push_format(key, values, |values| match values {
            [value] => present(*value).map(FormatValue::Integer),
            values => Some(FormatValue::Array(FormatArray::Integer(options(values)))),
        })
    }

    /// Adds a Float FORMAT field from the values of each sample in turn; no values leave it out.
    pub fn push_format_float(&mut self, key: &str, values: &[f32]) -> Result<(), String> {
        self.push_format(key, values, |values| match values {
            [value] => present(*value).map(FormatValue::Float),
            values => Some(FormatValue::Array(FormatArray::Float(options(values)))),
        })
    }

    /// Adds a FORMAT field, with an equal share of the values for each sample.
    fn push_format<V>(
        &mut self,
        key: &str,
        values: &[V],
        value: impl Fn(&[V]) -> Option<FormatValue>,
    ) -> Result<(), String> {
        declared("FORMAT", self.header.formats().contains_key(key), key)?;
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
        self.keys.as_mut().insert(key.to_string());
        let width = values.len() / samples;
        for (sample, values) in self.samples.iter_mut().zip(values.chunks_exact(width)) {
            sample.push(value(values));
        }
        Ok(())
    }

    /// Refuses an INFO field the header does not declare.
    fn check_info(&self, key: &str) -> Result<(), String> {
        declared("INFO", self.header.infos().contains_key(key), key)
    }

    /// Appends the record to `line` as one line of VCF text, ending with a newline; without a FORMAT
    /// field it has no FORMAT or sample columns.
    pub fn write(mut self, line: &mut Vec<u8>) -> Result<(), String> {
        if !self.keys.as_ref().is_empty() {
            *self.record.samples_mut() = Samples::new(self.keys, self.samples);
        }
        vcf::io::Writer::new(line)
            .write_variant_record(self.header, &self.record)
            .map_err(|error| {
                let mut cause: &dyn Error = &error;
                while let Some(source) = cause.source() {
                    cause = source;
                }
                format!(
                    "Could not write the record at {}:{}: {cause}!",
                    self.record.reference_sequence_name(),
                    self.record.variant_start().map_or(0, usize::from)
                )
            })
    }
}

/// Refuses a field the header does not declare.
fn declared(kind: &str, declared: bool, id: &str) -> Result<(), String> {
    if declared {
        Ok(())
    } else {
        Err(format!(
            "The {kind} field '{id}' is not declared in the VCF header!"
        ))
    }
}

/// A value, or `None` for the one that stands for a missing value.
fn present<V: Numeric>(value: V) -> Option<V> {
    (!value.is_missing()).then_some(value)
}

/// Each value, or `None` for the one that stands for a missing value.
fn options<V: Numeric>(values: &[V]) -> Vec<Option<V>> {
    values.iter().copied().map(present).collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

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
        [
            "##fileformat=VCFv4.2",
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
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
        ]
        .join("\n")
        .parse()
        .unwrap()
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
    fn test_record_writes_its_fields(#[case] fill: Fill, #[case] expected: &str) {
        let header = header();
        let mut record = Record::new(&header, "c2", 5, "A", "C").unwrap();
        fill(&mut record).unwrap();
        assert_eq!(record.rid(), 1);
        let mut line = Vec::new();
        record.write(&mut line).unwrap();
        assert_eq!(
            String::from_utf8(line).unwrap(),
            format!("c2\t5\t.\tA\tC\t.\t{expected}\n")
        );
    }

    #[test]
    fn test_record_writes_no_alt_allele_as_a_dot() {
        let header = header();
        let mut line = Vec::new();
        Record::new(&header, "c1", 3, "G", ".")
            .unwrap()
            .write(&mut line)
            .unwrap();
        assert_eq!(
            String::from_utf8(line).unwrap(),
            "c1\t3\t.\tG\t.\t.\t.\t.\n"
        );
    }

    #[rstest]
    #[case("R", 'R')]
    #[case("AMA", 'M')]
    #[case("y", 'y')]
    #[case("X", 'X')]
    fn test_record_refuses_a_ref_base_vcf_has_no_place_for(
        #[case] ref_allele: &str,
        #[case] base: char,
    ) {
        let header = header();
        assert!(Record::new(&header, "c1", 3, "ACGTNacgtn", "A").is_ok());
        assert_eq!(
            Record::new(&header, "c1", 3, ref_allele, "A").err(),
            Some(format!(
                "Expected REF bases A, C, G, T or N at c1:3, found '{base}' in {ref_allele}!"
            ))
        );
    }

    #[test]
    fn test_record_names_why_noodles_refuses_to_write_it() {
        let header: Header = "##fileformat=VCFv4.2\n##contig=<ID=c(1),length=10>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            .parse()
            .unwrap();
        let mut line = Vec::new();
        let error = Record::new(&header, "c(1)", 3, "A", "C")
            .unwrap()
            .write(&mut line)
            .unwrap_err();
        assert_eq!(
            error,
            "Could not write the record at c(1):3: invalid reference sequence name: c(1)!"
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
