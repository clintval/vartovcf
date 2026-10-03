//! A module for serialization-deserialization friendly VarDict/VarDictJava data types.
use std::clone::Clone;
use std::cmp::PartialEq;
use std::default::Default;
use std::error;
use std::fmt;
use std::ops::Range;
use std::str::FromStr;

use anyhow::Result;
use bio_types::genome::{AbstractInterval, Position};
use rust_htslib::bcf::Header;
use rust_htslib::bcf::record::{GenotypeAllele, Numeric};
use serde::{Deserialize, Serialize, de::Error};
use serde_with::DisplayFromStr;
use serde_with::serde_as;
use strum::EnumString;

const CARGO_PKG_NAME: &str = env!("CARGO_PKG_NAME");
const CARGO_PKG_VERSION: &str = env!("CARGO_PKG_VERSION");

/// The minimum allele frequency at which a variant call is genotyped homozygous alternate.
pub const MIN_HOM_ALT_AF: f32 = 0.8;

/// Deserialize a possibly infinite float into a <f32> or return a custom error. The floating point
/// number in VAR files must be expressed as a ratio for them to be true odds ratios.
///
/// Correctly serializes the following cases:
///
/// * `"Inf"`: floating point infinity with title-case string, will return `0`.
/// * `(0, 1)`: a floating point number in the range 0 to 1 (exclusive), will return `1 / n`
/// * All other values (including 0, 1, and values > 1): returned as-is
fn maybe_infinite_f32_odds_ratio<'de, D>(deserializer: D) -> Result<f32, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let string: &str = Deserialize::deserialize(deserializer)?;
    f32::from_str(string.to_lowercase().as_str())
        .map(|x| {
            if x.is_infinite() && x.is_sign_positive() {
                0.0
            } else if x > 0.0 && x < 1.0 {
                1.0 / x
            } else {
                x
            }
        })
        .map_err(D::Error::custom)
}

/// Deserialize encoded structural variant (SV) info or return a custom error.
///
/// The following cases are handled:
///
/// * `0`: `None`
/// * `#-#-#`: `Some(SvInfo)`
fn maybe_sv_info<'de, D>(deserializer: D) -> Result<Option<SvInfo>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let string: &str = Deserialize::deserialize(deserializer)?;
    match string {
        "0" => Ok(None),
        _ => Ok(Some(SvInfo::from_str(string).map_err(D::Error::custom)?)),
    }
}

/// Deserialize a duplication rate that may not exist or return a custom error.
///
/// The following cases are handled:
///
/// * `0`: `None`
/// * `#`: `Some(#)`
fn maybe_duplication_rate<'de, D>(deserializer: D) -> Result<Option<f32>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let string: &str = Deserialize::deserialize(deserializer)?;
    match string {
        "0" => Ok(None),
        _ => Ok(Some(f32::from_str(string).map_err(D::Error::custom)?)),
    }
}

/// An exception for when we cannot parse a string into a `PairBias`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParsePairBiasError;

impl error::Error for ParsePairBiasError {}

impl fmt::Display for ParsePairBiasError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

/// An exception for when we cannot parse a string into a `SvInfo`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseSvInfoError;

impl error::Error for ParseSvInfoError {}

impl fmt::Display for ParseSvInfoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

/// Enumeration of VarDict/VarDictJava strand bias statuses.
#[derive(Debug, Deserialize, EnumString, Eq, PartialEq, Serialize)]
pub enum StrandBias {
    /// There were 12 or fewer reads, all on one strand.
    #[strum(to_string = "0")]
    TooFewReads,
    /// Strand bias was detected.
    #[strum(to_string = "1")]
    Detected,
    /// Strand bias was undetected.
    #[strum(to_string = "2")]
    UnDetected,
}

impl fmt::Display for StrandBias {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

/// The strand bias status for a reference allele alternate allele pair.
#[derive(Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct PairBias {
    /// The reference allele strand bias status.
    pub reference: StrandBias,
    /// The alternate allele strand bias status.
    pub alternate: StrandBias,
}

impl Default for PairBias {
    /// The default paired bias status is no strand bias detected.
    fn default() -> Self {
        PairBias {
            reference: StrandBias::UnDetected,
            alternate: StrandBias::UnDetected,
        }
    }
}

impl fmt::Display for PairBias {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}", self.reference, self.alternate)
    }
}

impl FromStr for PairBias {
    type Err = ParsePairBiasError;

    /// Convert a string to a `StrandBias` status.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let items: Vec<&str> = s.split(';').collect();
        if items.len() != 2 {
            return Err(ParsePairBiasError);
        }
        let pair = match (items.first(), items.get(1)) {
            (Some(reference), Some(alternate)) => PairBias {
                reference: StrandBias::from_str(reference).map_err(|_| ParsePairBiasError)?,
                alternate: StrandBias::from_str(alternate).map_err(|_| ParsePairBiasError)?,
            },
            (_, _) => return Err(ParsePairBiasError),
        };
        Ok(pair)
    }
}

/// A container for structural variant (SV) information.
#[derive(Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct SvInfo {
    /// The number of split reads supporting the SV.
    pub supporting_split_reads: i32,
    /// The number of read pairs supporting the SV.
    pub supporting_pairs: i32,
    /// The number of clusters supporting the SV.
    pub supporting_clusters: i32,
}

impl Default for SvInfo {
    /// The default has all fields set to zero.
    fn default() -> Self {
        SvInfo {
            supporting_split_reads: 0,
            supporting_pairs: 0,
            supporting_clusters: 0,
        }
    }
}

impl fmt::Display for SvInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}-{}-{}",
            self.supporting_split_reads, self.supporting_pairs, self.supporting_clusters
        )
    }
}

impl FromStr for SvInfo {
    type Err = ParseSvInfoError;

    /// Convert a string to a `SvInfo`.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let items: Vec<&str> = s.split('-').collect();
        let sv_info = match (items.first(), items.get(1), items.get(2)) {
            (Some(split_reads), Some(pairs), Some(clusters)) => SvInfo {
                supporting_split_reads: split_reads.parse().map_err(|_| ParseSvInfoError)?,
                supporting_pairs: pairs.parse().map_err(|_| ParseSvInfoError)?,
                supporting_clusters: clusters.parse().map_err(|_| ParseSvInfoError)?,
            },
            (_, _, _) => return Err(ParseSvInfoError),
        };
        Ok(sv_info)
    }
}

/// A record of output from VarDict/VarDictJava run in tumor-only mode.
#[serde_as]
#[derive(Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct TumorOnlyVariant<'a> {
    /// The sample name as VarDict prints it: its -N value, or else one taken from the BAM file name.
    pub sample: &'a str,
    /// The name of the interval this variant call overlaps.
    pub interval_name: &'a str,
    /// The reference sequence name.
    pub contig: &'a str,
    /// The 1-based start of this variant call.
    pub start: u64,
    /// The 1-based inclusive end of this variant call.
    pub end: u64,
    /// The reference simple allele.
    pub ref_allele: &'a str,
    /// The alternate allele, simple or symbolic.
    pub alt_allele: &'a str,
    /// The total allele depth at this call locus.
    pub depth: i32,
    /// The total alternate depth at this call locus.
    pub alt_depth: i32,
    /// The number of forward reads supporting the reference call.
    pub ref_forward: i32,
    /// The number of reverse reads supporting the reference call.
    pub ref_reverse: i32,
    /// The number of forward reads supporting the alternate call.
    pub alt_forward: i32,
    /// The number of reverse reads supporting the alternate call.
    pub alt_reverse: i32,
    /// VarDict's genotype column, not a VCF genotype: the reference allele (or the leading ALT allele
    /// when the reference is below -f), a slash, then this call's allele, in VarDict's notation.
    pub gt: &'a str,
    /// The allele frequency of the alternate allele.
    pub af: f32,
    /// Strand bias status. That will take the values [0-2];[0-2] (_e.g._ "0;2", "2;1"). The first
    /// value refers to reads that support the reference allele, and the second to reads that
    /// support the variant allele.
    ///
    /// * `0`: there were 12 or fewer reads, all on one strand
    /// * `1`: strand bias was detected
    /// * `2`: strand bias was undetected
    #[serde_as(as = "DisplayFromStr")]
    pub strand_bias: PairBias,
    /// The mean 1-based distance from the variant to the nearer end of the aligned part of each read
    /// that supports the variant call.
    pub mean_position_in_read: f32,
    /// 1 when the supporting reads place the variant at two or more distances from the read end,
    /// otherwise 0; VarDict prints a flag here, not a standard deviation.
    pub stdev_position_in_read: f32,
    /// The mean base quality (Phred) of all bases that directly support the variant call.
    pub base_quality_mean: f32,
    /// 1 when the supporting reads give the variant two or more distinct per-read qualities,
    /// otherwise 0; VarDict prints a flag here, not a standard deviation.
    pub stdev_base_stdev: f32,
    /// The two-sided Fisher exact p-value that the ALT reads' forward/reverse split differs from the
    /// REF reads' split.
    pub strand_bias_p_value: f32,
    #[serde(deserialize_with = "maybe_infinite_f32_odds_ratio")]
    /// The odds ratio for strand bias.
    pub strand_bias_odds_ratio: f32,
    /// The arithmetic mean mapping quality of the reads that support the variant call, uncapped.
    pub mean_mapping_quality: f32,
    /// The ratio of high- to low-quality ALT reads split at VarDict's -q, with 0.5 standing in for no
    /// low-quality reads; not a signal-to-noise ratio.
    pub signal_to_noise: f32,
    /// Allele frequency calculated using only high quality bases. Lossy due to rounding.
    pub af_high_quality_bases: f32,
    /// The fraction of VarDict's depth made of reads it reassigned to this allele by realignment or
    /// MNV merging, not an adjusted allele frequency.
    pub af_adjusted: f32,
    /// The number of bases an insertion or deletion can slide toward the 3' end with the same
    /// haplotype; VarDict also prints a meaningless value here for SNVs and MNVs.
    pub num_bases_3_prime_shift_for_deletions: i32,
    /// The number of copies, possibly fractional, of the 1-6 base repeat unit VarDict found at the
    /// variant; sequence context, not microsatellite instability.
    pub microsatellite: f32,
    /// The length in bases (0-6) of the repeat unit counted in `microsatellite`.
    pub microsatellite_length: i32,
    /// The mean substitution mismatches per read across the reads that support the variant call.
    pub mean_mismatches_in_reads: f32,
    /// The number of high quality reads supporting the variant call.
    pub high_quality_variant_reads: i32,
    /// The number of high quality reads at the locus of the variant call.
    pub high_quality_total_reads: i32,
    /// 5-prime reference flanking sequence.
    pub flank_seq_5_prime: &'a str,
    /// 3-prime reference flanking sequence.
    pub flank_seq_3_prime: &'a str,
    /// The position formatted interval of the variant calling target.
    pub segment: &'a str,
    /// The type of variant this call is.
    pub variant_type: &'a str,
    /// The fraction of reads VarDict's -t removed as duplicates; VarDictJava 1.8.4 always prints 0,
    /// which parses as `None`.
    #[serde(default, deserialize_with = "maybe_duplication_rate")]
    pub duplication_rate: Option<f32>,
    /// The details of the structural variant.
    #[serde(default, deserialize_with = "maybe_sv_info")]
    pub sv_info: Option<SvInfo>,
    #[serde(default)]
    /// The number of bases VarDict moved this indel toward the -J CRISPR cut site; only with -J.
    pub distance_to_crispr_site: Option<i32>,
}

impl<'a> TumorOnlyVariant<'a> {
    /// Return the "AD" formatted VCF field for this record: the REF depth alone without an ALT allele,
    /// otherwise REF then ALT, with REF missing for Complex calls.
    pub fn ad_value(&self) -> Vec<i32> {
        self.allele_depths(self.ref_forward + self.ref_reverse, self.alt_depth)
    }

    /// Return a per-allele depth for this record: the REF depth alone without an ALT allele,
    /// otherwise REF then ALT, with REF missing for Complex calls.
    fn allele_depths(&self, ref_depth: i32, alt_depth: i32) -> Vec<i32> {
        if self.ref_allele == self.alt_allele {
            vec![ref_depth]
        } else if self.variant_type == "Complex" {
            vec![i32::missing(), alt_depth]
        } else {
            vec![ref_depth, alt_depth]
        }
    }

    /// Return the "ADF" formatted VCF field for this record: the forward-strand half of AD.
    pub fn adf_value(&self) -> Vec<i32> {
        self.allele_depths(self.ref_forward, self.alt_forward)
    }

    /// Return the "ADR" formatted VCF field for this record: the reverse-strand half of AD.
    pub fn adr_value(&self) -> Vec<i32> {
        self.allele_depths(self.ref_reverse, self.alt_reverse)
    }

    /// Return the "HICNT" formatted VCF field for this record: missing without an ALT allele.
    pub fn hicnt_value(&self) -> i32 {
        if self.ref_allele == self.alt_allele {
            i32::missing()
        } else {
            self.high_quality_variant_reads
        }
    }

    /// Return the "REALIGNED_FRAC_OF_DP" formatted VCF field for this record: missing without an ALT
    /// allele.
    pub fn realigned_frac_of_dp_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele {
            f32::missing()
        } else {
            self.af_adjusted
        }
    }

    /// Return the "MEAN_DIST_TO_READ_END" formatted VCF field for this record: missing without an
    /// ALT allele.
    pub fn mean_dist_to_read_end_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele {
            f32::missing()
        } else {
            self.mean_position_in_read
        }
    }

    /// Return the "ALT_READ_POS_VARIES" formatted VCF field for this record: missing without an ALT
    /// allele.
    pub fn alt_read_pos_varies_value(&self) -> i32 {
        if self.ref_allele == self.alt_allele {
            i32::missing()
        } else {
            self.stdev_position_in_read as i32
        }
    }

    /// Return the "MEAN_MAPQ" formatted VCF field for this record: missing without an ALT allele.
    pub fn mean_mapq_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele {
            f32::missing()
        } else {
            self.mean_mapping_quality
        }
    }

    /// Return the "STRAND_BIAS_FISHER_P" formatted VCF field for this record: missing without an ALT
    /// allele.
    pub fn strand_bias_fisher_p_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele {
            f32::missing()
        } else {
            self.strand_bias_p_value
        }
    }

    /// Return the "QMEAN" formatted VCF field for this record: missing without an ALT allele.
    pub fn qmean_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele {
            f32::missing()
        } else {
            self.base_quality_mean
        }
    }

    /// Return the "MEAN_MISMATCHES" formatted VCF field for this record: missing without an ALT allele.
    pub fn mean_mismatches_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele {
            f32::missing()
        } else {
            self.mean_mismatches_in_reads
        }
    }

    /// Return the "AF" formatted VCF field for this record: the ALT depth over DP at full precision,
    /// clamped to [0, 1], and missing without an ALT allele or without depth.
    pub fn af_value(&self) -> f32 {
        if self.ref_allele == self.alt_allele || self.depth <= 0 {
            return f32::missing();
        }
        (self.alt_depth as f64 / self.depth as f64).clamp(0.0, 1.0) as f32
    }

    /// Return the VCF-valid alternate allele for this record.
    pub fn alt_allele_for_vcf(&self) -> String {
        if self.ref_allele == self.alt_allele {
            String::from(".")
        } else {
            String::from(self.alt_allele)
        }
    }

    /// Return the "GT" formatted VCF field for this record: 0/0 without an ALT allele, otherwise 1/1
    /// at or above the minimum homozygous alternate allele frequency and 0/1 below it.
    pub fn gt_value(&self, min_hom_alt_af: f32) -> &[GenotypeAllele] {
        if self.ref_allele == self.alt_allele {
            &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
        } else if self.af_value() >= min_hom_alt_af {
            &[GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)]
        } else {
            &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]
        }
    }

    /// Return VarDict's variant class for the alleles this record writes, using VarDict's own rule
    /// (`Variant.varType`) on the final alleles: VarDict labels a Complex call before trimming it.
    pub fn variant_class(&self) -> Option<&'a str> {
        let (ref_allele, alt_allele) = (self.ref_allele, self.alt_allele);
        if ref_allele == alt_allele {
            return None;
        }
        if ref_allele.len() == 1 && alt_allele.len() == 1 {
            return Some("SNV");
        }
        if let Some(sv) = alt_allele
            .strip_prefix('<')
            .and_then(|rest| rest.strip_suffix('>'))
            .filter(|sv| sv.len() == 3)
        {
            return Some(sv);
        }
        let class = match (ref_allele.as_bytes().first(), alt_allele.as_bytes().first()) {
            (Some(r), Some(a)) if r == a => {
                if ref_allele.len() == 1 && alt_allele.starts_with(ref_allele) {
                    "Insertion"
                } else if alt_allele.len() == 1 && ref_allele.starts_with(alt_allele) {
                    "Deletion"
                } else {
                    "Complex"
                }
            }
            _ => "Complex",
        };
        Some(class)
    }

    /// Return the "INDEL_3P_SHIFT" formatted VCF field for this record: only for insertions and
    /// deletions, where VarDict's 3' shift describes the indel.
    pub fn indel_3p_shift_value(&self) -> Option<i32> {
        match self.variant_class() {
            Some("Insertion") | Some("Deletion") => {
                Some(self.num_bases_3_prime_shift_for_deletions)
            }
            _ => None,
        }
    }

    /// Return the "REPEAT_UNIT_COPIES" formatted VCF field for this record: absent without an ALT
    /// allele or when VarDict did not compute it, which it marks with 0.
    pub fn repeat_unit_copies_value(&self) -> Option<f32> {
        self.has_repeat_context().then_some(self.microsatellite)
    }

    /// Return the "REPEAT_UNIT_LEN" formatted VCF field for this record, under the same rule.
    pub fn repeat_unit_len_value(&self) -> Option<i32> {
        self.has_repeat_context()
            .then_some(self.microsatellite_length)
    }

    /// Whether VarDict computed the repeat context for this record's ALT allele.
    fn has_repeat_context(&self) -> bool {
        self.ref_allele != self.alt_allele && self.microsatellite != 0.0
    }

    /// Return the structural variant's soft-clipped and discordant read counts, only for a symbolic
    /// ALT allele; VarDict attaches its SV column to every row at the position.
    pub fn sv_read_counts(&self) -> Option<(i32, i32)> {
        if !self.alt_allele.starts_with('<') {
            return None;
        }
        let info = self.sv_info.as_ref();
        Some((
            info.map_or(0, |sv| sv.supporting_split_reads),
            info.map_or(0, |sv| sv.supporting_pairs),
        ))
    }

    /// Return the signed structural variant length from VarDict's own event length, which it writes
    /// after the slash in the genotype column: `-N` for DEL, `+N` for DUP and `<INVN>` for INV.
    pub fn sv_length(&self) -> Option<i32> {
        let event = self.gt.split_once('/')?.1;
        let (sign, rest) = match self.variant_type {
            "DEL" => (-1, event.strip_prefix('-')?),
            "DUP" => (1, event.strip_prefix('+')?),
            "INV" => (1, event.strip_prefix("<INV")?),
            _ => return None,
        };
        let digits = rest
            .find(|c: char| !c.is_ascii_digit())
            .map_or(rest, |end| &rest[..end]);
        digits.parse::<i32>().ok().map(|length| sign * length)
    }
}

impl<'a> AbstractInterval for TumorOnlyVariant<'a> {
    fn contig(&self) -> &str {
        self.contig
    }

    fn range(&self) -> Range<Position> {
        Range {
            start: self.start - 1,
            end: self.end,
        }
    }
}

/// Create a VCF header for VarDict/VarDictJava in tumor-only mode.
#[rustfmt::skip]
pub fn tumor_only_header(sample: &str) -> Header {
    let source = [CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-");
    let mut header = Header::default();
    header.push_sample(sample.as_bytes());
    header.remove_filter(b"PASS");
    header.push_record(format!("##source={source}").as_bytes());
    header.push_record(r#"##INFO=<ID=TYPE,Number=A,Type=String,Description="VarDict's class of the change from REF to ALT, from VarDict's own rule applied to the alleles in this record: SNV (one base to one base), Insertion (ALT is the single REF base followed by inserted bases), Deletion (REF is longer and ALT is its first base), Complex (every other change, including MNVs, which VarDict does not separate), or DEL, DUP or INV for a symbolic structural variant. Absent when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=INDEL_3P_SHIFT,Number=A,Type=Integer,Description="Number of bases this insertion or deletion can slide toward the 3' end (rightward on the forward strand) and still describe the same haplotype, counted within VarDict's 70-base window, so values top out near 70; POS plus this is the right-most equivalent position. Written only for insertions and deletions.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=REPEAT_UNIT_COPIES,Number=A,Type=Float,Description="Copies of the 1-6 bp tandem repeat unit next to the variant, counted by VarDict in the reference or ALT haplotype, whichever is larger (1 means no repeat). VarDict takes the unit from one side of the variant (for SNVs, the repeat starting at the following base), so a repeat on the other side can be missed; it also raises the count to an indel's 3' shift divided by its length when that is larger, which can make it fractional and leave REPEAT_UNIT_LEN describing a different repeat. Sequence context, not microsatellite instability. Rounded by VarDict to 3 decimals; absent when VarDict did not compute it.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=REPEAT_UNIT_LEN,Number=A,Type=Integer,Description="Length in bp (1-6) of the repeat unit counted in REPEAT_UNIT_COPIES; absent when VarDict did not compute it.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position, written only on records with a symbolic ALT allele: the last deleted base for DEL, the last inverted base for INV, and VarDict's end of the duplication for DUP, which can be off by one depending on how VarDict found it.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Signed length of the structural variant, written only on records with a symbolic ALT allele: VarDict's own event length, negative for DEL (the deleted bases), positive for DUP (the duplicated bases) and INV (the inverted bases).">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type, written only on records with a symbolic ALT allele: DEL, DUP or INV.">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=PASS,Description="The variant call has passed all filters and may be considered for downstream analysis.">"#.as_bytes());
    header.push_record(format!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="The genotype inferred from the allele frequency alone, since VarDict does not genotype: 1/1 when AF >= {MIN_HOM_ALT_AF}, 0/1 for any other call with an ALT allele, and 0/0 when there is no ALT allele.">"#).as_bytes());
    header.push_record(r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for REF then ALT as VarDict counts them: ALT is the reads carrying this allele and REF is the reads carrying the reference base at the variant's first internal base (the first deleted base for deletions, without insertion-carrying reads for insertions), or at the following base when an insertion at the same position raises DP. REF is missing for Complex calls, where VarDict counts only their first base and can count a read as both, and REF is the only value when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Reads on the forward strand (SAM flag 0x10 unset) supporting REF then ALT, the forward half of AD as VarDict counts it: REF is counted where AD's is (sometimes the following base), is missing for Complex calls and is the only value when there is no ALT allele. When VarDict runs with -u, overlapping mates count only through the reverse read.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Reads on the reverse strand (SAM flag 0x10 set) supporting REF then ALT, the reverse half of AD as VarDict counts it: REF is counted where AD's is (sometimes the following base), is missing for Complex calls and is the only value when there is no ALT allele. When VarDict runs with -u, overlapping mates count only through the reverse read.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=STRAND_BIAS_FISHER_P,Number=A,Type=Float,Description="Two-sided Fisher exact p-value that the ALT allele's forward/reverse read split differs from REF's, from VarDict's table of REF and ALT reads by SAM strand, not a test against 50:50; its REF counts are VarDict's, including the overcount on Complex calls. Rounded by VarDict to 5 decimals, so values below 0.000005 read 0. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth as VarDict counts it at the variant's first internal base (POS+1 for deletions), the denominator of FORMAT AF and usually of VarDict's own AF: reads with any base quality, reads whose deletion spans the base, reference-matching soft-clipped bases, reads VarDict reassigned by realignment, and N calls only under -K; overlapping mates count twice unless VarDict ran with -u or -UN, and an insertion at the same position can make it the following base's depth. REF and ALT depths need not sum to DP.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Fraction of reads carrying the ALT allele, AD[1] / DP, computed by vartovcf at full precision and clamped to [0, 1] because VarDict's ALT count can exceed DP; AD and DP keep the raw counts. VarDict's own AF column is rounded to 4 decimals and its denominator can differ from DP, for example at a position shared with an insertion. Missing when there is no ALT allele or DP is 0.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=HICNT,Number=A,Type=Integer,Description="ALT reads whose per-read variant quality, the one QMEAN averages (the base's quality for SNVs, the better flanking base's for deletions and so on), is at least VarDict's -q; mapping quality is not considered. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=REALIGNED_FRAC_OF_DP,Number=A,Type=Float,Description="Fraction of VarDict's depth made of reads it reassigned to this ALT allele by local realignment or MNV merging (VarDict's ExtraAF, which var2vcf_valid.pl calls ADJAF). Those reads are already counted in AD[1] and AF, so this is not an adjusted AF; divided by AF it gives roughly the share of the ALT reads that realignment contributed. Rounded by VarDict to 4 decimals. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=MEAN_DIST_TO_READ_END,Number=A,Type=Float,Description="Mean, over the reads carrying the ALT allele, of the 1-based distance from the variant to the nearer end of the read's aligned part, soft clips excluded. It is measured from the variant base for SNVs, the first inserted base for insertions and the first base after the gap for deletions; a Complex call is measured like the insertion or deletion it starts with, or from the last base of the block when it starts with a mismatch. Reads VarDict realigned out of soft clips contribute their clip length. Rounded by VarDict to 1 decimal. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ALT_READ_POS_VARIES,Number=A,Type=Integer,Description="1 when the reads carrying the ALT allele place it at two or more distinct distances from the read end, and 0 when every one has it at the same distance, which with several reads is a common sign of an artifact. VarDict also sets it to 1 whenever it reassigns reads to the allele, as in local realignment or MNV merging, so it is 1 whenever REALIGNED_FRAC_OF_DP is above 0; otherwise a call with one ALT read is always 0. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=QMEAN,Number=A,Type=Float,Description="Mean, over the reads carrying the ALT allele, of VarDict's per-read variant quality: the base's Phred quality for SNVs, the mean of the inserted bases for insertions, the higher of the two flanking bases for deletions, and the mean of the block's bases for Complex calls, with that flanking base averaged in when the block starts with a deletion. VarDict grows an MNV only through mismatches at least 5 above -q, so MNV values run high. Bases below -q are included and nothing is capped. Rounded by VarDict to 1 decimal. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=MEAN_MAPQ,Number=A,Type=Float,Description="Arithmetic mean mapping quality of the reads carrying the ALT allele, not the RMS over all reads that the VCF specification's MQ means; uncapped, so a MAPQ of 255 (unavailable) counts as 255. Rounded by VarDict to 1 decimal. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=MEAN_MISMATCHES,Number=A,Type=Float,Description="Mean, over the reads carrying the ALT allele, of each read's substitution mismatches: its NM tag minus every inserted or deleted base, and for Complex calls minus some of the mismatches VarDict folded into the allele after its first change. Includes the variant's own mismatch for SNVs and excludes soft clips; a read without an NM tag counts as 0, and reads with more than VarDict's -m mismatches (default 8) are not counted. Rounded by VarDict to 1 decimal. Missing when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=SV_SOFTCLIP_READS,Number=A,Type=Integer,Description="Reads soft-clipped at this structural variant's breakpoint whose clipped sequence VarDict matched to the other side; primary alignments only, as VarDict ignores supplementary (SA) records. Written only on records with a symbolic ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=SV_DISCORDANT_READS,Number=A,Type=Integer,Description="Discordant read records (unexpected insert size or orientation) in the clusters VarDict linked to this structural variant; each mate counts, so one pair can count twice. Written only on records with a symbolic ALT allele.">"#.as_bytes());
    header.push_record(r#"##ALT=<ID=DEL,Description="Deletion relative to the reference.">"#.as_bytes());
    header.push_record(r#"##ALT=<ID=DUP,Description="Region of elevated copy number relative to the reference.">"#.as_bytes());
    header.push_record(r#"##ALT=<ID=INV,Description="Inversion of reference sequence.">"#.as_bytes());
    header
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use rust_htslib::bcf::{Format, Read};
    use rust_htslib::bcf::{Reader as VcfReader, Writer as VcfWriter};
    use tempfile::NamedTempFile;

    use super::*;

    #[fixture]
    #[rustfmt::skip]
    fn variants() -> Vec<TumorOnlyVariant<'static>> {
        let inv_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        let dup_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        let del_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        vec![
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<INV>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<INV11063372>", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 58.0, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 33.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "INV", duplication_rate: None, sv_info: Some(inv_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<DUP>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<DUP11063372>", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 58.0, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 10.0, mean_mapping_quality: 33.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "DUP", duplication_rate: Some(0.001), sv_info: Some(dup_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "A", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 1, alt_reverse: 0, gt: "G/A", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 13.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, strand_bias_p_value: 0.34385, strand_bias_odds_ratio: 4.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1.0, microsatellite_length: 1, mean_mismatches_in_reads: 2.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "FLT3", contig: "chr13", start: 24684729, end: 28034141, ref_allele: "G", alt_allele: "<DEL>", depth: 7463, alt_depth: 2, ref_forward: 7463, ref_reverse: 0, alt_forward: 0, alt_reverse: 2, gt: "-3349412/-3349412", af: 0.0003, strand_bias: PairBias::from_str("0;0").unwrap(), mean_position_in_read: 60.5, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 2.0, mean_mapping_quality: 60.0, signal_to_noise: 4.0, af_high_quality_bases: 1.0000, af_adjusted: 0.0003, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 2, high_quality_total_reads: 2, flank_seq_5_prime: "TGCTGTAGTCTAATGATTCT", flank_seq_3_prime: "CAACGTAGAAGTACTCATTA", segment: "chr13:28033879-28034298", variant_type: "DEL", duplication_rate: None, sv_info: Some(del_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "T", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 0, alt_reverse: 1, gt: "G/T", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 28.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1.0, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713880, end: 114713880, ref_allele: "T", alt_allele: "A", depth: 8211, alt_depth: 1, ref_forward: 3001, ref_reverse: 5130, alt_forward: 1, alt_reverse: 0, gt: "T/A", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 18.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, strand_bias_p_value: 0.36916, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 2.0, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8132, flank_seq_5_prime: "CCTTCGCCTGTCCTCATGTA", flank_seq_3_prime: "TGGTCTCTCATGGCACTGTA", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
        ]
    }

    #[test]
    fn test_pair_bias_default() {
        let expected = PairBias {
            reference: StrandBias::UnDetected,
            alternate: StrandBias::UnDetected,
        };
        assert_eq!(PairBias::default(), expected);
    }

    #[rstest(
        int,
        expected,
        case("0", StrandBias::TooFewReads),
        case("1", StrandBias::Detected),
        case("2", StrandBias::UnDetected)
    )]
    fn test_strand_bias_from_str(int: &str, expected: StrandBias) {
        assert_eq!(StrandBias::from_str(int).expect("Parse failed!"), expected);
    }

    #[rstest(
        left => ["0", "1", "2"],
        right => ["0", "1", "2"]
    )]
    fn test_pair_bias_from_str(left: &str, right: &str) {
        assert!(PairBias::from_str(&format!("{};{}", left, right)).is_ok());
    }

    #[test]
    fn test_parse_pair_bias_err_display() {
        assert_eq!(&ParsePairBiasError.to_string(), "ParsePairBiasError");
    }

    #[test]
    fn test_pair_bias_from_str_err() {
        assert_eq!(PairBias::from_str("1"), Err(ParsePairBiasError));
        assert_eq!(PairBias::from_str("3;0"), Err(ParsePairBiasError));
        assert_eq!(PairBias::from_str("1;1;1"), Err(ParsePairBiasError));
    }

    #[test]
    fn test_sv_info_from_str_err() {
        assert_eq!(SvInfo::from_str("0-0"), Err(ParseSvInfoError));
    }

    #[test]
    fn test_parse_sv_info_err_display() {
        assert_eq!(&ParseSvInfoError.to_string(), "ParseSvInfoError");
    }

    #[test]
    fn test_sv_info_default() {
        let expected = SvInfo {
            supporting_split_reads: 0,
            supporting_pairs: 0,
            supporting_clusters: 0,
        };
        assert_eq!(SvInfo::default(), expected);
    }

    #[test]
    fn test_sv_info_to_string() {
        let sv_info = SvInfo {
            supporting_split_reads: 4,
            supporting_pairs: 5,
            supporting_clusters: 9,
        };
        assert_eq!(sv_info.to_string(), "4-5-9");
    }

    #[rstest]
    fn test_tumor_only_variant_ad_value(variants: Vec<TumorOnlyVariant>) {
        let expected = [
            vec![6579, 1],
            vec![6579, 1],
            vec![8046, 1],
            vec![7463, 2],
            vec![8046, 1],
            vec![8131, 1],
        ];
        for (variant, ad) in variants.iter().zip(expected.iter()) {
            assert_eq!(&variant.ad_value(), ad);
        }
    }

    #[rstest]
    #[case(1, 1_000_000, 1e-6)]
    #[case(9, 22, 9.0 / 22.0)]
    #[case(12, 10, 1.0)]
    #[case(0, 10, 0.0)]
    fn test_tumor_only_variant_af_value(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] alt_depth: i32,
        #[case] depth: i32,
        #[case] expected: f64,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_depth = alt_depth;
        variant.depth = depth;
        variant.af = 0.0;
        assert_eq!(variant.af_value(), expected as f32);
    }

    #[rstest]
    #[case("G", "A", 1, 0)]
    #[case("G", "G", 9, 10)]
    fn test_tumor_only_variant_af_value_is_missing(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] alt_depth: i32,
        #[case] depth: i32,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.alt_depth = alt_depth;
        variant.depth = depth;
        assert!(variant.af_value().is_missing());
    }

    #[rstest]
    #[case("G", "A", 12.0, 1, Some(12.0), Some(1))]
    #[case("GTT", "G", 1.571, 1, Some(1.571), Some(1))]
    #[case("G", "A", 0.0, 0, None, None)]
    #[case("G", "G", 12.0, 1, None, None)]
    fn test_tumor_only_variant_repeat_unit_values(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] copies: f32,
        #[case] unit_length: i32,
        #[case] expected_copies: Option<f32>,
        #[case] expected_length: Option<i32>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.microsatellite = copies;
        variant.microsatellite_length = unit_length;
        assert_eq!(variant.repeat_unit_copies_value(), expected_copies);
        assert_eq!(variant.repeat_unit_len_value(), expected_length);
    }

    #[rstest]
    fn test_tumor_only_variant_sv_read_counts(variants: Vec<TumorOnlyVariant<'static>>) {
        let variants: Vec<TumorOnlyVariant> = variants.into_iter().collect();
        assert_eq!(variants[0].sv_read_counts(), Some((1, 1)));
        assert_eq!(variants[2].sv_read_counts(), None);
    }

    #[rstest]
    fn test_tumor_only_variant_sv_read_counts_without_sv_info(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().next().unwrap();
        variant.sv_info = None;
        assert_eq!(variant.sv_read_counts(), Some((0, 0)));
    }

    #[rstest]
    #[case("GTT", "G", Some(3))]
    #[case("G", "GTT", Some(3))]
    #[case("G", "A", None)]
    #[case("GA", "AC", None)]
    #[case("A", "<DEL>", None)]
    #[case("G", "G", None)]
    fn test_tumor_only_variant_indel_3p_shift_value(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] expected: Option<i32>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.num_bases_3_prime_shift_for_deletions = 3;
        assert_eq!(variant.indel_3p_shift_value(), expected);
    }

    #[rstest]
    #[case("G", "A", "SNV", Some("SNV"))]
    #[case("G", "GTT", "Insertion", Some("Insertion"))]
    #[case("GTT", "G", "Deletion", Some("Deletion"))]
    #[case("GA", "AC", "Complex", Some("Complex"))]
    #[case("GAT", "GCC", "Complex", Some("Complex"))]
    #[case("A", "TA", "Complex", Some("Complex"))]
    #[case("A", "C", "Complex", Some("SNV"))]
    #[case("A", "<DEL>", "DEL", Some("DEL"))]
    #[case("G", "G", "", None)]
    fn test_tumor_only_variant_variant_class(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] variant_type: &'static str,
        #[case] expected: Option<&str>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.variant_type = variant_type;
        assert_eq!(variant.variant_class(), expected);
    }

    #[rstest]
    #[case("DEL", "-3349412/-3349412", Some(-3349412))]
    #[case("DUP", "A/+11063372", Some(11063372))]
    #[case("INV", "A/<INV11063372>", Some(11063372))]
    #[case("DEL", "G/G", None)]
    #[case("SNV", "G/A", None)]
    fn test_tumor_only_variant_sv_length(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] variant_type: &'static str,
        #[case] gt: &'static str,
        #[case] expected: Option<i32>,
    ) {
        let mut variant = variants.into_iter().nth(3).unwrap();
        variant.variant_type = variant_type;
        variant.gt = gt;
        assert_eq!(variant.sv_length(), expected);
    }

    #[rstest]
    fn test_tumor_only_variant_read_position_values_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.mean_position_in_read = 10.5;
        variant.stdev_position_in_read = 1.0;
        assert_eq!(variant.mean_dist_to_read_end_value(), 10.5);
        assert_eq!(variant.alt_read_pos_varies_value(), 1);
    }

    #[rstest]
    fn test_tumor_only_variant_read_position_values_are_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.mean_dist_to_read_end_value().is_missing());
        assert!(variant.alt_read_pos_varies_value().is_missing());
    }

    #[rstest]
    fn test_tumor_only_variant_realigned_frac_of_dp_value_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.af_adjusted = 0.1364;
        assert_eq!(variant.realigned_frac_of_dp_value(), 0.1364);
    }

    #[rstest]
    fn test_tumor_only_variant_realigned_frac_of_dp_value_is_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.realigned_frac_of_dp_value().is_missing());
    }

    #[rstest]
    fn test_tumor_only_variant_hicnt_value_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.high_quality_variant_reads = 7;
        assert_eq!(variant.hicnt_value(), 7);
    }

    #[rstest]
    fn test_tumor_only_variant_hicnt_value_is_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.hicnt_value().is_missing());
    }

    #[rstest]
    fn test_tumor_only_variant_mean_mismatches_value_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.mean_mismatches_in_reads = 1.5;
        assert_eq!(variant.mean_mismatches_value(), 1.5);
    }

    #[rstest]
    fn test_tumor_only_variant_mean_mismatches_value_is_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.mean_mismatches_value().is_missing());
    }

    #[rstest]
    fn test_tumor_only_variant_strand_bias_fisher_p_value_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.strand_bias_p_value = 0.34385;
        assert_eq!(variant.strand_bias_fisher_p_value(), 0.34385);
    }

    #[rstest]
    fn test_tumor_only_variant_strand_bias_fisher_p_value_is_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.strand_bias_fisher_p_value().is_missing());
    }

    #[rstest]
    fn test_tumor_only_variant_mean_mapq_value_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.mean_mapping_quality = 52.5;
        assert_eq!(variant.mean_mapq_value(), 52.5);
    }

    #[rstest]
    fn test_tumor_only_variant_mean_mapq_value_is_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.mean_mapq_value().is_missing());
    }

    #[rstest]
    fn test_tumor_only_variant_qmean_value_on_an_alt_call(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.base_quality_mean = 37.5;
        assert_eq!(variant.qmean_value(), 37.5);
    }

    #[rstest]
    fn test_tumor_only_variant_qmean_value_is_missing_without_an_alt(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_allele = variant.ref_allele;
        assert!(variant.qmean_value().is_missing());
    }

    #[rstest]
    #[case("G", "A", "SNV", vec![2766, 1], vec![5280, 0])]
    #[case("GA", "AC", "Complex", vec![i32::missing(), 1], vec![i32::missing(), 0])]
    #[case("G", "G", "", vec![2766], vec![5280])]
    fn test_tumor_only_variant_strand_depths(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] variant_type: &'static str,
        #[case] adf: Vec<i32>,
        #[case] adr: Vec<i32>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.variant_type = variant_type;
        assert_eq!(variant.adf_value(), adf);
        assert_eq!(variant.adr_value(), adr);
    }

    #[rstest]
    #[case("G", "A", "SNV", 1, vec![8046, 1])]
    #[case("G", "A", "SNV", 0, vec![8046, 0])]
    #[case("GA", "AC", "Complex", 9, vec![i32::missing(), 9])]
    #[case("G", "G", "", 0, vec![8046])]
    fn test_tumor_only_variant_ad_value_by_call(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] variant_type: &'static str,
        #[case] alt_depth: i32,
        #[case] expected: Vec<i32>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.variant_type = variant_type;
        variant.alt_depth = alt_depth;
        assert_eq!(variant.ad_value(), expected);
    }

    #[rstest]
    #[case("G", "A", 1, 10000, [0, 1])]
    #[case("G", "A", 25, 100, [0, 1])]
    #[case("G", "A", 7999, 10000, [0, 1])]
    #[case("G", "A", 8, 10, [1, 1])]
    #[case("G", "A", 12, 10, [1, 1])]
    #[case("G", "A", 1, 0, [0, 1])]
    #[case("G", "G", 9, 10, [0, 0])]
    fn test_tumor_only_variant_gt_value(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] alt_depth: i32,
        #[case] depth: i32,
        #[case] expected: [i32; 2],
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.alt_depth = alt_depth;
        variant.depth = depth;
        variant.af = 0.0;
        assert_eq!(
            variant.gt_value(MIN_HOM_ALT_AF),
            &expected.map(GenotypeAllele::Unphased)
        );
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_tumor_only_variant_interval(variants: Vec<TumorOnlyVariant>) {
        let expected = [
            ("chr12", Range { start: 112450447 - 1, end: 123513818 } ),
            ("chr12", Range { start: 112450447 - 1, end: 123513818 } ),
            ("chr1", Range { start: 114713883 - 1, end: 114713883 } ),
            ("chr13", Range { start: 24684729 - 1, end: 28034141 }),
            ("chr1", Range { start: 114713883 - 1, end: 114713883 } ),
            ("chr1", Range { start: 114713880 - 1, end: 114713880 } ),
        ];
        for (variant, (contig, range)) in variants.iter().zip(expected.iter()) {
            assert_eq!(&variant.contig(), contig);
            assert_eq!(&variant.range(), range);
        }
    }

    #[rstest]
    #[case("Inf", 0.0)]
    #[case("inf", 0.0)]
    #[case("0.5", 2.0)]
    #[case("0.25", 4.0)]
    #[case("0.1", 10.0)]
    #[case("10.0", 10.0)]
    #[case("2.0", 2.0)]
    #[case("1.0", 1.0)]
    #[case("0.0", 0.0)]
    fn test_maybe_infinite_f32_odds_ratio(#[case] input: &'static str, #[case] expected: f32) {
        use serde::Deserialize;
        use serde_test::{Token, assert_de_tokens};

        #[derive(Debug, Deserialize, PartialEq)]
        struct TestStruct {
            #[serde(deserialize_with = "maybe_infinite_f32_odds_ratio")]
            value: f32,
        }

        assert_de_tokens(
            &TestStruct { value: expected },
            &[
                Token::Struct {
                    name: "TestStruct",
                    len: 1,
                },
                Token::BorrowedStr("value"),
                Token::BorrowedStr(input),
                Token::StructEnd,
            ],
        );
    }

    #[test]
    fn test_tumor_only_header() {
        let header = tumor_only_header("dna00001");
        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let _ = VcfWriter::from_path(file.path(), &header, true, Format::Vcf).unwrap();
        let reader = VcfReader::from_path(file.path()).expect("Error opening tempfile!");
        let records = reader.header().header_records();
        let samples = reader.header().samples();
        assert_eq!(records.len(), 29);
        assert_eq!(samples.len(), 1);
        assert!(samples.iter().all(|&s| s == "dna00001".as_bytes()));
    }
}
