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

use crate::filter::FilterThresholds;
use crate::fisher::fisher_exact_greater;
use crate::hicnt::{hicnt_from_sn, parse_thousandths};
use rust_htslib::bcf::record::{GenotypeAllele, Numeric};
use serde::{Deserialize, Serialize, de::Error};

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

/// Deserialize VarDict's SN, which it prints with at most 3 decimals, exactly as thousandths.
fn thousandths<'de, D>(deserializer: D) -> Result<u64, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let string: &str = Deserialize::deserialize(deserializer)?;
    parse_thousandths(string).ok_or_else(|| {
        D::Error::custom(format!(
            "expected a number with at most 3 decimals, found '{string}'"
        ))
    })
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

/// A container for structural variant (SV) information.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
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

/// A record of output from VarDict/VarDictJava run in tumor-only mode, or one sample's half of a
/// tumor-normal record.
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
    /// VarDict's strand-bias flags, which no output field uses, so they are not parsed.
    pub strand_bias: &'a str,
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
    #[serde(deserialize_with = "maybe_duplication_rate")]
    pub duplication_rate: Option<f32>,
    /// The details of the structural variant.
    #[serde(deserialize_with = "maybe_sv_info")]
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

    /// Return the "REALIGNED_FRAC_OF_DP" formatted VCF field for this record: missing without ALT
    /// reads.
    pub fn realigned_frac_of_dp_value(&self) -> f32 {
        if !self.has_alt_reads() {
            f32::missing()
        } else {
            self.af_adjusted
        }
    }

    /// Return the "MEAN_DIST_TO_READ_END" formatted VCF field for this record: missing without
    /// ALT reads.
    pub fn mean_dist_to_read_end_value(&self) -> f32 {
        if !self.has_alt_reads() {
            f32::missing()
        } else {
            self.mean_position_in_read
        }
    }

    /// Return the "ALT_READ_POS_VARIES" formatted VCF field for this record: missing without ALT
    /// reads.
    pub fn alt_read_pos_varies_value(&self) -> i32 {
        if !self.has_alt_reads() {
            i32::missing()
        } else {
            self.stdev_position_in_read as i32
        }
    }

    /// Return the "MEAN_MAPQ" formatted VCF field for this record: missing without ALT reads.
    pub fn mean_mapq_value(&self) -> f32 {
        if !self.has_alt_reads() {
            f32::missing()
        } else {
            self.mean_mapping_quality
        }
    }

    /// Return the "STRAND_BIAS_FISHER_P" formatted VCF field for this record: missing without ALT
    /// reads.
    pub fn strand_bias_fisher_p_value(&self) -> f32 {
        if !self.has_alt_reads() {
            f32::missing()
        } else {
            self.strand_bias_p_value
        }
    }

    /// Return the "QMEAN" formatted VCF field for this record: missing without ALT reads.
    pub fn qmean_value(&self) -> f32 {
        if !self.has_alt_reads() {
            f32::missing()
        } else {
            self.base_quality_mean
        }
    }

    /// Return the "MEAN_MISMATCHES" formatted VCF field for this record: missing without ALT reads.
    pub fn mean_mismatches_value(&self) -> f32 {
        if !self.has_alt_reads() {
            f32::missing()
        } else {
            self.mean_mismatches_in_reads
        }
    }

    /// Whether this record has an ALT allele and at least one read carrying it.
    pub fn has_alt_reads(&self) -> bool {
        self.ref_allele != self.alt_allele && self.alt_depth > 0
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

    /// Return the "GT" formatted VCF field for this record: ./. without depth or ALT reads, 0/0 with
    /// depth but no ALT reads, otherwise 1/1 at or above the minimum homozygous alternate allele
    /// frequency and 0/1 below it.
    pub fn gt_value(&self, min_hom_alt_af: f32) -> &[GenotypeAllele] {
        if !self.has_alt_reads() && self.depth <= 0 {
            &[
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing,
            ]
        } else if !self.has_alt_reads() {
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

/// One sample's columns in a record of VarDict/VarDictJava output run on a tumor and normal pair.
#[derive(Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct SampleColumns<'a> {
    /// The total depth at this call locus.
    pub depth: i32,
    /// The number of reads supporting the alternate allele.
    pub alt_depth: i32,
    /// The number of forward reads supporting the reference allele.
    pub ref_forward: i32,
    /// The number of reverse reads supporting the reference allele.
    pub ref_reverse: i32,
    /// The number of forward reads supporting the alternate allele.
    pub alt_forward: i32,
    /// The number of reverse reads supporting the alternate allele.
    pub alt_reverse: i32,
    /// VarDict's genotype column, `0` when VarDict has no reads of this allele to describe in this
    /// sample, as in an empty block or one it filled from another allele's depth.
    pub gt: &'a str,
    /// VarDict's allele frequency of the alternate allele.
    pub af: f32,
    /// VarDict's strand-bias flags, which no output field uses, so they are not parsed.
    pub strand_bias: &'a str,
    /// The mean distance from the variant to the nearer end of each supporting read's aligned part.
    pub mean_position_in_read: f32,
    /// 1 when the supporting reads place the variant at two or more distances from the read end.
    pub stdev_position_in_read: f32,
    /// The mean per-read variant quality of the supporting reads.
    pub base_quality_mean: f32,
    /// 1 when the supporting reads give the variant two or more distinct per-read qualities.
    pub stdev_base_stdev: f32,
    /// The mean mapping quality of the supporting reads.
    pub mean_mapping_quality: f32,
    /// VarDict's SN, the ratio of high- to low-quality ALT reads split at its -q, with 0.5 standing in
    /// for no low-quality reads, held exactly in thousandths as VarDict prints it to 3 decimals.
    #[serde(deserialize_with = "thousandths")]
    pub signal_to_noise_thousandths: u64,
    /// Allele frequency calculated using only high quality bases.
    pub af_high_quality_bases: f32,
    /// The fraction of VarDict's depth made of reads it reassigned to this allele.
    pub af_adjusted: f32,
    /// The mean substitution mismatches per supporting read.
    pub mean_mismatches_in_reads: f32,
    /// The two-sided Fisher exact p-value of REF and ALT reads by strand.
    pub strand_bias_p_value: f32,
    /// The odds ratio for strand bias.
    #[serde(deserialize_with = "maybe_infinite_f32_odds_ratio")]
    pub strand_bias_odds_ratio: f32,
}

/// A record of output from VarDict/VarDictJava run on a tumor and normal pair with `-b "T|N"`.
#[derive(Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct TumorNormalVariant<'a> {
    /// VarDict's -N value, `tumor|normal` when given both names with a BED file of regions.
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
    /// The tumor's columns.
    #[serde(borrow)]
    pub tumor_columns: SampleColumns<'a>,
    /// The normal's columns.
    #[serde(borrow)]
    pub normal_columns: SampleColumns<'a>,
    /// The number of bases an insertion or deletion can slide toward the 3' end.
    pub num_bases_3_prime_shift_for_deletions: i32,
    /// The number of copies of the 1-6 base repeat unit VarDict found at the variant.
    pub microsatellite: f32,
    /// The length in bases (0-6) of the repeat unit counted in `microsatellite`.
    pub microsatellite_length: i32,
    /// 5-prime reference flanking sequence.
    pub flank_seq_5_prime: &'a str,
    /// 3-prime reference flanking sequence.
    pub flank_seq_3_prime: &'a str,
    /// The position formatted interval of the variant calling target.
    pub segment: &'a str,
    /// VarDict's tumor-normal label for this allele, such as StrongSomatic or Germline.
    pub status: &'a str,
    /// The type of variant this call is.
    pub variant_type: &'a str,
    /// The tumor's duplication rate, which VarDictJava 1.8.4 always prints as 0.
    #[serde(deserialize_with = "maybe_duplication_rate")]
    pub tumor_duplication_rate: Option<f32>,
    /// The tumor's structural variant read counts.
    #[serde(deserialize_with = "maybe_sv_info")]
    pub tumor_sv_info: Option<SvInfo>,
    /// The normal's duplication rate, which VarDictJava 1.8.4 always prints as 0.
    #[serde(deserialize_with = "maybe_duplication_rate")]
    pub normal_duplication_rate: Option<f32>,
    /// The normal's structural variant read counts.
    #[serde(deserialize_with = "maybe_sv_info")]
    pub normal_sv_info: Option<SvInfo>,
    /// VarDict's tumor-normal p-value, the smaller one-sided Fisher p-value rounded to 5 decimals.
    pub somatic_p_value: &'a str,
    /// VarDict's tumor-normal odds ratio.
    pub somatic_odds_ratio: &'a str,
}

impl<'a> TumorNormalVariant<'a> {
    /// Return the tumor's half of this record as a one-sample call.
    pub fn tumor_call(&self) -> TumorOnlyVariant<'a> {
        self.sample_call(
            &self.tumor_columns,
            &self.normal_columns,
            self.tumor_duplication_rate,
            &self.tumor_sv_info,
        )
    }

    /// Return the normal's half of this record as a one-sample call.
    pub fn normal_call(&self) -> TumorOnlyVariant<'a> {
        self.sample_call(
            &self.normal_columns,
            &self.tumor_columns,
            self.normal_duplication_rate,
            &self.normal_sv_info,
        )
    }

    /// Return one sample's columns with this record's shared columns as a one-sample call, with the
    /// high-quality ALT reads recovered from its SN; the other columns VarDict prints only in
    /// tumor-only mode keep their defaults.
    fn sample_call(
        &self,
        columns: &SampleColumns<'a>,
        other: &SampleColumns<'a>,
        duplication_rate: Option<f32>,
        sv_info: &Option<SvInfo>,
    ) -> TumorOnlyVariant<'a> {
        TumorOnlyVariant {
            sample: self.sample,
            interval_name: self.interval_name,
            contig: self.contig,
            start: self.start,
            end: self.end,
            ref_allele: self.ref_allele,
            alt_allele: self.alt_allele,
            depth: columns.depth,
            alt_depth: columns.alt_depth,
            ref_forward: columns.ref_forward,
            ref_reverse: columns.ref_reverse,
            alt_forward: columns.alt_forward,
            alt_reverse: columns.alt_reverse,
            gt: columns.gt,
            af: columns.af,
            strand_bias: columns.strand_bias,
            mean_position_in_read: columns.mean_position_in_read,
            stdev_position_in_read: columns.stdev_position_in_read,
            base_quality_mean: columns.base_quality_mean,
            stdev_base_stdev: columns.stdev_base_stdev,
            strand_bias_p_value: columns.strand_bias_p_value,
            strand_bias_odds_ratio: columns.strand_bias_odds_ratio,
            mean_mapping_quality: columns.mean_mapping_quality,
            signal_to_noise: columns.signal_to_noise_thousandths as f32 / 1000.0,
            af_high_quality_bases: columns.af_high_quality_bases,
            high_quality_variant_reads: self.high_quality_variant_reads(columns, other),
            af_adjusted: columns.af_adjusted,
            num_bases_3_prime_shift_for_deletions: self.num_bases_3_prime_shift_for_deletions,
            microsatellite: self.microsatellite,
            microsatellite_length: self.microsatellite_length,
            mean_mismatches_in_reads: columns.mean_mismatches_in_reads,
            flank_seq_5_prime: self.flank_seq_5_prime,
            flank_seq_3_prime: self.flank_seq_3_prime,
            segment: self.segment,
            variant_type: self.variant_type,
            duplication_rate,
            sv_info: sv_info.clone(),
            ..Default::default()
        }
    }

    /// Return a sample's ALT reads at or above VarDict's -q, which VarDict does not print for a
    /// tumor-normal pair, recovered from its ALT depth and SN: 0 without ALT reads, whatever SN
    /// says, and missing when no single count reproduces SN or when SN may be the other sample's.
    /// VarDict copies SN to the sample it fills in by subtraction when it re-calls a long indel on
    /// both samples' pooled reads and labels it Germline, setting that sample's PSTD and QSTD to 1.
    fn high_quality_variant_reads(&self, columns: &SampleColumns, other: &SampleColumns) -> i32 {
        let Ok(alt_depth) = u32::try_from(columns.alt_depth) else {
            return i32::missing();
        };
        if alt_depth == 0 {
            return 0;
        }
        let copied = self.status == "Germline"
            && columns.stdev_position_in_read == 1.0
            && columns.stdev_base_stdev == 1.0
            && columns.signal_to_noise_thousandths == other.signal_to_noise_thousandths;
        if copied {
            return i32::missing();
        }
        hicnt_from_sn(alt_depth, columns.signal_to_noise_thousandths)
            .and_then(|hicnt| i32::try_from(hicnt).ok())
            .unwrap_or(i32::missing())
    }

    /// Return the "TUMOR_NORMAL_FISHER_P" INFO field for this record: the one-sided Fisher exact
    /// p-value that the ALT reads make up more of the tumor's depth than of the normal's.
    pub fn tumor_normal_fisher_p_value(&self) -> f32 {
        let counts = |columns: &SampleColumns| {
            let alt = columns.alt_depth.max(0);
            (alt as u64, (columns.depth - alt).max(0) as u64)
        };
        let (tumor_alt, tumor_other) = counts(&self.tumor_columns);
        let (normal_alt, normal_other) = counts(&self.normal_columns);
        fisher_exact_greater(tumor_alt, tumor_other, normal_alt, normal_other) as f32
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

/// Create a VCF header for VarDict/VarDictJava in tumor-only mode, declaring the applied FILTER labels.
pub fn tumor_only_header(sample: &str, filters: &FilterThresholds) -> Header {
    vcf_header(&[sample], filters)
}

/// Create a VCF header for VarDict/VarDictJava run on a tumor and normal pair, declaring the applied
/// FILTER labels.
pub fn tumor_normal_header(tumor: &str, normal: &str, filters: &FilterThresholds) -> Header {
    vcf_header(&[tumor, normal], filters)
}

/// Create a VCF header for one tumor-only sample or a tumor and normal pair.
#[rustfmt::skip]
fn vcf_header(samples: &[&str], filters: &FilterThresholds) -> Header {
    let tumor_normal = samples.len() == 2;
    let source = [CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-");
    let mut header = Header::default();
    for sample in samples {
        header.push_sample(sample.as_bytes());
    }
    header.remove_filter(b"PASS");
    header.push_record(format!("##source={source}").as_bytes());
    header.push_record(r#"##INFO=<ID=TYPE,Number=A,Type=String,Description="VarDict's class of the change from REF to ALT, from VarDict's own rule applied to the alleles in this record: SNV (one base to one base), Insertion (ALT is the single REF base followed by inserted bases), Deletion (REF is longer and ALT is its first base), Complex (every other change, including MNVs, which VarDict does not separate), or DEL, DUP or INV for a symbolic structural variant. Absent when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=INDEL_3P_SHIFT,Number=A,Type=Integer,Description="Number of bases this insertion or deletion can slide toward the 3' end (rightward on the forward strand) and still describe the same haplotype, counted within VarDict's 70-base window, so values top out near 70; POS plus this is the right-most equivalent position. Written only for insertions and deletions.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=REPEAT_UNIT_COPIES,Number=A,Type=Float,Description="Copies of the 1-6 bp tandem repeat unit next to the variant, counted by VarDict in the reference or ALT haplotype, whichever is larger (1 means no repeat). VarDict takes the unit from one side of the variant (for SNVs, the repeat starting at the following base), so a repeat on the other side can be missed; it also raises the count to an indel's 3' shift divided by its length when that is larger, which can make it fractional and leave REPEAT_UNIT_LEN describing a different repeat. Sequence context, not microsatellite instability. Rounded by VarDict to 3 decimals; absent when VarDict did not compute it.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=REPEAT_UNIT_LEN,Number=A,Type=Integer,Description="Length in bp (1-6) of the repeat unit counted in REPEAT_UNIT_COPIES; absent when VarDict did not compute it.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End position, written only on records with a symbolic ALT allele: the last deleted base for DEL, the last inverted base for INV, and VarDict's end of the duplication for DUP, which can be off by one depending on how VarDict found it.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Signed length of the structural variant, written only on records with a symbolic ALT allele: VarDict's own event length, negative for DEL (the deleted bases), positive for DUP (the duplicated bases) and INV (the inverted bases).">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Structural variant type, written only on records with a symbolic ALT allele: DEL, DUP or INV.">"#.as_bytes());
    if tumor_normal {
        header.push_record(r#"##INFO=<ID=VARDICT_STATUS,Number=A,Type=String,Description="VarDict's tumor-normal label for this allele, from where the allele is found, its AF and VarDict's own call rules (-f, -r, -q and the others); not a statistical test. StrongSomatic: the normal has no reads of the allele; LikelySomatic: normal AF below -V or a single normal read; Germline: the normal's reads of the allele pass VarDict's call rules at a normal AF at or above -V, or at a tumor AF at or above -V when the top-ranked tumor allele at the position failed the call rules; AFDiff: normal AF at or above -V with more than one read, but the normal's reads fail VarDict's call rules; LikelyLOH: tumor AF above 1 minus -V while the normal's reads pass VarDict's call rules at an AF between 0.2 and 0.8, or a tumor AF below -V when the top-ranked tumor allele at the position failed the call rules; StrongLOH: only the normal has reads of the allele, and they pass VarDict's call rules; SampleSpecific: VarDict has no normal record at the position; Deletion: VarDict has no tumor record at the position. VarDict labels an allele both samples carry before it drops up to three low-quality normal reads of it as noise, so the normal's AD and AF can count fewer reads than the label saw, and it relabels an SNV that loses them StrongSomatic. It also re-calls long indels with few reads on both samples' pooled reads, which can relabel a call Germline, even with no normal reads under -p, which lifts its -f and -r limits.">"#.as_bytes());
        header.push_record(r#"##INFO=<ID=TUMOR_NORMAL_FISHER_P,Number=A,Type=Float,Description="One-sided Fisher exact p-value that this allele's read fraction is higher in the tumor than in the normal, from each sample's ALT reads (FORMAT AD[1]) against its other reads (FORMAT DP minus AD[1]), so it compares the two samples' FORMAT AF. Computed by vartovcf, not VarDict's SSF, which is the smaller of the two one-sided p-values rounded to 5 decimals. 1 when either sample has no depth; stored as a 32-bit float, so values below about 1e-45 read 0.">"#.as_bytes());
    }
    header.push_record(r#"##FILTER=<ID=PASS,Description="The variant call has passed all filters and may be considered for downstream analysis.">"#.as_bytes());
    for line in filters.header_lines() {
        header.push_record(line.as_bytes());
    }
    header.push_record(format!(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="The genotype inferred from the allele frequency alone, since VarDict does not genotype: 1/1 when AF >= {MIN_HOM_ALT_AF}, 0/1 for any other call with ALT reads, 0/0 when no read carries the ALT allele or there is no ALT allele, and ./. when DP is 0 too.">"#).as_bytes());
    header.push_record(r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for REF then ALT as VarDict counts them: ALT is the reads carrying this allele and REF is the reads carrying the reference base at the variant's first internal base (the first deleted base for deletions, without insertion-carrying reads for insertions), or at the following base when an insertion at the same position raises DP. REF is missing for Complex calls, where VarDict counts only their first base and can count a read as both, and REF is the only value when there is no ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Reads on the forward strand (SAM flag 0x10 unset) supporting REF then ALT, the forward half of AD as VarDict counts it: REF is counted where AD's is (sometimes the following base), is missing for Complex calls and is the only value when there is no ALT allele. When VarDict runs with -u, overlapping mates count only through the reverse read.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Reads on the reverse strand (SAM flag 0x10 set) supporting REF then ALT, the reverse half of AD as VarDict counts it: REF is counted where AD's is (sometimes the following base), is missing for Complex calls and is the only value when there is no ALT allele. When VarDict runs with -u, overlapping mates count only through the reverse read.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=STRAND_BIAS_FISHER_P,Number=A,Type=Float,Description="Two-sided Fisher exact p-value that the ALT allele's forward/reverse read split differs from REF's, from VarDict's table of REF and ALT reads by SAM strand, not a test against 50:50; its REF counts are VarDict's, including the overcount on Complex calls. Rounded by VarDict to 5 decimals, so values below 0.000005 read 0. Missing when no read carries the ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth as VarDict counts it at the variant's first internal base (POS+1 for deletions), the denominator of FORMAT AF and usually of VarDict's own AF: reads with any base quality, reads whose deletion spans the base, reference-matching soft-clipped bases, reads VarDict reassigned by realignment, and N calls only under -K; overlapping mates count twice unless VarDict ran with -u or -UN, and an insertion at the same position can make it the following base's depth. REF and ALT depths need not sum to DP.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=AF,Number=A,Type=Float,Description="Fraction of reads carrying the ALT allele, AD[1] / DP, computed by vartovcf at full precision and clamped to [0, 1] because VarDict's ALT count can exceed DP; AD and DP keep the raw counts. VarDict's own AF column is rounded to 4 decimals and its denominator can differ from DP, for example at a position shared with an insertion. Missing when there is no ALT allele or DP is 0.">"#.as_bytes());
    let hicnt = "ALT reads whose per-read variant quality, the one QMEAN averages (the base's quality for SNVs, the better flanking base's for deletions and so on), is at least VarDict's -q; mapping quality is not considered. Missing when there is no ALT allele.";
    let hicnt_tumor_normal = if tumor_normal {
        " VarDict does not print this count for a tumor-normal pair, so vartovcf recovers it from AD[1] and VarDict's SN (these reads over the ALT reads below -q, to 3 decimals): 0 when AD[1] is 0, missing unless exactly one count reproduces SN (several can once a thousand or more ALT reads fall below -q), and missing when SN may be the other sample's, as VarDict copies SN to a sample whose counts it fills in by subtraction when it re-calls a long indel on both samples' pooled reads and labels it Germline."
    } else {
        ""
    };
    header.push_record(format!(r#"##FORMAT=<ID=HICNT,Number=A,Type=Integer,Description="{hicnt}{hicnt_tumor_normal}">"#).as_bytes());
    header.push_record(r#"##FORMAT=<ID=REALIGNED_FRAC_OF_DP,Number=A,Type=Float,Description="Fraction of VarDict's depth made of reads it reassigned to this ALT allele by local realignment or MNV merging (VarDict's ExtraAF, which var2vcf_valid.pl calls ADJAF). Those reads are already counted in AD[1] and AF, so this is not an adjusted AF; divided by AF it gives roughly the share of the ALT reads that realignment contributed. Rounded by VarDict to 4 decimals. Missing when no read carries the ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=MEAN_DIST_TO_READ_END,Number=A,Type=Float,Description="Mean, over the reads carrying the ALT allele, of the 1-based distance from the variant to the nearer end of the read's aligned part, soft clips excluded. It is measured from the variant base for SNVs, the first inserted base for insertions and the first base after the gap for deletions; a Complex call is measured like the insertion or deletion it starts with, or from the last base of the block when it starts with a mismatch. Reads VarDict realigned out of soft clips contribute their clip length. Rounded by VarDict to 1 decimal. Missing when no read carries the ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ALT_READ_POS_VARIES,Number=A,Type=Integer,Description="1 when the reads carrying the ALT allele place it at two or more distinct distances from the read end, and 0 when every one has it at the same distance, which with several reads is a common sign of an artifact. VarDict also sets it to 1 whenever it reassigns reads to the allele, as in local realignment or MNV merging, so it is 1 whenever REALIGNED_FRAC_OF_DP is above 0, and for a sample of a tumor-normal pair whose counts it derived by subtracting the other sample from both samples' pooled reads; otherwise a call with one ALT read is always 0. Missing when no read carries the ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=QMEAN,Number=A,Type=Float,Description="Mean, over the reads carrying the ALT allele, of VarDict's per-read variant quality: the base's Phred quality for SNVs, the mean of the inserted bases for insertions, the higher of the two flanking bases for deletions, and the mean of the block's bases for Complex calls, with that flanking base averaged in when the block starts with a deletion. VarDict grows an MNV only through mismatches at least 5 above -q, so MNV values run high. Bases below -q are included and nothing is capped. Rounded by VarDict to 1 decimal. Missing when no read carries the ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=MEAN_MAPQ,Number=A,Type=Float,Description="Arithmetic mean mapping quality of the reads carrying the ALT allele, not the RMS over all reads that the VCF specification's MQ means; uncapped, so a MAPQ of 255 (unavailable) counts as 255. Rounded by VarDict to 1 decimal. Missing when no read carries the ALT allele.">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=MEAN_MISMATCHES,Number=A,Type=Float,Description="Mean, over the reads carrying the ALT allele, of each read's substitution mismatches: its NM tag minus every inserted or deleted base, and for Complex calls minus some of the mismatches VarDict folded into the allele after its first change. Includes the variant's own mismatch for SNVs and excludes soft clips; a read without an NM tag counts as 0, and reads with more than VarDict's -m mismatches (default 8) are not counted. Rounded by VarDict to 1 decimal. Missing when no read carries the ALT allele.">"#.as_bytes());
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
    use rust_htslib::bcf::header::{TagLength, TagType};
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
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<INV>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<INV11063372>", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 58.0, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 33.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "INV", duplication_rate: None, sv_info: Some(inv_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<DUP>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<DUP11063372>", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 58.0, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 10.0, mean_mapping_quality: 33.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "DUP", duplication_rate: Some(0.001), sv_info: Some(dup_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "A", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 1, alt_reverse: 0, gt: "G/A", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 13.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, strand_bias_p_value: 0.34385, strand_bias_odds_ratio: 4.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1.0, microsatellite_length: 1, mean_mismatches_in_reads: 2.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "FLT3", contig: "chr13", start: 24684729, end: 28034141, ref_allele: "G", alt_allele: "<DEL>", depth: 7463, alt_depth: 2, ref_forward: 7463, ref_reverse: 0, alt_forward: 0, alt_reverse: 2, gt: "-3349412/-3349412", af: 0.0003, strand_bias: "0;0", mean_position_in_read: 60.5, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 2.0, mean_mapping_quality: 60.0, signal_to_noise: 4.0, af_high_quality_bases: 1.0000, af_adjusted: 0.0003, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 2, high_quality_total_reads: 2, flank_seq_5_prime: "TGCTGTAGTCTAATGATTCT", flank_seq_3_prime: "CAACGTAGAAGTACTCATTA", segment: "chr13:28033879-28034298", variant_type: "DEL", duplication_rate: None, sv_info: Some(del_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "T", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 0, alt_reverse: 1, gt: "G/T", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 28.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1.0, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713880, end: 114713880, ref_allele: "T", alt_allele: "A", depth: 8211, alt_depth: 1, ref_forward: 3001, ref_reverse: 5130, alt_forward: 1, alt_reverse: 0, gt: "T/A", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 18.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, strand_bias_p_value: 0.36916, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 2.0, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8132, flank_seq_5_prime: "CCTTCGCCTGTCCTCATGTA", flank_seq_3_prime: "TGGTCTCTCATGGCACTGTA", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
        ]
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
    fn test_tumor_only_variant_alt_read_summaries_are_missing_without_alt_reads(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.alt_depth = 0;
        assert!(variant.realigned_frac_of_dp_value().is_missing());
        assert!(variant.mean_dist_to_read_end_value().is_missing());
        assert!(variant.alt_read_pos_varies_value().is_missing());
        assert!(variant.mean_mapq_value().is_missing());
        assert!(variant.strand_bias_fisher_p_value().is_missing());
        assert!(variant.qmean_value().is_missing());
        assert!(variant.mean_mismatches_value().is_missing());
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
    #[case("G", "A", 0, 10, [GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)])]
    #[case("G", "A", 0, 0, [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing])]
    #[case("G", "G", 0, 0, [GenotypeAllele::UnphasedMissing, GenotypeAllele::UnphasedMissing])]
    fn test_tumor_only_variant_gt_value_without_alt_reads(
        variants: Vec<TumorOnlyVariant<'static>>,
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] alt_depth: i32,
        #[case] depth: i32,
        #[case] expected: [GenotypeAllele; 2],
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        variant.ref_allele = ref_allele;
        variant.alt_allele = alt_allele;
        variant.alt_depth = alt_depth;
        variant.depth = depth;
        assert_eq!(variant.gt_value(MIN_HOM_ALT_AF), &expected);
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
        let header = tumor_only_header("dna00001", &FilterThresholds::default());
        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let _ = VcfWriter::from_path(file.path(), &header, true, Format::Vcf).unwrap();
        let reader = VcfReader::from_path(file.path()).expect("Error opening tempfile!");
        let records = reader.header().header_records();
        let samples = reader.header().samples();
        assert_eq!(records.len(), 29);
        assert_eq!(samples.len(), 1);
        assert!(samples.iter().all(|&s| s == "dna00001".as_bytes()));
    }

    fn tumor_normal_records() -> Vec<csv::StringRecord> {
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path("tests/calls.tumor-normal.var")
            .expect("Cannot open the tumor-normal fixture!")
            .records()
            .collect::<Result<_, _>>()
            .expect("Cannot read the tumor-normal fixture!")
    }

    #[test]
    fn test_tumor_normal_variant_reads_both_samples_columns() {
        let records = tumor_normal_records();
        let row: TumorNormalVariant = records[0].deserialize(None).unwrap();
        assert_eq!(
            (
                row.sample,
                row.contig,
                row.start,
                row.ref_allele,
                row.alt_allele
            ),
            ("T|N", "cA", 200, "A", "C")
        );
        assert_eq!(
            (row.tumor_columns.depth, row.tumor_columns.alt_depth),
            (30, 10)
        );
        assert_eq!(
            (row.normal_columns.depth, row.normal_columns.alt_depth),
            (30, 0)
        );
        assert_eq!(row.normal_columns.gt, "A/A");
        assert_eq!((row.status, row.variant_type), ("StrongSomatic", "SNV"));
        assert_eq!((row.microsatellite, row.microsatellite_length), (1.0, 1));
        assert_eq!(row.segment, "cA:150-360");
    }

    #[test]
    fn test_tumor_normal_variant_calls_hold_each_samples_columns() {
        let records = tumor_normal_records();
        let row: TumorNormalVariant = records[0].deserialize(None).unwrap();
        let (tumor, normal) = (row.tumor_call(), row.normal_call());
        assert_eq!(tumor.ad_value(), vec![20, 10]);
        assert_eq!(normal.ad_value(), vec![30, 0]);
        assert_eq!(tumor.adf_value(), vec![10, 5]);
        assert_eq!(tumor.qmean_value(), 30.0);
        assert_eq!(tumor.mean_dist_to_read_end_value(), 22.9);
        assert!(normal.qmean_value().is_missing());
        assert!(normal.mean_dist_to_read_end_value().is_missing());
        assert_eq!(
            normal.gt_value(MIN_HOM_ALT_AF),
            &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)]
        );
        for call in [&tumor, &normal] {
            assert_eq!((call.contig, call.start, call.end), ("cA", 200, 200));
            assert_eq!((call.ref_allele, call.alt_allele), ("A", "C"));
            assert_eq!(call.repeat_unit_copies_value(), Some(1.0));
            assert_eq!(call.variant_class(), Some("SNV"));
        }
    }

    fn hicnt_records() -> Vec<csv::StringRecord> {
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path("tests/calls.tumor-normal.hicnt.var")
            .expect("Cannot open the HICNT fixture!")
            .records()
            .collect::<Result<_, _>>()
            .expect("Cannot read the HICNT fixture!")
    }

    #[rstest]
    #[case::two_counts_give_the_same_sn(0, None, Some(0))]
    #[case::alt_depth_one_above_hicnt_plus_locnt(1, Some(4), Some(0))]
    #[case::sn_copied_to_a_sample_filled_in_by_subtraction(2, None, None)]
    #[case::normal_sn_from_the_reference_reads(3, Some(9), Some(0))]
    fn test_tumor_normal_variant_recovers_hicnt(
        #[case] index: usize,
        #[case] tumor: Option<i32>,
        #[case] normal: Option<i32>,
    ) {
        let records = hicnt_records();
        let row: TumorNormalVariant = records[index].deserialize(None).unwrap();
        let hicnt = |call: TumorOnlyVariant| Some(call.hicnt_value()).filter(|h| !h.is_missing());
        assert_eq!(hicnt(row.tumor_call()), tumor);
        assert_eq!(hicnt(row.normal_call()), normal);
    }

    #[test]
    fn test_tumor_normal_variant_reads_sn_exactly() {
        let records = tumor_normal_records();
        let row: TumorNormalVariant = records[3].deserialize(None).unwrap();
        assert_eq!(row.tumor_columns.signal_to_noise_thousandths, 20_000);
        assert_eq!(row.normal_columns.signal_to_noise_thousandths, 2_000);
        let records = hicnt_records();
        let row: TumorNormalVariant = records[1].deserialize(None).unwrap();
        assert_eq!(row.tumor_columns.signal_to_noise_thousandths, 1_333);
    }

    #[test]
    fn test_tumor_normal_variant_refuses_sn_with_more_than_3_decimals() {
        let records = tumor_normal_records();
        let mut fields: Vec<&str> = records[0].iter().collect();
        fields[21] = "0.3333";
        let record = csv::StringRecord::from(fields);
        let error = record.deserialize::<TumorNormalVariant>(None).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("expected a number with at most 3 decimals, found '0.3333'"),
            "{error}"
        );
    }

    #[rstest]
    #[case(0, 0.0003985065657050687)]
    #[case(8, 1.0)]
    #[case(12, 1.0)]
    #[case(13, 1.0)]
    #[case(16, 0.0006576798232705035)]
    fn test_tumor_normal_variant_fisher_p_value(#[case] index: usize, #[case] expected: f64) {
        let records = tumor_normal_records();
        let row: TumorNormalVariant = records[index].deserialize(None).unwrap();
        assert_eq!(row.tumor_normal_fisher_p_value(), expected as f32);
    }

    #[test]
    fn test_tumor_normal_header() {
        let header = tumor_normal_header("T", "N", &FilterThresholds::default());
        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let _ = VcfWriter::from_path(file.path(), &header, true, Format::Vcf).unwrap();
        let reader = VcfReader::from_path(file.path()).expect("Error opening tempfile!");
        let samples = reader.header().samples();
        assert_eq!(samples, vec![b"T".as_slice(), b"N".as_slice()]);
        assert_eq!(
            reader.header().info_type(b"VARDICT_STATUS").unwrap(),
            (TagType::String, TagLength::AltAlleles)
        );
        assert!(reader.header().info_type(b"TUMOR_NORMAL_FISHER_P").is_ok());
        assert_eq!(
            reader.header().format_type(b"HICNT").unwrap(),
            (TagType::Integer, TagLength::AltAlleles)
        );
        assert!(reader.header().format_type(b"QMEAN").is_ok());
    }
}
