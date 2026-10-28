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

use crate::filter::FilterThresholds;
use crate::fisher::{fisher_exact_greater, fisher_exact_two_sided};
use crate::hicnt::{hicnt_from_sn, parse_thousandths};
use crate::mnv::{LEAK_TAIL, LeakyMnvFix};
use crate::vcf::{GenotypeAllele, Numeric};
use noodles::vcf::Header;
use noodles::vcf::header::FileFormat;
use noodles::vcf::header::record::Value;
use noodles::vcf::header::record::value::Map;
use noodles::vcf::header::record::value::map::{
    AlternativeAllele, Filter, Format, Info, format, info,
};
use serde::{Deserialize, Serialize, de::Error};

const CARGO_PKG_NAME: &str = env!("CARGO_PKG_NAME");
const CARGO_PKG_VERSION: &str = env!("CARGO_PKG_VERSION");

/// The minimum allele frequency at which a variant call is genotyped homozygous alternate.
pub const MIN_HOM_ALT_AF: f32 = 0.8;

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
#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
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
    /// The fraction of reads VarDict's -t removed as duplicates; VarDictJava 1.8.4 of the joachimwolff fork always prints 0,
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

    /// Return the "STRAND_BIAS_FISHER_P" formatted VCF field for this record: the two-sided Fisher
    /// exact p-value of the REF and ALT reads by strand, and missing without ALT reads.
    pub fn strand_bias_fisher_p_value(&self) -> f32 {
        if !self.has_alt_reads() {
            return f32::missing();
        }
        let count = |reads: i32| reads.max(0) as u64;
        fisher_exact_two_sided(
            count(self.ref_forward),
            count(self.ref_reverse),
            count(self.alt_forward),
            count(self.alt_reverse),
        ) as f32
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
        variant_class(self.ref_allele, self.alt_allele)
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

/// Return VarDict's variant class (`Variant.varType`) for a REF and ALT allele: `None` without an ALT
/// allele, SNV, Insertion, Deletion, Complex, or the type of a symbolic structural variant.
pub fn variant_class<'s>(ref_allele: &'s str, alt_allele: &'s str) -> Option<&'s str> {
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
    /// The tumor's duplication rate, which VarDictJava 1.8.4 of the joachimwolff fork always prints as 0.
    #[serde(deserialize_with = "maybe_duplication_rate")]
    pub tumor_duplication_rate: Option<f32>,
    /// The tumor's structural variant read counts.
    #[serde(deserialize_with = "maybe_sv_info")]
    pub tumor_sv_info: Option<SvInfo>,
    /// The normal's duplication rate, which VarDictJava 1.8.4 of the joachimwolff fork always prints as 0.
    #[serde(deserialize_with = "maybe_duplication_rate")]
    pub normal_duplication_rate: Option<f32>,
    /// The normal's structural variant read counts.
    #[serde(deserialize_with = "maybe_sv_info")]
    pub normal_sv_info: Option<SvInfo>,
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
        tumor_normal_fisher_p_value(&self.tumor_call(), &self.normal_call())
    }
}

/// Return the "TUMOR_NORMAL_FISHER_P" INFO field for a tumor and normal call: the one-sided Fisher
/// exact p-value that the ALT reads make up more of the tumor's depth than of the normal's.
pub fn tumor_normal_fisher_p_value(tumor: &TumorOnlyVariant, normal: &TumorOnlyVariant) -> f32 {
    let counts = |call: &TumorOnlyVariant| {
        let alt = call.alt_depth.max(0);
        (alt as u64, (call.depth - alt).max(0) as u64)
    };
    let (tumor_alt, tumor_other) = counts(tumor);
    let (normal_alt, normal_other) = counts(normal);
    fisher_exact_greater(tumor_alt, tumor_other, normal_alt, normal_other) as f32
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

/// The INFO description of the count of partial rows folded into an MNV, which also says how a
/// partial is recognized, by its counts at this leak rate or by its reads, and how each field of
/// such a record is combined. By counts the description is the one written before reads could decide.
fn mnv_partials_description(fix: &LeakyMnvFix) -> String {
    match fix {
        LeakyMnvFix::Counts(leak_rate) => format!(
            r#"Number of VarDict rows folded into this MNV as its partials: SNVs and shorter MNVs that carry its ALT bases over their own span, which VarDict writes for reads that end, have a no-call (N) or a base below its -q plus 5 inside the MNV. A row can be folded only when, in every sample with reads of it, it has fewer ALT reads (c) than the MNV (n) and P(X >= c) is {} or more for X ~ Binomial(n + c, {}); at that rate a partial of 1, 2, 5 or 10 reads needs an MNV of at least {}, {}, {} or {} reads. Written only on records that absorbed partials. On such a record, AD[1], ADF[1] and ADR[1] add the partials' reads, and HICNT adds them too or is missing when any row's is; DP adds the reads of partials starting after the MNV's first base as VarDict wrote it, as VarDict's own MNV merging does, while a partial starting at that base is folded only when its reads fit in DP beside the MNV's, so AD[1] never exceeds DP by more than VarDict's own row did; AF, GT, FILTER, STRAND_BIAS_FISHER_P and TUMOR_NORMAL_FISHER_P follow from these; QMEAN, MEAN_MAPQ, MEAN_MISMATCHES and MEAN_DIST_TO_READ_END are the rows' values weighted by their ALT reads, each as VarDict rounded it and measured at that row's own bases, and the weighted mean is not rounded; ALT_READ_POS_VARIES is 1 when any row's is or the rows' MEAN_DIST_TO_READ_END differ; REALIGNED_FRAC_OF_DP counts the partials' reads as reassigned and is not rounded; VARDICT_STATUS is VarDict's label for the MNV's own row."#,
            LEAK_TAIL,
            leak_rate.value(),
            leak_rate.fewest_mnv_reads(1),
            leak_rate.fewest_mnv_reads(2),
            leak_rate.fewest_mnv_reads(5),
            leak_rate.fewest_mnv_reads(10),
        ),
        LeakyMnvFix::Reads { options, .. } => format!(
            r#"Number of VarDict rows folded into this MNV, whole or in part, as its partials: SNVs and shorter MNVs that carry its ALT bases over their own span, which VarDict writes for reads that end, have a no-call (N) or a base below its -q plus 5 inside the MNV. Each read VarDict counted toward the row was followed in the BAM VarDict read to the MNV's other changed bases, with its mate's bases called into one: a read whose template carries the MNV's ALT base at every one of them at quality {q} or more is the MNV's own molecule and is folded; one whose template carries a base other than the ALT at quality {q} or more at any of them is another allele's and stays; one with a no-call, a base below quality {q} or no coverage at some of them is consistent with leakage and is folded only where the MNV has reads in that sample, no read of the row in any sample is another allele's, and, in a sample after the tumor, the tumor keeps no read of the row. A sample's reads are split between the MNV and the row only when the reads found are exactly the reads VarDict counted; otherwise the row moves whole when every read found would move, and stays whole when not. A row with reads that stay is written with only them, with FORMAT MNV_MOVED_READS, and a row with none is not written. Reads were taken as VarDict takes them: mapping quality {mapq} or more, none of the flags {flags:#x} and no supplementary record, at most {mismatches} mismatches beside indels, and an MNV grown through a changed base only at quality {grow} or more, bridging up to {gap} matching bases. Written only on records that absorbed partials. On such a record, AD[1], ADF[1] and ADR[1] add the partials' reads, and HICNT adds them too or is missing when any row's is; DP adds the reads of partials starting after the MNV's first base as VarDict wrote it, as VarDict's own MNV merging does, while a partial starting at that base is folded only when its reads fit in DP beside the MNV's, so AD[1] never exceeds DP by more than VarDict's own row did; AF, GT, FILTER, STRAND_BIAS_FISHER_P and TUMOR_NORMAL_FISHER_P follow from these; QMEAN, MEAN_MAPQ, MEAN_MISMATCHES and MEAN_DIST_TO_READ_END are the rows' values weighted by their ALT reads, each as VarDict rounded it and measured at that row's own bases, or measured on the reads that moved when a row was split, and the weighted mean is not rounded; ALT_READ_POS_VARIES is 1 when any row's is or the rows' MEAN_DIST_TO_READ_END differ; REALIGNED_FRAC_OF_DP counts the partials' reads as reassigned and is not rounded; VARDICT_STATUS is VarDict's label for the MNV's own row."#,
            q = options.base_quality,
            mapq = options.min_mapping_quality,
            flags = options.exclude_flags,
            mismatches = options.max_mismatches,
            grow = options.base_quality + 5.0,
            gap = options.mnv_gap,
        ),
    }
}

/// The FORMAT description of the ALT reads a partial row lost to the MNV it was folded into in part.
const MNV_MOVED_READS_DESCRIPTION: &str = r#"ALT reads moved from this partial row to the MNV it leaked from (that record's INFO MNV_PARTIALS counts this row), as the reads decided; written only on a partial row that kept some reads. In a sample with reads moved, AD[1], ADF[1], ADR[1] and HICNT count the reads that stayed, and QMEAN, MEAN_MAPQ, MEAN_MISMATCHES, MEAN_DIST_TO_READ_END and ALT_READ_POS_VARIES are measured on them and not rounded; in a sample with none moved they are VarDict's; DP, AD[0] and REALIGNED_FRAC_OF_DP are VarDict's, and AF, GT, FILTER and the Fisher tests follow from the reads that stayed."#;

/// The FORMAT description of the ALT reads the folded partial rows added to each sample.
const MNV_PARTIAL_READS_DESCRIPTION: &str = r#"ALT reads the partial rows folded into this MNV (INFO MNV_PARTIALS) added to AD[1]; written only on records that absorbed partials."#;

/// Create a VCF header for VarDict/VarDictJava in tumor-only mode, declaring the applied FILTER labels,
/// the fields that describe the fold when MNV partials are folded, and the 3-mer context when asked.
pub fn tumor_only_header(
    sample: &str,
    filters: &FilterThresholds,
    fix: Option<&LeakyMnvFix>,
    context_3mer: bool,
) -> Header {
    vcf_header(&[sample], filters, fix, context_3mer)
}

/// Create a VCF header for VarDict/VarDictJava run on a tumor and normal pair, declaring the applied
/// FILTER labels, the fields that describe the fold when MNV partials are folded, and the 3-mer
/// context when asked.
pub fn tumor_normal_header(
    tumor: &str,
    normal: &str,
    filters: &FilterThresholds,
    fix: Option<&LeakyMnvFix>,
    context_3mer: bool,
) -> Header {
    vcf_header(&[tumor, normal], filters, fix, context_3mer)
}

/// Create a VCF header for one tumor-only sample or a tumor and normal pair.
#[rustfmt::skip]
fn vcf_header(samples: &[&str], filters: &FilterThresholds, fix: Option<&LeakyMnvFix>, context_3mer: bool) -> Header {
    let tumor_normal = samples.len() == 2;
    let source = [CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-");
    let mut header = Header::builder()
        .set_file_format(FileFormat::new(4, 2))
        .set_sample_names(samples.iter().map(|sample| sample.to_string()).collect())
        .build();
    let key = "source".parse().expect("source is a header key");
    header.insert(key, Value::from(source)).expect("source is an unstructured header record");
    add_info(&mut header, "TYPE", info::Number::AlternateBases, info::Type::String, r#"VarDict's class of the change from REF to ALT, from VarDict's own rule applied to the alleles in this record: SNV (one base to one base), Insertion (ALT is the single REF base followed by inserted bases), Deletion (REF is longer and ALT is its first base), Complex (every other change, including MNVs, which VarDict does not separate), or DEL, DUP or INV for a symbolic structural variant. Absent when there is no ALT allele."#);
    add_info(&mut header, "INDEL_3P_SHIFT", info::Number::AlternateBases, info::Type::Integer, r#"Number of bases this insertion or deletion can slide toward the 3' end (rightward on the forward strand) and still describe the same haplotype, counted by vartovcf in the reference as VarDict counts its Shift3, but without the window that capped VarDict's at 70 for insertions and 71 for deletions; POS plus this is the right-most equivalent position. Written only for insertions and deletions."#);
    add_info(&mut header, "REPEAT_UNIT_COPIES", info::Number::AlternateBases, info::Type::Float, r#"Copies of the 1-6 bp tandem repeat unit next to the variant, counted by VarDict in the reference or ALT haplotype, whichever is larger (1 means no repeat). VarDict takes the unit from one side of the variant (for SNVs, the repeat starting at the following base), so a repeat on the other side can be missed, and for a row it labels Complex it can take the unit as for an SNV even when the call normalizes to an insertion or deletion; it also raises the count to an indel's 3' shift divided by its length when that is larger, which can make it fractional and leave REPEAT_UNIT_LEN describing a different repeat. Sequence context, not microsatellite instability. Rounded by VarDict to 3 decimals; absent when VarDict did not compute it."#);
    add_info(&mut header, "REPEAT_UNIT_LEN", info::Number::AlternateBases, info::Type::Integer, r#"Length in bp (1-6) of the repeat unit counted in REPEAT_UNIT_COPIES; absent when VarDict did not compute it."#);
    add_info(&mut header, "END", info::Number::Count(1), info::Type::Integer, r#"End position, written only on records with a symbolic ALT allele: the last deleted base for DEL, the last inverted base for INV, and VarDict's end of the duplication for DUP, which can be off by one depending on how VarDict found it."#);
    add_info(&mut header, "SVLEN", info::Number::Count(1), info::Type::Integer, r#"Signed length of the structural variant, written only on records with a symbolic ALT allele: VarDict's own event length, negative for DEL (the deleted bases), positive for DUP (the duplicated bases) and INV (the inverted bases)."#);
    add_info(&mut header, "SVTYPE", info::Number::Count(1), info::Type::String, r#"Structural variant type, written only on records with a symbolic ALT allele: DEL, DUP or INV."#);
    if tumor_normal {
        add_info(&mut header, "VARDICT_STATUS", info::Number::AlternateBases, info::Type::String, r#"VarDict's tumor-normal label for this allele, from where the allele is found, its AF and VarDict's own call rules (-f, -r, -q and the others); not a statistical test. StrongSomatic: the normal has no reads of the allele; LikelySomatic: normal AF below -V or a single normal read; Germline: the normal's reads of the allele pass VarDict's call rules at a normal AF at or above -V, or at a tumor AF at or above -V when the top-ranked tumor allele at the position failed the call rules; AFDiff: normal AF at or above -V with more than one read, but the normal's reads fail VarDict's call rules; LikelyLOH: tumor AF above 1 minus -V while the normal's reads pass VarDict's call rules at an AF between 0.2 and 0.8, or a tumor AF below -V when the top-ranked tumor allele at the position failed the call rules; StrongLOH: only the normal has reads of the allele, and they pass VarDict's call rules; SampleSpecific: VarDict has no normal record at the position; Deletion: VarDict has no tumor record at the position. VarDict labels an allele both samples carry before it drops up to three low-quality normal reads of it as noise, so the normal's AD and AF can count fewer reads than the label saw, and it relabels an SNV that loses them StrongSomatic. It also re-calls long indels with few reads on both samples' pooled reads, which can relabel a call Germline, even with no normal reads under -p, which lifts its -f and -r limits."#);
        add_info(&mut header, "TUMOR_NORMAL_FISHER_P", info::Number::AlternateBases, info::Type::Float, r#"One-sided Fisher exact p-value that this allele's read fraction is higher in the tumor than in the normal, from each sample's ALT reads (FORMAT AD[1]) against its other reads (FORMAT DP minus AD[1]), so it compares the two samples' FORMAT AF. Computed by vartovcf, not VarDict's SSF, which is the smaller of the two one-sided p-values rounded to 5 decimals. 1 when either sample has no depth; stored as a 32-bit float, so values below about 1e-45 read 0."#);
    }
    if let Some(fix) = fix {
        add_info(&mut header, "MNV_PARTIALS", info::Number::Count(1), info::Type::Integer, &mnv_partials_description(fix));
    }
    if context_3mer {
        add_info(&mut header, "CONTEXT_3mer", info::Number::Count(1), info::Type::String, r#"The reference bases at POS-1, POS and POS+1, uppercase even where the reference is soft-masked, with N for a base past either end of the contig. POS is the record's own after normalization, so the middle base of an insertion or deletion is its anchor base. Not reverse-complemented to a pyrimidine middle base, and an IUPAC ambiguity code is kept as the reference spells it, even where REF writes the first base it stands for."#);
    }
    add_filter(&mut header, "PASS", r#"The call failed none of the FILTER labels, which test the tumor (or only) sample; a call whose tumor sample has no reads of the ALT allele keeps no FILTER (.)."#);
    for (id, description) in filters.header_filters() {
        add_filter(&mut header, id, &description);
    }
    add_format(&mut header, "GT", format::Number::Count(1), format::Type::String, &format!(r#"The genotype inferred from the allele frequency alone, since VarDict does not genotype: 1/1 when AF >= {MIN_HOM_ALT_AF}, 0/1 for any other call with ALT reads, 0/0 when no read carries the ALT allele or there is no ALT allele, and ./. when DP is 0 too."#));
    add_format(&mut header, "AD", format::Number::ReferenceAlternateBases, format::Type::Integer, r#"Read depth for REF then ALT as VarDict counts them: ALT is the reads carrying this allele and REF is the reads carrying the reference base at the variant's first internal base (the first deleted base for deletions, without insertion-carrying reads for insertions), or at the following base when an insertion at the same position raises DP. REF is missing for calls VarDict labels Complex, MNVs included, whatever their TYPE after trimming, as VarDict counts only their first base and can count a read as both, and REF is the only value when there is no ALT allele."#);
    add_format(&mut header, "ADF", format::Number::ReferenceAlternateBases, format::Type::Integer, r#"Reads on the forward strand (SAM flag 0x10 unset) supporting REF then ALT, the forward half of AD as VarDict counts it: REF is counted where AD's is (sometimes the following base), is missing where AD's is and is the only value when there is no ALT allele. When VarDict runs with -u, overlapping mates count only through the reverse read."#);
    add_format(&mut header, "ADR", format::Number::ReferenceAlternateBases, format::Type::Integer, r#"Reads on the reverse strand (SAM flag 0x10 set) supporting REF then ALT, the reverse half of AD as VarDict counts it: REF is counted where AD's is (sometimes the following base), is missing where AD's is and is the only value when there is no ALT allele. When VarDict runs with -u, overlapping mates count only through the reverse read."#);
    add_format(&mut header, "STRAND_BIAS_FISHER_P", format::Number::AlternateBases, format::Type::Float, r#"Two-sided Fisher exact p-value that the ALT allele's forward/reverse read split differs from REF's, from VarDict's table of REF and ALT reads by SAM strand, not a test against 50:50; its REF counts are VarDict's, including the overcount on Complex calls. Computed by vartovcf at full precision rather than read from VarDictJava's -fisher column, which rounds it to 5 decimals; stored as a 32-bit float, so values below about 1e-45 read 0. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "DP", format::Number::Count(1), format::Type::Integer, r#"Read depth as VarDict counts it at the variant's first internal base where VarDict placed the variant, which is POS+1 for a deletion that normalization did not move left, the denominator of FORMAT AF and usually of VarDict's own AF: reads with any base quality, reads whose deletion spans the base, reference-matching soft-clipped bases, reads VarDict reassigned by realignment, and N calls only under -K; overlapping mates count twice unless VarDict ran with -u or -UN, and an insertion at the same position can make it the following base's depth. REF and ALT depths need not sum to DP."#);
    add_format(&mut header, "AF", format::Number::AlternateBases, format::Type::Float, r#"Fraction of reads carrying the ALT allele, AD[1] / DP, computed by vartovcf at full precision and clamped to [0, 1] because VarDict's ALT count can exceed DP; AD and DP keep the raw counts. VarDict's own AF column is rounded to 4 decimals and its denominator can differ from DP, for example at a position shared with an insertion. Missing when there is no ALT allele or DP is 0."#);
    let hicnt = "ALT reads whose per-read variant quality, the one QMEAN averages (the base's quality for SNVs, the better flanking base's for deletions and so on), is at least VarDict's -q; mapping quality is not considered. Missing when there is no ALT allele.";
    let hicnt_tumor_normal = if tumor_normal {
        " VarDict does not print this count for a tumor-normal pair, so vartovcf recovers it from AD[1] and VarDict's SN (these reads over the ALT reads below -q, to 3 decimals): 0 when AD[1] is 0, missing unless exactly one count reproduces SN (several can once a thousand or more ALT reads fall below -q), and missing when SN may be the other sample's, as VarDict copies SN to a sample whose counts it fills in by subtraction when it re-calls a long indel on both samples' pooled reads and labels it Germline."
    } else {
        ""
    };
    add_format(&mut header, "HICNT", format::Number::AlternateBases, format::Type::Integer, &format!("{hicnt}{hicnt_tumor_normal}"));
    add_format(&mut header, "REALIGNED_FRAC_OF_DP", format::Number::AlternateBases, format::Type::Float, r#"Fraction of VarDict's depth made of reads it reassigned to this ALT allele by local realignment or MNV merging (VarDict's ExtraAF, which var2vcf_valid.pl calls ADJAF). Those reads are already counted in AD[1] and AF, so this is not an adjusted AF; divided by AF it gives roughly the share of the ALT reads that realignment contributed. Rounded by VarDict to 4 decimals, except on a record with INFO MNV_PARTIALS, where it is recomputed and not rounded. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "MEAN_DIST_TO_READ_END", format::Number::AlternateBases, format::Type::Float, r#"Mean, over the reads carrying the ALT allele, of the 1-based distance from the variant to the nearer end of the read's aligned part, soft clips excluded. It is measured from the variant base for SNVs, the first inserted base for insertions and the first base after the gap for deletions; a Complex call is measured like the insertion or deletion it starts with, or from the last base of the block when it starts with a mismatch. Reads VarDict realigned out of soft clips contribute their clip length. Rounded by VarDict to 1 decimal, except on a record with INFO MNV_PARTIALS, whose weighted mean is not rounded. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "ALT_READ_POS_VARIES", format::Number::AlternateBases, format::Type::Integer, r#"1 when the reads carrying the ALT allele place it at two or more distinct distances from the read end, and 0 when every one has it at the same distance, which with several reads is a common sign of an artifact. VarDict also sets it to 1 whenever it reassigns reads to the allele, as in local realignment or MNV merging, so it is 1 whenever REALIGNED_FRAC_OF_DP is above 0 on a record without INFO MNV_PARTIALS, while a folded record is 1 only when a folded row's is or the rows' MEAN_DIST_TO_READ_END differ, and for a sample of a tumor-normal pair whose counts it derived by subtracting the other sample from both samples' pooled reads; otherwise a call with one ALT read is always 0. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "QMEAN", format::Number::AlternateBases, format::Type::Float, r#"Mean, over the reads carrying the ALT allele, of VarDict's per-read variant quality: the base's Phred quality for SNVs, the mean of the inserted bases for insertions, the higher of the two flanking bases for deletions, and the mean of the block's bases for Complex calls, with that flanking base averaged in when the block starts with a deletion. VarDict grows an MNV only through mismatches at least 5 above -q, so MNV values run high. Bases below -q are included and nothing is capped. Rounded by VarDict to 1 decimal, except on a record with INFO MNV_PARTIALS, whose weighted mean is not rounded. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "MEAN_MAPQ", format::Number::AlternateBases, format::Type::Float, r#"Arithmetic mean mapping quality of the reads carrying the ALT allele, not the RMS over all reads that the VCF specification's MQ means; uncapped, so a MAPQ of 255 (unavailable) counts as 255. Rounded by VarDict to 1 decimal, except on a record with INFO MNV_PARTIALS, whose weighted mean is not rounded. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "MEAN_MISMATCHES", format::Number::AlternateBases, format::Type::Float, r#"Mean, over the reads carrying the ALT allele, of each read's substitution mismatches: its NM tag minus every inserted or deleted base, and for Complex calls minus some of the mismatches VarDict folded into the allele after its first change. Includes the variant's own mismatch for SNVs and excludes soft clips; a read without an NM tag counts as 0, and reads with more than VarDict's -m mismatches (default 8) are not counted. Rounded by VarDict to 1 decimal, except on a record with INFO MNV_PARTIALS, whose weighted mean is not rounded. Missing when no read carries the ALT allele."#);
    add_format(&mut header, "SV_SOFTCLIP_READS", format::Number::AlternateBases, format::Type::Integer, r#"Reads soft-clipped at this structural variant's breakpoint whose clipped sequence VarDict matched to the other side; primary alignments only, as VarDict ignores supplementary (SA) records. Written only on records with a symbolic ALT allele."#);
    add_format(&mut header, "SV_DISCORDANT_READS", format::Number::AlternateBases, format::Type::Integer, r#"Discordant read records (unexpected insert size or orientation) in the clusters VarDict linked to this structural variant; each mate counts, so one pair can count twice. Written only on records with a symbolic ALT allele."#);
    if fix.is_some() {
        add_format(&mut header, "MNV_PARTIAL_READS", format::Number::AlternateBases, format::Type::Integer, MNV_PARTIAL_READS_DESCRIPTION);
    }
    if matches!(fix, Some(LeakyMnvFix::Reads { .. })) {
        add_format(&mut header, "MNV_MOVED_READS", format::Number::AlternateBases, format::Type::Integer, MNV_MOVED_READS_DESCRIPTION);
    }
    add_alternative_allele(&mut header, "DEL", r#"Deletion relative to the reference."#);
    add_alternative_allele(&mut header, "DUP", r#"Region of elevated copy number relative to the reference."#);
    add_alternative_allele(&mut header, "INV", r#"Inversion of reference sequence."#);
    header
}

/// Declares an INFO field.
fn add_info(
    header: &mut Header,
    id: &str,
    number: info::Number,
    ty: info::Type,
    description: &str,
) {
    let info = Map::<Info>::new(number, ty, description);
    header.infos_mut().insert(id.to_string(), info);
}

/// Declares a FORMAT field.
fn add_format(
    header: &mut Header,
    id: &str,
    number: format::Number,
    ty: format::Type,
    description: &str,
) {
    let format = Map::<Format>::new(number, ty, description);
    header.formats_mut().insert(id.to_string(), format);
}

/// Declares a FILTER label.
fn add_filter(header: &mut Header, id: &str, description: &str) {
    let filter = Map::<Filter>::new(description);
    header.filters_mut().insert(id.to_string(), filter);
}

/// Declares a symbolic ALT allele.
fn add_alternative_allele(header: &mut Header, id: &str, description: &str) {
    let alternative_allele = Map::<AlternativeAllele>::new(description);
    header
        .alternative_alleles_mut()
        .insert(id.to_string(), alternative_allele);
}

#[cfg(test)]
mod tests {
    use noodles::vcf;
    use noodles::vcf::header::record::value::map::{format, info};
    use pretty_assertions::assert_eq;
    use rstest::*;
    use tempfile::NamedTempFile;

    use super::*;

    /// Writes a header to a file and reads it back, with its text.
    fn read_back(header: &Header, file: &NamedTempFile) -> (vcf::Header, String) {
        std::fs::write(file.path(), crate::vcf::header_bytes(header).unwrap()).unwrap();
        let header = vcf::io::reader::Builder::default()
            .build_from_path(file.path())
            .expect("Error opening tempfile!")
            .read_header()
            .unwrap();
        (header, std::fs::read_to_string(file.path()).unwrap())
    }

    #[fixture]
    #[rustfmt::skip]
    fn variants() -> Vec<TumorOnlyVariant<'static>> {
        let inv_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        let dup_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        let del_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        vec![
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<INV>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<INV11063372>", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 58.0, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, mean_mapping_quality: 33.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "INV", duplication_rate: None, sv_info: Some(inv_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<DUP>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<DUP11063372>", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 58.0, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, mean_mapping_quality: 33.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "DUP", duplication_rate: Some(0.001), sv_info: Some(dup_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "A", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 1, alt_reverse: 0, gt: "G/A", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 13.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1.0, microsatellite_length: 1, mean_mismatches_in_reads: 2.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "FLT3", contig: "chr13", start: 24684729, end: 28034141, ref_allele: "G", alt_allele: "<DEL>", depth: 7463, alt_depth: 2, ref_forward: 7463, ref_reverse: 0, alt_forward: 0, alt_reverse: 2, gt: "-3349412/-3349412", af: 0.0003, strand_bias: "0;0", mean_position_in_read: 60.5, stdev_position_in_read: 1.0, base_quality_mean: 90.0, stdev_base_stdev: 1.0, mean_mapping_quality: 60.0, signal_to_noise: 4.0, af_high_quality_bases: 1.0000, af_adjusted: 0.0003, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0.0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 2, high_quality_total_reads: 2, flank_seq_5_prime: "TGCTGTAGTCTAATGATTCT", flank_seq_3_prime: "CAACGTAGAAGTACTCATTA", segment: "chr13:28033879-28034298", variant_type: "DEL", duplication_rate: None, sv_info: Some(del_info), distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "T", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 0, alt_reverse: 1, gt: "G/T", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 28.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1.0, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713880, end: 114713880, ref_allele: "T", alt_allele: "A", depth: 8211, alt_depth: 1, ref_forward: 3001, ref_reverse: 5130, alt_forward: 1, alt_reverse: 0, gt: "T/A", af: 0.0001, strand_bias: "2;0", mean_position_in_read: 18.0, stdev_position_in_read: 0.0, base_quality_mean: 90.0, stdev_base_stdev: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2.0, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 2.0, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8132, flank_seq_5_prime: "CCTTCGCCTGTCCTCATGTA", flank_seq_3_prime: "TGGTCTCTCATGGCACTGTA", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
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
        let variant = variants.into_iter().nth(2).unwrap();
        let p = variant.strand_bias_fisher_p_value() as f64;
        assert!((p - 0.343_854_852_740_151_5).abs() < 1e-7, "{p}");
    }

    #[rstest]
    fn test_tumor_only_variant_strand_bias_fisher_p_value_is_computed_from_the_strand_counts(
        variants: Vec<TumorOnlyVariant<'static>>,
    ) {
        let mut variant = variants.into_iter().nth(2).unwrap();
        (variant.ref_forward, variant.ref_reverse) = (500, 500);
        (variant.alt_forward, variant.alt_reverse) = (2, 60);
        let p = variant.strand_bias_fisher_p_value() as f64;
        assert!((p - 3.328_906_454_641_75e-15).abs() < 1e-21, "{p}");
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

    #[test]
    fn test_header_is_written_as_vcf_4_2_with_the_source_last() {
        let header = tumor_only_header("dna00001", &FilterThresholds::default(), None, false);
        let text = String::from_utf8(crate::vcf::header_bytes(&header).unwrap()).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "##fileformat=VCFv4.2");
        assert_eq!(
            lines[lines.len() - 2],
            format!("##source={CARGO_PKG_NAME}-{CARGO_PKG_VERSION}")
        );
    }

    #[test]
    fn test_tumor_only_header() {
        let header = tumor_only_header("dna00001", &FilterThresholds::default(), None, false);
        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let (header, text) = read_back(&header, &file);
        let records: Vec<&str> = text.lines().filter(|line| line.starts_with("##")).collect();
        let samples: Vec<&[u8]> = header.sample_names().iter().map(|s| s.as_bytes()).collect();
        assert_eq!(records.len(), 29);
        assert_eq!(samples.len(), 1);
        assert!(samples.iter().all(|&s| s == "dna00001".as_bytes()));
    }

    /// Reads the rows of a tumor-normal fixture written with `-fisher`, without its Fisher columns.
    fn paired_records(path: &str) -> Vec<csv::StringRecord> {
        let columns = crate::Columns::new(crate::Layout::TumorNormal, true);
        csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(path)
            .expect("Cannot open the tumor-normal fixture!")
            .records()
            .map(|record| {
                let mut record = record.expect("Cannot read the tumor-normal fixture!");
                columns.drop_fisher_columns(&mut record, &mut csv::StringRecord::new());
                record
            })
            .collect()
    }

    fn tumor_normal_records() -> Vec<csv::StringRecord> {
        paired_records("tests/calls.tumor-normal.var")
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
        paired_records("tests/calls.tumor-normal.hicnt.var")
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
        let header = tumor_normal_header("T", "N", &FilterThresholds::default(), None, false);
        let file = NamedTempFile::new().expect("Cannot create temporary file!");
        let (header, _) = read_back(&header, &file);
        let samples: Vec<&[u8]> = header.sample_names().iter().map(|s| s.as_bytes()).collect();
        assert_eq!(samples, vec![b"T".as_slice(), b"N".as_slice()]);
        let info_type = |id: &str| {
            header
                .infos()
                .get(id)
                .map(|info| (info.ty(), info.number()))
        };
        let format_type = |id: &str| {
            header
                .formats()
                .get(id)
                .map(|format| (format.ty(), format.number()))
        };
        assert_eq!(
            info_type("VARDICT_STATUS").unwrap(),
            (info::Type::String, info::Number::AlternateBases)
        );
        assert!(info_type("TUMOR_NORMAL_FISHER_P").is_some());
        assert_eq!(
            format_type("HICNT").unwrap(),
            (format::Type::Integer, format::Number::AlternateBases)
        );
        assert!(format_type("QMEAN").is_some());
    }

    #[rstest]
    #[case(false, None)]
    #[case(true, Some((info::Type::String, info::Number::Count(1))))]
    fn test_header_declares_the_3mer_context_only_when_asked(
        #[case] context_3mer: bool,
        #[case] expected: Option<(info::Type, info::Number)>,
    ) {
        let filters = FilterThresholds::default();
        for header in [
            tumor_only_header("dna00001", &filters, None, context_3mer),
            tumor_normal_header("T", "N", &filters, None, context_3mer),
        ] {
            let declared = header
                .infos()
                .get("CONTEXT_3mer")
                .map(|info| (info.ty(), info.number()));
            assert_eq!(declared, expected);
        }
    }
}
