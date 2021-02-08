//! A module for serialization-deserialization friendly VarDict/VarDictJava data types.
use std::clone::Clone;
use std::cmp::PartialEq;
use std::default::Default;
use std::error;
use std::fmt;
use std::ops::Range;
use std::str::FromStr;
use std::string::ToString;

use anyhow::Result;
use bio_types::genome::{AbstractInterval, Position};
use rust_htslib::bcf::Header;
use serde::{de::Error, Deserialize, Serialize};
use strum::EnumString;

const CARGO_PKG_NAME: &str = env!("CARGO_PKG_NAME");
const CARGO_PKG_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Deserialize a possibly infinite float into a <f32> or return a custom error.
///
/// The following cases are handled:
///
/// * `"Inf"`: floating point infinity
/// * `"-Inf"`: floating point negative infinity
/// * `<other>`: a non-infinite floating point number
fn maybe_infinite_f32<'de, D>(deserializer: D) -> Result<f32, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)
        .expect("Could not parse a maybe-infinite (Inf, -Inf) floating point number.");
    if s == "Inf" {
        Ok(f32::INFINITY)
    } else if s == "-Inf" {
        Ok(f32::NEG_INFINITY)
    } else {
        f32::from_str(&s).map_err(D::Error::custom)
    }
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
    let s: &str = Deserialize::deserialize(deserializer)
        .expect("Could not parse an encoded structural variant info.");
    if s == "0" {
        Ok(None)
    } else {
        Ok(Some(SvInfo::from_str(s).map_err(D::Error::custom)?))
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
    let s: &str = Deserialize::deserialize(deserializer)
        .expect("Could not parse an encoded structural variant info.");
    if s == "0" {
        Ok(None)
    } else {
        Ok(Some(f32::from_str(s).map_err(D::Error::custom)?))
    }
}

/// Deserialize VarDict/VarDictJava strand bias status into an enumeration of strand bias. The
/// incoming value takes the values [0-2];[0-2] (_e.g._ "0;2", "2;1"). The first value refers to
/// reads that support the reference allele, and the second to reads that support the variant
/// allele.
///
/// * `0`: there were too few reads to say otherwise (less than 12 for the sum of forward and reverse reads)
/// * `1`: strand bias was detected
/// * `2`: no strand bias was undetected
fn to_strand_bias_status<'de, D>(deserializer: D) -> Result<PairBias, D::Error>
where
    D: serde::Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)
        .expect("Could not parse the paired strand bias status.");
    PairBias::from_str(s).map_err(D::Error::custom)
}

/// An exception for when we cannot parse a string into a `PairBias`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParsePairBiasError;

impl error::Error for ParsePairBiasError {}

impl fmt::Display for ParsePairBiasError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ParsePairBiasError")
    }
}

/// An exception for when we cannot parse a string into a `SvInfo`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseSvInfoError;

impl error::Error for ParseSvInfoError {}

impl fmt::Display for ParseSvInfoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ParseSvInfoError")
    }
}

/// Enumeration of VarDict/VarDictJava strand bias statuses.
#[derive(Debug, Deserialize, EnumString, Eq, PartialEq, Serialize)]
pub enum StrandBias {
    /// There were too few reads to say otherwise (less than 12 for the sum of forward and reverse reads).
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
        write!(f, "{:?}", self)
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

impl ToString for PairBias {
    fn to_string(&self) -> String {
        format!("{}:{}", self.reference, self.alternate)
    }
}

impl FromStr for PairBias {
    type Err = ParsePairBiasError;

    /// Convert a string to a `StrandBias` status.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let items: Vec<&str> = s.split(';').collect();
        let pair = match (items.get(0), items.get(1)) {
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

impl ToString for SvInfo {
    fn to_string(&self) -> String {
        format!(
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
        let sv_info = match (items.get(0), items.get(1), items.get(2)) {
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
#[derive(Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct TumorOnlyVariant<'a> {
    /// Sample name (with whitespace translated to underscores).
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
    /// The number of alternate reads supporting the alternate call.
    pub alt_reverse: i32,
    /// The call's genotype.
    pub gt: &'a str,
    /// The allele frequency of the alternate allele.
    pub af: f32,
    /// Strand bias status. That will take the values [0-2];[0-2] (_e.g._ "0;2", "2;1"). The first
    /// value refers to reads that support the reference allele, and the second to reads that
    /// support the variant allele.
    ///
    /// * `0`: there were too few reads to say otherwise (less than 12 for the sum of forward and reverse reads)
    /// * `1`: strand bias was detected
    /// * `2`: no strand bias was undetected
    #[serde(deserialize_with = "to_strand_bias_status")]
    pub strand_bias: PairBias,
    /// The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads
    /// that support the variant call.
    pub mean_position_in_read: f32,
    /// The standard deviation of the distance to the nearest 5 or 3 prime read end (whichever is
    /// closer) in all reads that support the variant call.
    pub stdev_position_in_read: f32,
    /// The mean base quality (Phred) of all bases that directly support the variant call.
    pub mean_base_quality: f32,
    /// The standard deviation of the base quality (Phred)) of all bases that directly support
    /// the variant call.
    pub stdev_base_quality: f32,
    /// The Fisher test p-value for if you should reject the hypothesis that there is no strand
    /// bias. Non-multiple hypothesis test corrected.
    pub strand_bias_p_value: f32,
    #[serde(deserialize_with = "maybe_infinite_f32")]
    /// The odds ratio for strand bias.
    pub strand_bias_odds_ratio: f32,
    /// The mean mapping quality (Phred) of all reads that directly support the variant call.
    pub mean_mapping_quality: f32,
    /// The signal to noise ratio.
    pub signal_to_noise: i32,
    /// Allele frequency calculated using only high quality bases. Lossy due to rounding.
    pub af_high_quality_bases: f32,
    /// Adjusted allele frequency for indels due to local realignment. Lossy due to rounding.
    pub af_adjusted: f32,
    /// The number of bases to be shifted 3-prime for deletions due to alternative alignment(s).
    pub num_bases_3_prime_shift_for_deletions: i32,
    /// Whether the variant call is in a microsatellite (MSI) or not. Greater than 1 indicates MSI.
    pub microsatellite: i32,
    /// The length of the microSatellite in base pairs of reference genome.
    pub microsatellite_length: i32,
    /// The length of the microSatellite in base pairs of reference genome.
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
    /// The duplication rate, if this call is a duplication.
    #[serde(default, deserialize_with = "maybe_duplication_rate")]
    pub duplication_rate: Option<f32>,
    /// The details of the structural variant as a 0 or as a triplet of ints separated with "=".
    #[serde(default, deserialize_with = "maybe_sv_info")]
    pub sv_info: Option<SvInfo>,
    #[serde(default)]
    /// The distance in reference genome base pairs to the nearest CRISPR-site (CRISPR-mode only).
    pub distance_to_crispr_site: Option<i32>,
}

impl<'a> TumorOnlyVariant<'a> {
    /// Return the "AD" formatted VCF field for this record.
    pub fn ad_value(&self) -> String {
        if self.alt_depth == 0 {
            format!("{}", self.depth)
        } else {
            format!("{},{}", self.depth - self.alt_depth, self.alt_depth)
        }
    }

    /// Return the "ALD" formatted VCF field for this record.
    pub fn ald_value(&self) -> String {
        format!("{},{}", self.alt_forward, self.alt_reverse)
    }

    /// The length of this variant from the perspective of the reference coordinate system.
    pub fn length(&self) -> i32 {
        if self.variant_type == "Deletion" || self.variant_type == "DEL" {
            -((self.end - self.start) as i32)
        } else {
            (self.end - self.start) as i32
        }
    }

    /// Return the "RD" formatted VCF field for this record.
    pub fn rd_value(&self) -> String {
        format!("{},{}", self.ref_forward, self.ref_reverse)
    }
}

impl<'a> AbstractInterval for TumorOnlyVariant<'a> {
    fn contig(&self) -> &str {
        &self.contig
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
    header.push_record(format!("##source={}", source).as_bytes());
    header.push_record(r#"##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand bias flags (UnDetected, Detected, TooFewReads) in the format `reference`:`alternate`.">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Strand bias in reads that support the reference call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Strand bias in reads that support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=PMEAN,Number=1,Type=Float,Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=PSTD,Number=1,Type=Float,Description="The standard deviation of the distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=QUAL,Number=1,Type=Float,Description="The mean base quality (Phred) of all bases that directly support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=QSTD,Number=1,Type=Float,Description="The standard deviation of the base quality (Phred)) of all bases that directly support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SBF,Number=1,Type=Float,Description="The Fisher test p-value for if you should reject the hypothesis that there is no strand bias. Non-multiple hypothesis test corrected">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="The odds ratio for strand bias for this variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=MQ,Number=1,Type=Float,Description="The mean mapping quality (Phred) of all reads that directly support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SN,Number=1,Type=Float,Description="The signal to noise ratio for this variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=HIAF,Number=1,Type=Float,Description="Allele frequency calculated using only high quality bases. Lossy due to rounding">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=ADJAF,Number=1,Type=Float,Description="Adjusted allele frequency for indels due to local realignment. Lossy due to rounding">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description="The number of bases to be shifted 3-prime for deletions due to alternative alignment(s)">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=MSI,Number=1,Type=Float,Description="Whether the variant call is in a microsatellite (MSI) or not. Greater than 1 indicates MSI">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=MSILEN,Number=1,Type=Float,Description="The length of the microSatellite in base pairs of reference genome">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=NM,Number=1,Type=Float,Description="The mean mismatches within all reads that directly support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=LSEQ,Number=1,Type=String,Description="5-prime reference flanking sequence">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=RSEQ,Number=1,Type=String,Description="3-prime reference flanking sequence">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=HICNT,Number=1,Type=Integer,Description="The number of high quality reads supporting the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=HICOV,Number=1,Type=Integer,Description="The number of high quality reads at the locus of the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description="The number of split reads supporting the variant call if this call is a structural variant">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description="The number of paired-end reads supporting the variant call if this call is a structural variant">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The structural variant type (DEL DUP FUS INS INV), if this call is a structural variant">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="The length of stuctural variant in base pairs of reference genome, if this call is a structural variant">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="The duplication rate, if this call is a duplication">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=PASS,Description="The variant call has passed all filters and may be considered for downstream analysis">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=q22.5,Description="The mean base quality (Phred) of all bases that directly support this variant call is below 22.5">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=Q10,Description="The mean mapping quality (Phred) in reads that suppr 10">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=MSI12,Description="The variant call is in a microsatellite region with 12 non-monomer MSI or 13 monomer MSI">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=NM8.0,Description="The mean mismatches in reads that support the variant call is >= 8.0, and might be a false positive or contamination">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=InGap,Description="The variant call is in a deletion gap, and might be a false positive">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=InIns,Description="The variant call was found to be adjacent to an insertion variant">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=Cluster0bp,Description="At least two variant calls are within 0 base pairs from each other in the reference sequence coordinate system">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=LongMSI,Description="The variant call is flanked by a long A/T stretch (>=14 base pairs)">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="The genotype for this sample for this variant call ">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=VD,Number=1,Type=Integer,Description="The variant allele depth at this location">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="The total allele depth at this location which potentially includes No-calls">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="The allelic depths for the REF and ALT alleles">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="The number of variant call forward and reverse reads">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=RD,Number=2,Type=Integer,Description="The number of reference forward and reverse reads">"#.as_bytes());
    header
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use anyhow::Result;
    use csv::ReaderBuilder;
    use rstest::*;
    use rust_htslib::bcf::{Format, HeaderRecord, Read};
    use rust_htslib::bcf::{Reader as VcfReader, Writer as VcfWriter};
    use tempfile::NamedTempFile;

    use super::*;

    #[fixture]
    #[rustfmt::skip]
    fn variants() -> Vec<TumorOnlyVariant<'static>> {
        let sv_info = SvInfo { supporting_split_reads: 1, supporting_pairs: 1, supporting_clusters: 1 };
        vec!(
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "A", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 1, alt_reverse: 0, gt: "G/A", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 13.0, stdev_position_in_read: 0.0, mean_base_quality: 90.0, stdev_base_quality: 0.0, strand_bias_p_value: 0.34385, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 60.0, signal_to_noise: 2, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1, microsatellite_length: 1, mean_mismatches_in_reads: 2.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713883, end: 114713883, ref_allele: "G", alt_allele: "T", depth: 8104, alt_depth: 1, ref_forward: 2766, ref_reverse: 5280, alt_forward: 0, alt_reverse: 1, gt: "G/T", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 28.0, stdev_position_in_read: 0.0, mean_base_quality: 90.0, stdev_base_quality: 0.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: f32::INFINITY, mean_mapping_quality: 60.0, signal_to_noise: 2, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 1, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8048, flank_seq_5_prime: "TCGCCTGTCCTCATGTATTG", flank_seq_3_prime: "TCTCTCATGGCACTGTACTC", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "NRAS-Q61", contig: "chr1", start: 114713880, end: 114713880, ref_allele: "T", alt_allele: "A", depth: 8211, alt_depth: 1, ref_forward: 3001, ref_reverse: 5130, alt_forward: 1, alt_reverse: 0, gt: "T/A", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 18.0, stdev_position_in_read: 0.0, mean_base_quality: 90.0, stdev_base_quality: 0.0, strand_bias_p_value: 0.36916, strand_bias_odds_ratio: f32::NEG_INFINITY, mean_mapping_quality: 60.0, signal_to_noise: 2, af_high_quality_bases: 0.0001, af_adjusted: 0.0, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 2, microsatellite_length: 1, mean_mismatches_in_reads: 1.0, high_quality_variant_reads: 1, high_quality_total_reads: 8132, flank_seq_5_prime: "CCTTCGCCTGTCCTCATGTA", flank_seq_3_prime: "TGGTCTCTCATGGCACTGTA", segment: "chr1:114713749-114713988", variant_type: "SNV", duplication_rate: None, sv_info: None, distance_to_crispr_site: None },
            TumorOnlyVariant { sample: "dna00001", interval_name: "PTPN11", contig: "chr12", start: 112450447, end: 123513818, ref_allele: "A", alt_allele: "<INV>", depth: 6775, alt_depth: 1, ref_forward: 3991, ref_reverse: 2588, alt_forward: 1, alt_reverse: 0, gt: "A/<INV11063372>", af: 0.0001, strand_bias: PairBias::from_str("2;0").unwrap(), mean_position_in_read: 58.0, stdev_position_in_read: 1.0, mean_base_quality: 90.0, stdev_base_quality: 1.0, strand_bias_p_value: 1.0, strand_bias_odds_ratio: 0.0, mean_mapping_quality: 33.0, signal_to_noise: 2, af_high_quality_bases: 0.0002, af_adjusted: 0.0001, num_bases_3_prime_shift_for_deletions: 0, microsatellite: 0, microsatellite_length: 0, mean_mismatches_in_reads: 0.0, high_quality_variant_reads: 1, high_quality_total_reads: 6582, flank_seq_5_prime: "GAACATCACGGGCAATTAAA", flank_seq_3_prime: "GGGACCTAGATTTTAAGAGA", segment: "chr12:112450168-112450587", variant_type: "INV", duplication_rate: None, sv_info: Some(sv_info), distance_to_crispr_site: None },
        )
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
        assert_eq!(&format!("{}", ParsePairBiasError), "ParsePairBiasError");
    }

    #[test]
    fn test_parse_sv_info_err_display() {
        assert_eq!(&format!("{}", ParseSvInfoError), "ParseSvInfoError");
    }

    #[test]
    fn test_pair_bias_from_str_err() {
        assert_eq!(PairBias::from_str("3;0"), Err(ParsePairBiasError));
    }

    #[test]
    fn test_sv_info_from_str_err() {
        assert_eq!(SvInfo::from_str("0-0"), Err(ParseSvInfoError));
    }

    #[rstest]
    fn test_tumor_only_variant_ad_value(variants: Vec<TumorOnlyVariant>) {
        let expected = vec!["8103,1", "8103,1", "8210,1"];
        for (variant, ad) in variants.iter().zip(expected.iter()) {
            assert_eq!(&variant.ad_value(), ad);
        }
    }

    #[rstest]
    fn test_tumor_only_variant_ald_value(variants: Vec<TumorOnlyVariant>) {
        let expected = vec!["1,0", "0,1", "1,0"];
        for (variant, ad) in variants.iter().zip(expected.iter()) {
            assert_eq!(&variant.ald_value(), ad);
        }
    }

    #[rstest]
    fn test_tumor_only_variant_rd_value(variants: Vec<TumorOnlyVariant>) {
        let expected = vec!["2766,5280", "2766,5280", "3001,5130"];
        for (variant, ad) in variants.iter().zip(expected.iter()) {
            assert_eq!(&variant.rd_value(), ad);
        }
    }

    #[rstest]
    #[rustfmt::skip]
    fn test_tumor_only_variant_interval(variants: Vec<TumorOnlyVariant>) {
        let expected = vec![
            ("chr1", Range { start: 114713883 - 1, end: 114713883 } ),
            ("chr1", Range { start: 114713883 - 1, end: 114713883 } ),
            ("chr1", Range { start: 114713880 - 1, end: 114713880 } ),
        ];
        for (variant, (contig, range)) in variants.iter().zip(expected.iter()) {
            assert_eq!(&variant.contig(), contig);
            assert_eq!(&variant.range(), range);
        }
    }

    #[rstest]
    fn test_maybe_infinite_f32(
        variants: Vec<TumorOnlyVariant>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let input = PathBuf::from("tests/nras.var");

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(&input)
            .unwrap_or_else(|_| panic!("Could not open a reader for file path: {:?}", &input));

        let mut index = 0;
        let mut carry = csv::StringRecord::new();

        while reader
            .read_record(&mut carry)
            .expect("Failed to read the TumorOnlyVariant record.")
        {
            let variant: TumorOnlyVariant = carry
                .deserialize(None)
                .expect("Failed to deserialize the TumorOnlyVariant record.");
            assert_eq!(variants[index], variant);
            index += 1;
        }

        Ok(())
    }

    #[test]
    fn test_tumor_only_header() {
        let header = tumor_only_header("dna00001");
        let file = NamedTempFile::new().expect("Cannot create temporary file.");
        let _ = VcfWriter::from_path(&file.path(), &header, true, Format::VCF).unwrap();
        let reader = VcfReader::from_path(&file.path()).expect("Error opening tempfile.");
        let records = reader.header().header_records();
        assert_eq!(records.len(), 43);

        if let Some(HeaderRecord::Generic { ref key, ref value }) = records.get(2) {
            assert_eq!(key, &"source");
            assert_eq!(value, &[CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-"));
        } else {
            panic!("Expected source header record (tool, version) at this index.")
        }

        let samples = reader.header().samples();
        assert_eq!(samples.len(), 1);
        assert!(samples.iter().all(|&s| s == "dna00001".as_bytes()));
    }
}
