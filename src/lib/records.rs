//! A module for serialization-deserialization friendly VarDict/VarDictJava data types.
use std::ops::Range;
use std::str::FromStr;

use bio_types::genome::{AbstractInterval, Position};
use serde::{de::Error, Deserialize};

/// A record of output from VarDict/VarDictJava run in tumor-only mode.
#[derive(Debug, Deserialize)]
pub struct TumorOnlyVariant {
    sample: String,
    interval_name: String,
    contig: String,
    start: u64,
    end: u64,
    ref_allele: String,
    alt_allele: String,
    depth: u32,
    alt_depth: u32,
    ref_forward: u32,
    ref_reverse: u32,
    alt_forward: u32,
    alt_reverse: u32,
    gt: String,
    af: f32,
    // Semi-colon separated integers, could be a struct
    strand_bias: String,
    mean_position_in_read: f32,
    stdev_position_in_read: f32,
    mean_base_quality: f32,
    stdev_base_quality: f32,
    strand_bias_p_value: f32,
    #[serde(deserialize_with = "maybe_infinite_f32")]
    strand_bias_odds_ratio: f32,
    mean_mapping_quality: f32,
    signal_to_noise: u32,
    af_high_quality_bases: f32,
    af_adjusted: f32,
    num_bases_3_prime_shift_for_deletions: u32,
    microsatellite: u32,
    microsatellite_length: u32,
    mean_mismatches_in_reads: f32,
    high_quality_variant_reads: u32,
    high_quality_total_reads: u32,
    flank_seq_5_prime: String,
    flank_seq_3_prime: String,
    // Position format, unused in VCF.
    segment: String,
    variant_type: String,
    duplication_rate: String,
    // Either a zero, or triple of ints separated with "-".
    // We can deserialize this further into a struct of SV details.
    sv_details: String,
    #[serde(default)]
    distance_to_crispr_site: Option<String>,
}

impl AbstractInterval for TumorOnlyVariant {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn range(&self) -> Range<Position> {
        Range { start: self.start.clone(), end: self.end.clone() }
    }
}

/// A record of output from VarDict/VarDictJava run in amplicon-aware mode.
#[derive(Debug, Deserialize)]
struct AmpliconVariant {
    sample: String,
    interval_name: String,
    contig: String,
    start: u64,
    end: u64,
    ref_allele: String,
    alt_allele: String,
    depth: u32,
    alt_depth: u32,
    ref_forward: u32,
    ref_reverse: u32,
    alt_forward: u32,
    alt_reverse: u32,
    gt: String,
    af: f32,
    strand_bias: String,
    mean_position_in_read: f32,
    stdev_position_in_read: f32,
    mean_base_quality: f32,
    stdev_base_quality: f32,
    strand_bias_p_value: f32,
    #[serde(deserialize_with = "maybe_infinite_f32")]
    strand_bias_odds_ratio: f32,
    mean_mapping_quality: f32,
    signal_to_noise: u32,
    af_high_quality_bases: f32,
    af_adjusted: f32,
    num_bases_3_prime_shift_for_deletions: u32,
    microsatellite: u32,
    microsatellite_length: u32,
    mean_mismatches_in_reads: f32,
    high_quality_variant_reads: u32,
    high_quality_total_reads: u32,
    flank_seq_5_prime: String,
    flank_seq_3_prime: String,
    // Position format, unused in VCF.
    segment: String,
    variant_type: String,
    num_amplicons_supporting_variant: u32,
    total_amplicons_overlapping: u32,
    num_amplicons_rare: u32,
    // A zero or one, need a custom deserializer.
    top_variant_in_amplicon_does_not_match: bool,
}

impl AbstractInterval for AmpliconVariant {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn range(&self) -> Range<Position> {
        Range { start: self.start.clone(), end: self.end.clone() }
    }
}

/// Deserialize a possibly infinite float into a <f32> or return a custom error.
///
/// The following cases are handled:
///
/// * `"Inf"`: floating point infinity
/// * `"-Inf"`: floating point negative infinity
/// * `<other>`: a non-infinite floating point number
fn maybe_infinite_f32<'de, D>(deserializer: D) -> Result<f32, D::Error>
    where D: serde::Deserializer<'de> {
    let s: &str = Deserialize::deserialize(deserializer)?;
    if s == "Inf" { Ok(f32::INFINITY) }
    else if s == "-Inf" { Ok(f32::NEG_INFINITY) }
    else { f32::from_str(&s).map_err(D::Error::custom) }
}
