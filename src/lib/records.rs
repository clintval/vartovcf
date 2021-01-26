//! A module for serialization-deserialization friendly VarDict/VarDictJava data types.
use std::ops::Range;
use std::str::FromStr;

use bio_types::genome::{AbstractInterval, Position};
use linear_map::LinearMap;
use rust_htslib::bcf::{Header, HeaderRecord};
use rust_htslib::bcf::header::HeaderRecord::Generic;
use serde::{de::Error, Deserialize};

const CARGO_PKG_NAME: &str = env!("CARGO_PKG_NAME");
const CARGO_PKG_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Build a line of key-value pair data for the VCF header record.
fn hdr_rec_format_key_values(values: &LinearMap<String, String>) -> String {
    let mut pairs = String::new();
    for (key, value) in values.iter() { pairs.push_str(&format!(",{}={}", key, value)) };
    pairs
}

/// Render a header record into a string for writing to a VCF.
fn hdr_rec_as_str(rec: &HeaderRecord) -> String {
    match rec { // TODO: Is this really how you build the HeaderRecords?
        HeaderRecord::Contig     {key, values} => format!("##FILTER=<ID={},{}>", key, hdr_rec_format_key_values(values)),
        HeaderRecord::Filter     {key, values} => format!("##FILTER=<ID={},{}>", key, hdr_rec_format_key_values(values)),
        HeaderRecord::Format     {key, values} => format!("##FILTER=<ID={},{}>", key, hdr_rec_format_key_values(values)),
        HeaderRecord::Generic    {key, value}  => format!("##{}={}", key, value),
        HeaderRecord::Info       {key, values} => format!("##FILTER=<ID={},{}>", key, hdr_rec_format_key_values(values)),
        HeaderRecord::Structured {key, values} => format!("##FILTER=<ID={},{}>", key, hdr_rec_format_key_values(values)),
    }
}

/// Create a VCF header for VarDict/VarDictJava in tumor-only mode.
pub fn tumor_only_header(sample: String) -> Header {
    let source = [CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-");
    let mut header = Header::default();
    header.push_sample(&sample.into_bytes());
    let records: Vec<HeaderRecord> = vec![
        Generic { key: String::from("fileFormat"), value: String::from("VCFv4.2") },
        Generic { key: String::from("source"), value: source }
    ];
    for record in &records { header.push_record(hdr_rec_as_str(&record).as_bytes()); }
    header
}

/// A record of output from VarDict/VarDictJava run in tumor-only mode.
#[derive(Debug, Deserialize)]
pub struct TumorOnlyVariant<'a> {
    sample: &'a str,
    interval_name: &'a str,
    contig: &'a str,
    start: u64,
    end: u64,
    ref_allele: &'a str,
    alt_allele: &'a str,
    depth: u32,
    alt_depth: u32,
    ref_forward: u32,
    ref_reverse: u32,
    alt_forward: u32,
    alt_reverse: u32,
    gt: &'a str,
    af: f32,
    // Semi-colon separated integers, could be a struct
    strand_bias: &'a str,
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
    flank_seq_5_prime: &'a str,
    flank_seq_3_prime: &'a str,
    // Position format, unused in VCF.
    segment: &'a str,
    variant_type: &'a str,
    duplication_rate: &'a str,
    // Either a zero, or triple of ints separated with "-".
    // We can deserialize this further into a struct of SV details.
    sv_details: &'a str,
    #[serde(default)]
    distance_to_crispr_site: Option<u32>,
}

impl<'a> AbstractInterval for TumorOnlyVariant<'a> {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn range(&self) -> Range<Position> {
        Range { start: self.start, end: self.end }
    }
}

/// A record of output from VarDict/VarDictJava run in amplicon-aware mode.
#[derive(Debug, Deserialize)]
struct AmpliconVariant<'a> {
    sample: &'a str,
    interval_name: &'a str,
    contig: &'a str,
    start: u64,
    end: u64,
    ref_allele: &'a str,
    alt_allele: &'a str,
    depth: u32,
    alt_depth: u32,
    ref_forward: u32,
    ref_reverse: u32,
    alt_forward: u32,
    alt_reverse: u32,
    gt: &'a str,
    af: f32,
    strand_bias: &'a str,
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
    flank_seq_5_prime: &'a str,
    flank_seq_3_prime: &'a str,
    // Position format, unused in VCF.
    segment: &'a str,
    variant_type: &'a str,
    num_amplicons_supporting_variant: u32,
    total_amplicons_overlapping: u32,
    num_amplicons_rare: u32,
    // A zero or one, need a custom deserializer.
    top_variant_in_amplicon_does_not_match: bool,
}

impl<'a> AbstractInterval for AmpliconVariant<'a> {
    fn contig(&self) -> &str {
        &self.contig
    }
    fn range(&self) -> Range<Position> {
        Range { start: self.start, end: self.end }
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
