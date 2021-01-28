//! A module for serialization-deserialization friendly VarDict/VarDictJava data types.
use std::ops::Range;
use std::str::FromStr;

use anyhow::Result;
use bio_types::genome::{AbstractInterval, Position};
use rust_htslib::bcf::Header;
use serde::{de::Error, Deserialize};

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
        .expect("Could not deserialize a maybe-infinite (Inf, -Inf) floating point number.");
    if s == "Inf" {
        Ok(f32::INFINITY)
    } else if s == "-Inf" {
        Ok(f32::NEG_INFINITY)
    } else {
        f32::from_str(&s).map_err(D::Error::custom)
    }
}

/// A record of output from VarDict/VarDictJava run in tumor-only mode.
#[derive(Debug, Deserialize)]
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
    pub depth: u32,
    /// The total alternate depth at this call locus.
    pub alt_depth: u32,
    /// The number of forward reads supporting the reference call.
    pub ref_forward: u32,
    /// The number of reverse reads supporting the reference call.
    pub ref_reverse: u32,
    /// The number of forward reads supporting the alternate call.
    pub alt_forward: u32,
    /// The number of alternate reads supporting the alternate call.
    pub alt_reverse: u32,
    /// The call's genotype.
    pub gt: &'a str,
    /// The allele frequency of the alternate allele.
    pub af: f32,
    /// Strand bias status. That will take the values [0-2];[0-2] (_e.g._ "0;2", "2;1"). The first
    /// value refers to reads that support the reference allele, and the second to reads that
    /// support the variant allele.
    ///
    /// * `0`: small total count of reads (less than 12 for the sum of forward and reverse reads)
    /// * `1`: strand bias
    /// * `2`: no strand bias
    pub strand_bias: &'a str,
    /// The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads
    /// that support the variant call.
    pub mean_position_in_read: f32,
    /// The standard deviation of the distance to the nearest 5 or 3 prime read end (whichever is
    /// closer) in all reads that support the variant call.
    pub stdev_position_in_read: f32,
    /// The mean base quality (phred) of all bases that directly support the variant call.
    pub mean_base_quality: f32,
    /// The standard deviation of the base quality (phred)) of all bases that directly support
    /// the variant call.
    pub stdev_base_quality: f32,
    /// The Fisher test p-value for if you should reject the hypothesis that there is no strand
    /// bias. Non multiple hypothesis test corrected.
    pub strand_bias_p_value: f32,
    #[serde(deserialize_with = "maybe_infinite_f32")]
    /// The odds ratio for strand bias.
    pub strand_bias_odds_ratio: f32,
    /// The mean mapping quality (phred) of all reads that directly support the variant call.
    pub mean_mapping_quality: f32,
    /// The signal to noise ratio.
    pub signal_to_noise: u32,
    /// Allele frequency calculated using only high quality bases. Lossy due to rounding.
    pub af_high_quality_bases: f32,
    /// Adjusted allele frequency for indels due to local realignment. Lossy due to rounding.
    pub af_adjusted: f32,
    /// The number of bases to be shifted 3-prime for deletions due to alternative alignment(s).
    pub num_bases_3_prime_shift_for_deletions: u32,
    /// Whether the variant call is in a microsatellite (MSI) or not. Greater than 1 indicates MSI.
    pub microsatellite: u32,
    /// The length of the microSatellite in base pairs of reference genome.
    pub microsatellite_length: u32,
    /// The length of the microSatellite in base pairs of reference genome.
    pub mean_mismatches_in_reads: f32,
    /// The number of high quality reads supporting the variant call.
    pub high_quality_variant_reads: u32,
    /// The number of high quality reads at the locus of the variant call.
    pub high_quality_total_reads: u32,
    /// 5-prime reference flanking sequence.
    pub flank_seq_5_prime: &'a str,
    /// 3-prime reference flanking sequence.
    pub flank_seq_3_prime: &'a str,
    /// The position formatted interval of the variant calling target.
    pub segment: &'a str,
    /// The type of variant this call is.
    pub variant_type: &'a str,
    /// The duplication rate, if this call is a duplication.
    pub duplication_rate: &'a str,
    /// The details of the structural variant as a 0 or as a triplet of ints separated with "=".
    // We can deserialize this further into a struct of SV details.
    pub sv_details: &'a str,
    #[serde(default)]
    /// The distance in reference genome base pairs to the nearest CRISPR-site (CRISPR-mode only).
    pub distance_to_crispr_site: Option<u32>,
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
pub fn tumor_only_header(sample: String) -> Header {
    let source = [CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-");
    let mut header = Header::default();
    header.push_sample(&sample.into_bytes());
    header.remove_filter(b"PASS");
    header.push_record(format!("##source={}", source).as_bytes());
    header.push_record(r#"##INFO=<ID=SAMPLE,Number=1,Type=String,Description="Sample name (with whitespace translated to underscores)">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=BIAS,Number=1,Type=String,Description="Strand bias flags, see VarDictJava documentation for more information">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=REFBIAS,Number=1,Type=String,Description="Strand bias in reads that support the reference call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=VARBIAS,Number=1,Type=String,Description="Strand bias in reads that support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=PMEAN,Number=1,Type=Float,Description="The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=PSTD,Number=1,Type=Float,Description="The standard deviation of the distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=QUAL,Number=1,Type=Float,Description="The mean base quality (phred) of all bases that directly support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=QSTD,Number=1,Type=Float,Description="The standard deviation of the base quality (phred)) of all bases that directly support the variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SBF,Number=1,Type=Float,Description="The Fisher test p-value for if you should reject the hypothesis that there is no strand bias. Non multiple hypothesis test corrected">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description="The odds ratio for strand bias for this variant call">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=MQ,Number=1,Type=Float,Description="The mean mapping quality (phred) of all reads that directly support the variant call">"#.as_bytes());
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
    header.push_record(r#"##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The structural variant type (INV DUP DEL INS FUS), if this call is a structural variant">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="The length of stuctural variant in base pairs of reference genome, if this call is a structural variant">"#.as_bytes());
    header.push_record(r#"##INFO=<ID=DUPRATE,Number=1,Type=Float,Description="The duplication rate, if this call is a duplication">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=PASS,Description="The variant call has passed all filters and may be considered for downstream analysis">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=q22.5,Description="The mean base quality (phred) of all bases that directly support this variant call is below 22.5">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=Q10,Description="The mean mapping quality (phred) in reads that suppr 10">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=MSI12,Description="The variant call is in a microsatellite region with 12 non-monomer MSI or 13 monomer MSI">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=NM8.0,Description="The mean mismatches in reads that support the variant call is >= 8.0, and might be a false positive or contamination">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=InGap,Description="The variant call is in a deletion gap, and might be a false positive">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=InIns,Description="The variant call was found to be adjacent to an insertion variant">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=Cluster0bp,Description="At least two variant calls are within 0 base pairs from each other in the reference sequence coordinate system">"#.as_bytes());
    header.push_record(r#"##FILTER=<ID=LongMSI,Description="The variant call is flanked by a long A/T stretch (>=14 base pairs)">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=GT,Number=1,Type=String,Description="The genotype for this sample for this variant call ">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=DP,Number=1,Type=Integer,Description="The total allele depth at this location which potentially includes No-calls">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=VD,Number=1,Type=Integer,Description="The variant allele depth at this location">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=AD,Number=R,Type=Integer,Description="The allelic depths for the REF and ALT alleles">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=RD,Number=2,Type=Integer,Description="The number of reference forward and reverse reads">"#.as_bytes());
    header.push_record(r#"##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="The number of variant call forward and reverse reads">"#.as_bytes());
    header
}

#[cfg(test)]
mod tests {
    use rust_htslib::bcf::Writer;
    use rust_htslib::bcf::{Format, HeaderRecord, Read, Reader};
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_tumor_only_header() {
        let header = tumor_only_header("DNA00001".into());
        let file = NamedTempFile::new().expect("Cannot create temporary file.");
        let _ = Writer::from_path(&file.path(), &header, true, Format::VCF).unwrap();
        let reader = Reader::from_path(&file.path()).expect("Error opening tempfile.");
        let records = reader.header().header_records();
        assert_eq!(records.len(), 44);
        match records[2] {
            HeaderRecord::Generic { ref key, ref value } => {
                assert_eq!(key, &"source");
                assert_eq!(value, &[CARGO_PKG_NAME, CARGO_PKG_VERSION].join("-"));
            }
            _ => {
                panic!(
                    "Expected source header record (tool, version), but found: {:?}",
                    records[2]
                );
            }
        }
        let samples = reader.header().samples();
        assert_eq!(samples.len(), 1);
        assert!(samples.iter().all(|&s| s == "DNA00001".as_bytes()));
    }
}
