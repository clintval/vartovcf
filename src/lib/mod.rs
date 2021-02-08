//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]

use std::error;
use std::fmt::Debug;
use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::Result;
use csv::ReaderBuilder;
use log::*;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::Format;
use rust_htslib::bcf::Writer as VcfWriter;
use strum::{EnumString, EnumVariantNames, ToString as EnumToString};

use crate::fai::{reference_contigs_to_vcf_header, reference_path_to_vcf_header};
use crate::io::has_gzip_ext;
use crate::progress_logger::{ProgressLogger, RecordLogger, DEFAULT_LOG_EVERY};
use crate::record::tumor_only_header;
use crate::record::TumorOnlyVariant;

pub mod fai;
pub mod io;
pub mod progress_logger;
pub mod record;

/// The default log level for the `vartovcf` tool.
pub const DEFAULT_LOG_LEVEL: &str = "info";

/// The valid structural variation (SV) type values for the `SVTYPE` FORMAT field.
pub const VALID_SV_TYPES: &[&str] = &["BND", "CNV", "DEL", "DUP", "INS", "INV"];

/// Namespace for path parts and extensions.
pub mod path {

    /// The Gzip extension.
    pub const GZIP_EXTENSION: &str = "gz";
}

/// The variant calling modes for VarDict/VarDictJava.
#[derive(Clone, Copy, Debug, EnumString, EnumToString, EnumVariantNames, PartialEq, PartialOrd)]
pub enum VarDictMode {
    /// The amplicon variant calling mode.
    Amplicon,
    /// The tumor-normal variant calling mode.
    TumorNormal,
    /// The tumor-only variant calling mode.
    TumorOnly,
}

/// Runs the tool `vartovcf` on an input VAR file and writes the records to an output VCF file.
///
/// # Arguments
///
/// * `input` - The input VAR file or stream
/// * `output` - The output VCF file or stream
/// * `fasta` - The reference sequence FASTA file, must be indexed
/// * `sample` - The sample name
/// * `mode` - The variant calling modes for VarDict/VarDictJava
///
/// # Returns
///
/// Returns the result of the execution with an integer exit code for success (0).
///
pub fn vartovcf<I, R>(
    input: I,
    output: Option<PathBuf>,
    fasta: R,
    sample: &str,
    mode: &VarDictMode,
) -> Result<i32, Box<dyn error::Error>>
where
    I: Read,
    R: AsRef<Path> + Debug,
{
    assert_eq!(
        mode,
        &VarDictMode::TumorOnly,
        "The only mode currently supported is [TumorOnly]."
    );
    let mut header = tumor_only_header(&sample);
    reference_contigs_to_vcf_header(&fasta, &mut header);
    reference_path_to_vcf_header(&fasta, &mut header)
        .expect("Could not add the FASTA file path to the VCF header");

    let mut progress = ProgressLogger::new("processed", "variant records", DEFAULT_LOG_EVERY);

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(input);

    let mut writer = match output {
        Some(output) => VcfWriter::from_path(&output, &header, !has_gzip_ext(&output), Format::VCF),
        None => VcfWriter::from_stdout(&header, true, Format::VCF),
    }
    .expect("Could not build a VCF writer.");

    let mut carry = csv::StringRecord::new();
    let mut variant = writer.empty_record();

    while reader.read_record(&mut carry)? {
        let var: TumorOnlyVariant = carry.deserialize(None)?;
        if var.sample != sample {
            let message = format!("Expected sample '{}' found '{}'", sample, var.sample);
            progress.emit()?;
            return Err(message.into());
        };

        let rid = writer.header().name2rid(var.contig.as_bytes()).unwrap();
        let genotypes = &[
            GenotypeAllele::UnphasedMissing,
            GenotypeAllele::UnphasedMissing,
        ];

        variant.set_rid(Some(rid));
        variant.set_pos(var.start as i64 - 1);
        variant.set_alleles(&[var.ref_allele.as_bytes(), var.alt_allele.as_bytes()])?;
        variant.push_genotypes(genotypes).unwrap();

        variant.push_info_float(b"ADJAF", &[var.af_adjusted])?;
        variant.push_info_string(b"BIAS", &[var.strand_bias.to_string().as_bytes()])?;
        if let Some(duplication_rate) = var.duplication_rate {
            variant.push_info_float(b"DUPRATE", &[duplication_rate])?;
        }
        variant.push_info_float(b"HIAF", &[var.af_high_quality_bases])?;
        variant.push_info_integer(b"HICNT", &[var.high_quality_variant_reads])?;
        variant.push_info_integer(b"HICOV", &[var.high_quality_total_reads])?;
        variant.push_info_string(b"LSEQ", &[var.flank_seq_5_prime.as_bytes()])?;
        variant.push_info_float(b"MQ", &[var.mean_mapping_quality])?;
        variant.push_info_integer(b"MSI", &[var.microsatellite])?;
        variant.push_info_integer(b"MSILEN", &[var.microsatellite_length])?;
        variant.push_info_float(b"NM", &[var.mean_mismatches_in_reads])?;
        variant.push_info_float(b"ODDRATIO", &[var.strand_bias_odds_ratio])?;
        variant.push_info_float(b"PMEAN", &[var.mean_position_in_read])?;
        variant.push_info_float(b"PSTD", &[var.stdev_position_in_read])?;
        variant.push_info_float(b"QSTD", &[var.stdev_base_quality])?;
        variant.push_info_float(b"QUAL", &[var.mean_base_quality])?;
        // variant.push_info_integer(b"REFBIAS", &[var.])?;
        variant.push_info_string(b"RSEQ", &[var.flank_seq_3_prime.as_bytes()])?;
        variant.push_info_float(b"SBF", &[var.strand_bias_p_value])?;
        variant.push_info_integer(b"SHIFT3", &[var.num_bases_3_prime_shift_for_deletions])?;
        variant.push_info_integer(b"SN", &[var.signal_to_noise])?;

        // NB: If you do not explicitly clear the buffers, you'll end up with stale references.
        if let Some(sv_info) = &var.sv_info {
            variant.push_info_integer(b"SPLITREAD", &[sv_info.supporting_split_reads])?;
            variant.push_info_integer(b"SPANPAIR", &[sv_info.supporting_pairs])?;
        } else {
            variant.clear_info_integer(b"SPLITREAD")?;
            variant.clear_info_integer(b"SPANPAIR")?;
        }
        if VALID_SV_TYPES.contains(&var.variant_type) {
            variant.push_info_integer(b"SVLEN", &[var.length()])?;
            variant.push_info_string(b"SVTYPE", &[var.variant_type.as_bytes()])?;
        } else {
            variant.clear_info_integer(b"SVLEN")?;
            variant.clear_info_integer(b"SVTYPE")?;
        }
        // variant.push_info_integer(b"VARBIAS", &[var.])?;

        variant.push_format_integer(b"VD", &[var.alt_depth])?;
        variant.push_format_integer(b"DP", &[var.depth])?;
        variant.push_format_string(b"AD", &[var.ad_value().as_bytes()])?;
        variant.push_format_string(b"ALD", &[var.ald_value().as_bytes()])?;
        variant.push_format_string(b"RD", &[var.rd_value().as_bytes()])?;

        writer.write(&variant)?;
        progress.observe()?;
    }

    progress.emit()?;

    Ok(0)
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use anyhow::Result;
    use file_diff::diff;
    use tempfile::NamedTempFile;

    use super::VarDictMode::TumorOnly;
    use super::*;

    #[test]
    fn test_vartovcf_run() -> Result<(), Box<dyn std::error::Error>> {
        let sample = "dna00001";
        let input = BufReader::new(File::open("tests/nras.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file.");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            &sample,
            &TumorOnly,
        )?;
        assert_eq!(exit, 0);
        assert!(diff(&output.path().to_str().unwrap(), "tests/nras.vcf"));
        Ok(())
    }

    #[test]
    fn test_when_incorrect_sample() {
        let sample = "XXXXXXXX";
        let input = BufReader::new(File::open("tests/nras.var").unwrap());
        let output = NamedTempFile::new().expect("Cannot create temporary file.");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            &sample,
            &TumorOnly,
        );
        assert!(result.is_err());
    }
}
