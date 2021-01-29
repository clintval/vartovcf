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

use crate::fai::{contigs_to_vcf_header, fai_file, vcf_contig_header_records};
use crate::record::tumor_only_header;
use crate::record::TumorOnlyVariant;

pub mod fai;
pub mod io;
pub mod record;

/// The default log level for the `vartovcf` tool.
pub const DEFAULT_LOG_LEVEL: &str = "info";

/// Namespace for path parts and extensions.
mod path {

    /// The Gzip extension.
    pub const GZIP_EXTENSION: &str = "gz";
}

/// Runs the tool `vartovcf` on an input VAR file and writes the records to an output VCF file.
///
/// # Arguments
///
/// * `input` - The input VAR file or stream
/// * `output` - The output VCF file or stream
/// * `fasta` - The reference sequence FASTA file, must be indexed
/// * `sample` - The sample name
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
) -> Result<i32, Box<dyn error::Error>>
where
    I: Read,
    R: AsRef<Path> + Debug,
{
    let contigs = vcf_contig_header_records(fai_file(fasta))?;
    let header = contigs_to_vcf_header(&contigs, tumor_only_header(sample));
    let plain_text = if let Some(path) = &output {
        !io::has_gzip_ext(&path)
    } else {
        true
    };

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(input);

    let mut writer = match output {
        Some(output) => VcfWriter::from_path(&output, &header, plain_text, Format::VCF),
        None => VcfWriter::from_stdout(&header, plain_text, Format::VCF),
    }
    .expect("Could not build a VCF writer.");

    let mut count = 0;
    let mut carry = csv::StringRecord::new();
    let mut variant = writer.empty_record();

    while reader.read_record(&mut carry)? {
        let var: TumorOnlyVariant = carry.deserialize(None)?;
        let rid = writer.header().name2rid(var.contig.as_bytes()).unwrap();

        // TODO: Set these correctly based on parsing the phasing information.
        let alleles = &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)];

        variant.set_rid(Some(rid));
        variant.set_pos(var.start as i64 - 1);

        variant.push_info_float(b"PMEAN", &[var.mean_position_in_read])?;
        variant.push_info_float(b"PSTD", &[var.stdev_position_in_read])?;
        // variant.push_info_integer(b"BIAS", &[var.])?;
        // variant.push_info_integer(b"REFBIAS", &[var.])?;
        // variant.push_info_integer(b"VARBIAS", &[var.])?;
        variant.push_info_float(b"QUAL", &[var.mean_base_quality])?;
        variant.push_info_float(b"QSTD", &[var.stdev_base_quality])?;
        variant.push_info_float(b"SBF", &[var.strand_bias_p_value])?;
        variant.push_info_float(b"ODDRATIO", &[var.strand_bias_odds_ratio])?;
        variant.push_info_float(b"MQ", &[var.mean_mapping_quality])?;
        variant.push_info_integer(b"SN", &[var.signal_to_noise])?;
        variant.push_info_float(b"HIAF", &[var.af_high_quality_bases])?;
        variant.push_info_float(b"ADJAF", &[var.af_adjusted])?;
        variant.push_info_integer(b"SHIFT3", &[var.num_bases_3_prime_shift_for_deletions])?;
        variant.push_info_integer(b"MSI", &[var.microsatellite])?;
        variant.push_info_integer(b"MSILEN", &[var.microsatellite_length])?;
        variant.push_info_float(b"NM", &[var.mean_mismatches_in_reads])?;
        variant.push_info_string(b"LSEQ", &[var.flank_seq_5_prime.as_bytes()])?;
        variant.push_info_string(b"RSEQ", &[var.flank_seq_3_prime.as_bytes()])?;
        variant.push_info_integer(b"HICNT", &[var.high_quality_variant_reads])?;
        variant.push_info_integer(b"HICOV", &[var.high_quality_total_reads])?;
        // variant.push_info_integer(b"SPLITREAD", &[var.])?;
        // variant.push_info_integer(b"SPANPAIR", &[var.])?;
        // variant.push_info_integer(b"SVTYPE", &[var.])?;
        // variant.push_info_integer(b"SVLEN", &[var.])?; // Ensure is negative for deletion
        variant.push_info_float(b"DUPRATE", &[var.duplication_rate])?;

        variant.push_genotypes(alleles).unwrap();
        variant.set_alleles(&[var.ref_allele.as_bytes(), var.alt_allele.as_bytes()])?;

        writer.write(&variant)?;
        count += 1;
    }

    info!("Processed {} variant records", count);
    Ok(0)
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use anyhow::Result;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_vartovcf_run() -> Result<(), Box<dyn std::error::Error>> {
        let sample = "DNA00001";
        let input = BufReader::new(File::open("tests/nras.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file.");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(input, Some(output.path().into()), &reference, &sample)?;
        assert_eq!(exit, 0);
        // TODO: Check output against a plain text file, line-by-line?
        Ok(())
    }
}
