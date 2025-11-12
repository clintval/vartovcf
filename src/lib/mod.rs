//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]

use anyhow::Result;
use csv::ReaderBuilder;
use log::*;
use proglog::ProgLogBuilder;
use rust_htslib::bcf::Format;
use rust_htslib::bcf::Writer as VcfWriter;
use std::collections::HashSet;
use std::error;
use std::fmt::Debug;
use std::io::Read;
use std::path::{Path, PathBuf};
use strum::Display;
use strum::{EnumString, VariantNames};

use crate::fai::{fasta_contigs_to_vcf_header, fasta_path_to_vcf_header};
use crate::io::has_gzip_ext;
use crate::record::TumorOnlyVariant;
use crate::record::tumor_only_header;

pub mod fai;
pub mod io;
pub mod record;

/// The valid structural variation (SV) type values for the `SVTYPE` FORMAT field.
pub const VALID_SV_TYPES: &[&str] = &["BND", "CNV", "DEL", "DUP", "INS", "INV"];

/// Namespace for path parts and extensions.
pub mod path {

    /// The Gzip extension.
    pub const GZIP_EXTENSION: &str = "gz";
}

/// The variant calling modes for VarDict/VarDictJava.
#[derive(Clone, Copy, Debug, Display, EnumString, VariantNames, PartialEq, PartialOrd)]
pub enum VarDictMode {
    /// The amplicon variant calling mode.
    //Amplicon,
    /// The tumor-normal variant calling mode.
    //TumorNormal,
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
/// * `skip_non_variants` - Skip non-variant sites (where ref_allele == alt_allele)
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
    skip_non_variants: bool,
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

    let mut header = tumor_only_header(sample);

    fasta_contigs_to_vcf_header(&fasta, &mut header);
    fasta_path_to_vcf_header(&fasta, &mut header).expect("Adding FASTA path to header failed!");

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(input);

    let mut writer = match output {
        Some(output) => VcfWriter::from_path(&output, &header, !has_gzip_ext(&output), Format::Vcf),
        None => VcfWriter::from_stdout(&header, true, Format::Vcf),
    }
    .expect("Could not build a VCF writer!");

    let progress = ProgLogBuilder::new()
        .name("main")
        .verb("Processed")
        .noun("variant records")
        .unit(100_000)
        .build();

    let mut carry = csv::StringRecord::new();
    let mut variant = writer.empty_record();
    let mut seen: HashSet<String> = HashSet::new();

    while reader.read_record(&mut carry)? {
        if carry.iter().collect::<Vec<&str>>()[5].is_empty() {
            continue; // If the 5th field is empty, it's a record we need to avoid deserializing.
        }

        let var: TumorOnlyVariant = carry
            .deserialize(None)
            .expect("Could not deserialize record!");

        // Skip non-variant sites if the flag is set
        if skip_non_variants && var.ref_allele == var.alt_allele {
            continue;
        }

        let key = format!(
            "{}-{}-{}-{}",
            var.contig, var.start, var.ref_allele, var.alt_allele
        );
        if !seen.insert(key) {
            continue; // Skip this record if we have seen this variant before.
        }

        if var.sample != sample {
            let message = format!("Expected sample '{}' found '{}'!", sample, var.sample);
            return Err(message.into());
        };

        let rid = writer.header().name2rid(var.contig.as_bytes()).unwrap();

        variant.set_rid(Some(rid));
        variant.set_pos(var.start as i64 - 1);
        variant.set_alleles(&[
            var.ref_allele.as_bytes(),
            var.alt_allele_for_vcf().as_bytes(),
        ])?;

        if var.alt_depth == 0 {
            variant.set_qual(0.0)
        } else {
            let qual = (var.alt_depth as f32).ln() / 2.0_f32.ln() * var.base_quality_mean;
            variant.set_qual(qual)
        }

        variant.push_info_integer(b"END", &[var.end as i32])?;

        if VALID_SV_TYPES.contains(&var.variant_type) {
            variant.push_info_integer(b"SVLEN", &[var.length()])?;
            variant.push_info_string(b"SVTYPE", &[var.variant_type.as_bytes()])?;
        } else {
            // NB: Without clearing the fields, you'll end up with stale references.
            variant.clear_info_integer(b"SVLEN")?;
            variant.clear_info_integer(b"SVTYPE")?;
        }

        variant.push_genotypes(var.gt_value(0.25))?;
        variant.push_format_integer(b"AD", &var.ad_value())?;
        variant.push_format_integer(b"DP", &[var.depth])?;
        variant.push_format_integer(b"VD", &[var.alt_depth])?;

        writer.write(&variant)?;
        progress.record();
    }

    Ok(0)
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use anyhow::Result;
    use file_diff::diff;
    use pretty_assertions::assert_eq;
    use tempfile::NamedTempFile;

    use super::VarDictMode::TumorOnly;
    use super::*;

    #[test]
    fn test_vartovcf_run() -> Result<(), Box<dyn std::error::Error>> {
        let sample = "dna00001";
        let input = BufReader::new(File::open("tests/calls.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            &sample,
            &TumorOnly,
            false,
        )?;
        assert_eq!(exit, 0);
        assert!(diff(&output.path().to_str().unwrap(), "tests/calls.vcf"));
        Ok(())
    }

    #[test]
    fn test_when_incorrect_sample() {
        let sample = "XXXXXXXX";
        let input = BufReader::new(File::open("tests/calls.var").unwrap());
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            &sample,
            &TumorOnly,
            false,
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_skip_non_variants_flag() -> Result<(), Box<dyn std::error::Error>> {
        let sample = "dna00001";
        let input = BufReader::new(File::open("tests/calls.g.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            &sample,
            &TumorOnly,
            true,
        )?;
        assert_eq!(exit, 0);

        // Read the output and verify no non-variant sites exist
        use rust_htslib::bcf::{Read, Reader as VcfReader};
        let mut reader = VcfReader::from_path(&output.path()).expect("Error opening output file!");

        let mut record_count = 0;
        for record_result in reader.records() {
            let record = record_result?;
            let alleles = record.alleles();

            // Verify that REF != ALT (i.e., no non-variant sites)
            // ALT should not be "." which represents non-variants
            assert_ne!(alleles.len(), 0, "Record should have alleles");
            if alleles.len() >= 2 {
                let alt = alleles[1];
                assert_ne!(
                    alt, b".",
                    "Non-variant site found when --skip-non-variants is true"
                );
            }
            record_count += 1;
        }

        // Verify we got some records (not all were filtered)
        assert!(record_count > 0, "No records found in output");
        Ok(())
    }
}
