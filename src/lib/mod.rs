//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]

use ahash::AHashSet;
use anyhow::Result;
use csv::{Reader, ReaderBuilder, StringRecord};
use proglog::ProgLogBuilder;
use rust_htslib::bcf::Format;
use rust_htslib::bcf::Writer as VcfWriter;
use rust_htslib::bcf::record::Numeric;
use std::error;
use std::fmt::Debug;
use std::io::Read;
use std::path::{Path, PathBuf};
use strum::Display;
use strum::{EnumString, VariantNames};

use crate::fai::{fasta_contigs_to_vcf_header, fasta_path_to_vcf_header};
use crate::io::has_gzip_ext;
use crate::record::MIN_HOM_ALT_AF;
use crate::record::TumorOnlyVariant;
use crate::record::tumor_only_header;

pub mod fai;
pub mod io;
pub mod record;

/// The structural variant types VarDict writes, which are the values of the `SVTYPE` INFO field.
pub const VALID_SV_TYPES: &[&str] = &["DEL", "DUP", "INV"];

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

    // Read the first record before the writer emits the header so a refused stream writes nothing.
    let mut carry = StringRecord::new();
    let mut peeked = read_fisher_record(&mut reader, &mut carry)?;

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

    let mut seen: AHashSet<String> = AHashSet::new();

    while std::mem::take(&mut peeked) || read_fisher_record(&mut reader, &mut carry)? {
        if carry.get(5).is_none_or(|f| f.is_empty()) {
            continue; // If the 5th field is empty, it's a record we need to avoid deserializing.
        }

        let var: TumorOnlyVariant = carry.deserialize(None).map_err(describe_parse_error)?;

        if skip_non_variants && var.ref_allele == var.alt_allele {
            continue;
        }

        let key = format!(
            "{}-{}-{}-{}-{}",
            var.contig, var.start, var.end, var.ref_allele, var.alt_allele
        );
        if !seen.insert(key) {
            continue; // Skip this record if we have seen this variant before.
        }

        if var.sample != sample {
            let message = format!("Expected sample '{}' found '{}'!", sample, var.sample);
            return Err(message.into());
        };

        let rid = writer.header().name2rid(var.contig.as_bytes()).unwrap();

        let mut variant = writer.empty_record();
        variant.set_rid(Some(rid));
        variant.set_pos(var.start as i64 - 1);
        variant.set_alleles(&[
            var.ref_allele.as_bytes(),
            var.alt_allele_for_vcf().as_bytes(),
        ])?;

        variant.set_qual(f32::missing());

        if let Some(class) = var.variant_class() {
            variant.push_info_string(b"TYPE", &[class.as_bytes()])?;
        }

        if VALID_SV_TYPES.contains(&var.variant_type) {
            variant.push_info_integer(b"END", &[var.end as i32])?;
            let sv_length = var.sv_length().ok_or_else(|| {
                format!(
                    "{}:{}: cannot read the {} length from the genotype column '{}'",
                    var.contig, var.start, var.variant_type, var.gt
                )
            })?;
            variant.push_info_integer(b"SVLEN", &[sv_length])?;
            variant.push_info_string(b"SVTYPE", &[var.variant_type.as_bytes()])?;
        }

        variant.push_genotypes(var.gt_value(MIN_HOM_ALT_AF))?;
        variant.push_format_integer(b"AD", &var.ad_value())?;
        variant.push_format_integer(b"ADF", &var.adf_value())?;
        variant.push_format_integer(b"ADR", &var.adr_value())?;
        variant.push_format_integer(b"DP", &[var.depth])?;
        variant.push_format_float(b"AF", &[var.af_value()])?;
        variant.push_format_integer(b"HICNT", &[var.hicnt_value()])?;
        variant.push_format_float(b"QMEAN", &[var.qmean_value()])?;
        variant.push_format_float(b"MEAN_MISMATCHES", &[var.mean_mismatches_value()])?;

        writer.write(&variant)?;
        progress.record();
    }

    Ok(0)
}

/// Reads the next record, refusing it unless its segment sits in column 35 as with `--fisher`.
fn read_fisher_record<R: Read>(
    reader: &mut Reader<R>,
    record: &mut StringRecord,
) -> Result<bool, Box<dyn error::Error>> {
    if !reader.read_record(record)? {
        return Ok(false);
    }
    let contig = record.get(2).unwrap_or_default();
    let is_segment = |index: usize| {
        record
            .get(index)
            .and_then(|field| {
                field
                    .strip_prefix(contig)?
                    .strip_prefix(':')?
                    .split_once('-')
            })
            .is_some_and(|(start, end)| start.parse::<u64>().is_ok() && end.parse::<u64>().is_ok())
    };
    if is_segment(34) {
        return Ok(true);
    }
    let line = record.position().map_or(0, |position| position.line());
    let message = if is_segment(32) {
        format!(
            "Expected the strand-bias p-value and odds-ratio columns on line {line}: run VarDictJava with --fisher!"
        )
    } else {
        format!("Expected a {contig}:start-end segment in column 35 on line {line}!")
    };
    Err(message.into())
}

/// Names the 1-based line and column of a record that could not be deserialized.
fn describe_parse_error(error: csv::Error) -> String {
    match error.kind() {
        csv::ErrorKind::Deserialize {
            pos: Some(pos),
            err,
        } => match err.field() {
            Some(index) => format!(
                "Could not parse column {} on line {}: {}!",
                index + 1,
                pos.line(),
                err.kind()
            ),
            None => format!("Could not parse line {}: {}!", pos.line(), err.kind()),
        },
        _ => error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use anyhow::Result;
    use file_diff::diff;
    use pretty_assertions::assert_eq;
    use rstest::rstest;
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
            sample,
            &TumorOnly,
            false,
        )?;
        assert_eq!(exit, 0);
        assert!(diff(output.path().to_str().unwrap(), "tests/calls.vcf"));
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
            sample,
            &TumorOnly,
            false,
        );
        assert!(result.is_err());
    }

    #[rstest]
    #[case("tests/calls.no-fisher.reference.var")]
    #[case("tests/calls.no-fisher.variant.var")]
    fn test_refuses_rows_without_fisher_columns(
        #[case] path: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let input = BufReader::new(File::open(path)?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            "dna00001",
            &TumorOnly,
            false,
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Expected the strand-bias p-value and odds-ratio columns on line 1: run VarDictJava with --fisher!"
        );
        assert_eq!(output.as_file().metadata()?.len(), 0);
        Ok(())
    }

    #[test]
    fn test_malformed_row_is_an_error() -> Result<(), Box<dyn std::error::Error>> {
        let input = BufReader::new(File::open("tests/calls.malformed.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            "dna00001",
            &TumorOnly,
            false,
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Could not parse column 8 on line 2: invalid digit found in string!"
        );
        Ok(())
    }

    #[test]
    fn test_refuses_rows_without_a_segment() {
        let input = "dna00001\tBRINP3\tchr1\t190098265\t190098265\tA\tA\n".as_bytes();
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            Some(output.path().into()),
            &reference,
            "dna00001",
            &TumorOnly,
            false,
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Expected a chr1:start-end segment in column 35 on line 1!"
        );
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
            sample,
            &TumorOnly,
            true,
        )?;
        assert_eq!(exit, 0);

        // Read the output and verify no non-variant sites exist
        use rust_htslib::bcf::{Read, Reader as VcfReader};
        let mut reader = VcfReader::from_path(output.path()).expect("Error opening output file!");

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
