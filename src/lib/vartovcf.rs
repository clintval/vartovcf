//! A module for running VarDict/VarDictJava data conversions.
use std::error;
use std::fmt::Debug;
use std::fs;
use std::path::{Path, PathBuf};

use csv::ReaderBuilder;
use log::*;
use rust_htslib::bcf::Format;
use rust_htslib::bcf::Writer as VcfWriter;

use crate::io;
use crate::record::tumor_only_header;
use crate::record::TumorOnlyVariant;
use rust_htslib::bcf::record::GenotypeAllele;

/// Runs the tool `vartovcf` on an input VAR file and writes the records to an output VCF file.
///
/// # Arguments
///
/// * `input` - The input VAR file or stream
/// * `output` - The output VCF file or stream
///
/// # Returns
///
/// Returns the result of the execution with an integer exit code for success (0).
pub fn run<I, O>(input: I, output: O, sample: String) -> Result<i32, Box<dyn error::Error>>
    where I: AsRef<Path> + Debug, O: AsRef<Path> + Debug {
    let input: PathBuf = fs::canonicalize(input)?;
    info!("Input file:  {:?}", input);
    info!("Output file: {:?}", output);

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(input)?;

    let header      = tumor_only_header(sample);
    let plain_text  = !io::has_gzip_ext(&output);
    let mut carry   = csv::StringRecord::new();
    let mut writer  = VcfWriter::from_path(output, &header, plain_text, Format::VCF)?;
    let mut variant = writer.empty_record();
    let mut count   = 0;

    while reader.read_record(&mut carry)? {
        let var: TumorOnlyVariant = carry.deserialize(None)?;
        let rid = writer.header().name2rid(var.contig.as_bytes()).unwrap();
        let alleles = &[GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)];
        variant.set_rid(Some(rid));
        variant.set_pos(var.start as i64 - 1);
        variant.set_alleles(&[var.ref_allele.as_bytes(), var.alt_allele.as_bytes()])?;
        variant.push_genotypes(alleles).unwrap();
        writer.write(&variant)?;
        count += 1;
    }

    info!("Processed {} variant records", count);
    Ok(0)
}
