//! A module for running VarDict/VarDictJava data conversions.
use std::error;
use std::fmt::Debug;
use std::fs;
use std::path::{Path, PathBuf};

use csv::ReaderBuilder;
use log::*;

use crate::records::TumorOnlyVariant;

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
#[must_use]
pub fn run<P>(input: P, output: P) -> Result<i32, Box<dyn error::Error>>
    where P: AsRef<Path> + Debug {
    let input: PathBuf = fs::canonicalize(input)?;
    info!("Input file:  {:?}", input);
    info!("Output file: {:?}", output);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(input)?;

    // TODO: Prepare writers over the Path, but understand they may be standard streams

    let mut count: isize = 0;
    for result in reader.deserialize() {
        let record: TumorOnlyVariant = result?;
        println!("{:?}", record); // TODO: Write to buffered writer.
        count += 1;               // TODO: Implement a ProgressLogger like in fgbio/htsjdk?
    }
    info!("Processed {} variant records", count);
    Ok(0)
}
