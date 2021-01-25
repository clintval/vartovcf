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
pub fn run<P>(input: P, output: P) -> Result<i32, Box<dyn error::Error>>
    where P: AsRef<Path> + Debug {
    let input: PathBuf = fs::canonicalize(input)?;
    info!("Input file:  {:?}", input);
    info!("Output file: {:?}", output);
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(input)?;

    let mut count: isize = 0;
    for result in reader.deserialize() {
        let record: TumorOnlyVariant = result?;
        println!("{:?}", record);
        count += 1;
    }
    info!("Processed {} variant records", count);
    Ok(0)
}
