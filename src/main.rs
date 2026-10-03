//! Convert variants from VarDict/VarDictJava into VCF v4.2 format.
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use clap::Parser;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use env_logger::Env;
use log::*;
use strum::VariantNames;

use vartovcflib::{VarDictMode, vartovcf};

#[derive(Clone, Debug, Parser)]
#[command(version, about)]
struct Opt {
    /// The indexed FASTA reference sequence file
    #[arg(short, long)]
    reference: PathBuf,

    /// The input sample name, must match input data stream
    #[arg(short, long)]
    sample: String,

    /// Input VAR file or stream [default: /dev/stdin]
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output VCF file or stream [default: /dev/stdout]
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Variant calling mode.
    #[arg(
        short,
        long,
        default_value = "TumorOnly",
        value_parser = PossibleValuesParser::new(VarDictMode::VARIANTS).try_map(|mode| mode.parse::<VarDictMode>()),
    )]
    mode: VarDictMode,

    /// Skip non-variant sites (where ref_allele == alt_allele)
    #[arg(long)]
    skip_non_variants: bool,
}

/// Main binary entrypoint.
fn main() -> Result<(), Error> {
    let env = Env::default().default_filter_or("info");
    // Usage errors exit 1 instead of clap's default of 2.
    let opt = Opt::try_parse().unwrap_or_else(|err| {
        let _ = err.print();
        process::exit(if err.use_stderr() { 1 } else { 0 })
    });

    env_logger::Builder::from_env(env).init();

    let input: Box<dyn Read> = match &opt.input {
        Some(path) if path.to_str().unwrap() != "-" => {
            info!("Input file: {path:?}");
            Box::new(BufReader::new(File::open(path)?))
        }
        _ => {
            info!("Input stream: STDIN");
            Box::new(std::io::stdin())
        }
    };

    match &opt.output {
        Some(output) if output.to_str().unwrap() != "-" => info!("Output file: {output:?}"),
        _ => info!("Output stream: STDOUT"),
    }

    match vartovcf(
        input,
        opt.output,
        &opt.reference,
        &opt.sample,
        &opt.mode,
        opt.skip_non_variants,
    ) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => panic!("{}", except),
    }
}
