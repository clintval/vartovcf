//! Convert variants from VarDict/VarDictJava into VCF v4.2 format.
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use env_logger::Env;
use log::*;
use structopt::StructOpt;
use strum::VariantNames;

use vartovcflib::{vartovcf, VarDictMode};

#[derive(Clone, Debug, StructOpt)]
#[structopt(
    setting = structopt::clap::AppSettings::ColoredHelp,
    setting = structopt::clap::AppSettings::DeriveDisplayOrder,
    rename_all = "kebab-case",
    about
)]
struct Opt {
    /// The indexed FASTA reference sequence file
    #[structopt(short = "r", long = "--reference", parse(from_os_str))]
    reference: PathBuf,

    /// The input sample name, must match input data stream
    #[structopt(short = "s", long = "--sample")]
    sample: String,

    /// Input VAR file or stream [default: /dev/stdin]
    #[structopt(short = "i", long = "--input", parse(from_os_str))]
    input: Option<PathBuf>,

    /// Output VCF file or stream [default: /dev/stdout]
    #[structopt(short = "o", long = "--output", parse(from_os_str))]
    output: Option<PathBuf>,

    /// Variant calling mode.
    #[structopt(short = "m", long = "--mode", default_value = "TumorOnly", possible_values = &VarDictMode::VARIANTS)]
    mode: VarDictMode,
}

/// Main binary entrypoint.
#[cfg(not(tarpaulin_include))]
fn main() -> Result<(), Error> {
    let env = Env::default().default_filter_or(vartovcflib::DEFAULT_LOG_LEVEL);
    let opt = Opt::from_args();

    env_logger::Builder::from_env(env).init();

    let input: Box<dyn Read> = match &opt.input {
        Some(path) if path.to_str().unwrap() != "-" => {
            info!("Input file: {:?}", path);
            Box::new(BufReader::new(File::open(path)?))
        }
        _ => {
            info!("Input stream: STDIN");
            Box::new(std::io::stdin())
        }
    };

    match &opt.output {
        Some(output) if output.to_str().unwrap() != "-" => info!("Output file: {:?}", output),
        _ => info!("Output stream: STDOUT"),
    }

    match vartovcf(input, opt.output, &opt.reference, &opt.sample, &opt.mode) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => panic!("{}", except),
    }
}
