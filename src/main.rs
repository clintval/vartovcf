//! Convert variants from VarDict/VarDictJava into VCF format, fast.
//! Only output from tumor-only calling is currently supported.
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process;

use anyhow::{Error, Result};
use env_logger::Env;
use log::*;
use structopt::StructOpt;

use vartovcflib::vartovcf;

#[derive(Clone, Debug, StructOpt)]
#[structopt(rename_all = "kebab-case", about)]
struct Opt {
    /// The sample name
    #[structopt(short = "s", long = "--sample")]
    sample: String,

    /// The indexed reference sequence file
    #[structopt(short = "r", long = "--reference", parse(from_os_str))]
    reference: PathBuf,

    /// Input VAR file or stream, defaults to /dev/stdin
    #[structopt(short = "i", long = "--input", parse(from_os_str))]
    input: Option<PathBuf>,

    /// Output VCF file or stream, defaults to /dev/stdout
    #[structopt(short = "o", long = "--output", parse(from_os_str))]
    output: Option<PathBuf>,
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

    match vartovcf(input, opt.output, &opt.reference, &opt.sample) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => panic!("{}", except),
    }
}
