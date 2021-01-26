//! Convert variants from VarDict/VarDictJava into VCF format, fast.
use std::path::PathBuf;
use std::process;

use env_logger::Env;
use structopt::StructOpt;

use vartovcflib::io;
use vartovcflib::vartovcf;

#[derive(Clone, StructOpt, Debug)]
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
fn main() {
    let env    = Env::default().default_filter_or(vartovcflib::DEFAULT_LOG_LEVEL);
    let opt    = Opt::from_args();
    let input  = opt.input.unwrap_or_else(io::stdin);
    let output = opt.output.unwrap_or_else(io::stdout);

    env_logger::Builder::from_env(env).init();
    match vartovcf::run(&input, &output,  opt.reference, opt.sample) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except)   => panic!("{}", except),
    }
}
