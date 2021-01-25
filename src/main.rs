use std::path::PathBuf;
use std::process;

use env_logger::Env;
use structopt::StructOpt;

use vartovcf::{self, io};

#[derive(Clone, StructOpt, Debug)]
#[structopt(rename_all = "kebab-case", about)]
struct Opt {

    /// Input VAR file or stream, defaults to /dev/stdin
    #[structopt(short = "i", long = "--input", parse(from_os_str))]
    input: Option<PathBuf>,

    /// Output VCF file or stream, defaults to /dev/stdout
    #[structopt(short = "o", long = "--output", parse(from_os_str))]
    output: Option<PathBuf>,
}

/// Main binary entrypoint.
fn main() {
    let env    = Env::default().default_filter_or(io::DEFAULT_LOG_LEVEL);
    let opt    = Opt::from_args();
    let input  = opt.input.unwrap_or(io::stdin());
    let output = opt.output.unwrap_or(io::stdout());

    env_logger::Builder::from_env(env).init();
    match vartovcf::run(&input, &output) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except)   => panic!("{}", except),
    }
}
