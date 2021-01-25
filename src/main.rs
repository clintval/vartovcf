use std::path::PathBuf;
use std::process;

use env_logger::Env;
use structopt::StructOpt;

use vartovcf;

#[derive(StructOpt, Debug)]
#[structopt(rename_all = "kebab-case", about)]
struct Opt {

    /// Input VAR file or stream
    #[structopt(short = "i", long = "--input", parse(from_os_str))]
    input: PathBuf,

    /// Output VCF file or stream
    #[structopt(short = "o", long = "--output", parse(from_os_str))]
    output: PathBuf,
}

/// Main binary entrypoint.
fn main() {
    let env: Env = Env::default().default_filter_or("info");
    env_logger::Builder::from_env(env).init();
    let opt = Opt::from_args();
    match vartovcf::run(&opt.input, &opt.output) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => panic!("{}", except),
    }
}
