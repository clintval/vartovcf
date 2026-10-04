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

use vartovcflib::filter::FilterThresholds;
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

    /// Label calls NEAR_READ_END when their mean distance to the nearer read end is below this; not applied when not given
    #[arg(long, value_name = "MIN_MEAN_DIST")]
    filter_near_read_end: Option<f32>,

    /// Label calls LOW_MEAN_MAPQ when the mean mapping quality of their ALT reads is below this; not applied when not given
    #[arg(long, value_name = "MIN_MEAN_MAPQ")]
    filter_low_mean_mapq: Option<f32>,

    /// Label one-base indels HOMOPOLYMER_INDEL in homopolymers of at least this many copies; not applied when not given
    #[arg(long, value_name = "MIN_COPIES")]
    filter_homopolymer_indel: Option<f32>,

    /// Only label HOMOPOLYMER_INDEL when AF is below this; no AF limit when not given
    #[arg(long, value_name = "MAX_AF", requires = "filter_homopolymer_indel")]
    filter_homopolymer_indel_max_af: Option<f32>,

    /// Label one-unit indels TANDEM_REPEAT_INDEL in 2-6 bp tandem repeats of at least this many copies; not applied when not given
    #[arg(long, value_name = "MIN_COPIES")]
    filter_tandem_repeat_indel: Option<f32>,

    /// Only label TANDEM_REPEAT_INDEL when AF is below this; no AF limit when not given
    #[arg(long, value_name = "MAX_AF", requires = "filter_tandem_repeat_indel")]
    filter_tandem_repeat_indel_max_af: Option<f32>,

    /// Label calls HIGH_MEAN_MISMATCHES when their ALT reads average more substitution mismatches than this; not applied when not given
    #[arg(long, value_name = "MAX_MEAN_MISMATCHES")]
    filter_high_mean_mismatches: Option<f32>,
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
        &FilterThresholds {
            near_read_end: opt.filter_near_read_end,
            low_mean_mapq: opt.filter_low_mean_mapq,
            homopolymer_indel: opt.filter_homopolymer_indel,
            homopolymer_indel_max_af: opt.filter_homopolymer_indel_max_af,
            tandem_repeat_indel: opt.filter_tandem_repeat_indel,
            tandem_repeat_indel_max_af: opt.filter_tandem_repeat_indel_max_af,
            high_mean_mismatches: opt.filter_high_mean_mismatches,
        },
    ) {
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => {
            error!("{except}");
            process::exit(1)
        }
    }
}
