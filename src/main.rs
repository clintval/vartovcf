//! Convert variants from VarDict/VarDictJava into VCF v4.2 format.
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::process;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize};

use anyhow::{Error, Result};
use clap::Parser;
use env_logger::Env;
use log::*;
use signal_hook::consts::{SIGINT, SIGTERM};
use signal_hook::flag;

use vartovcflib::filter::FilterThresholds;
use vartovcflib::mnv::LeakRate;
use vartovcflib::{Options, vartovcf};

#[derive(Clone, Debug, Parser)]
#[command(version, about)]
struct Opt {
    /// The indexed FASTA reference VarDict was run with, which every record is checked and normalized against
    #[arg(short, long)]
    reference: PathBuf,

    /// The tumor (or only) sample name, or `tumor|normal` for both; read from the input when not given, checked against tumor-only input, and used in place of a tumor-normal input's, with a warning, when it differs
    #[arg(short, long)]
    sample: Option<String>,

    /// The matched normal sample name; read from tumor-normal input that names it as `tumor|normal` and required otherwise, and used in place of the input's, with a warning, when it differs
    #[arg(long)]
    normal_sample: Option<String>,

    /// Input VAR file or stream [default: /dev/stdin]
    #[arg(short, long)]
    input: Option<PathBuf>,

    /// Output VCF file or stream, bgzipped when it ends in .gz; a bgzipped regular file is then indexed beside it as <output>.tbi, or as <output>.csi when a reference contig is longer than the 2^29 bases a tabix index can address [default: /dev/stdout]
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Directory for temporary files: the header while it is formatted and the records sorting spills. They are removed when vartovcf exits, fails or is stopped by a SIGINT or SIGTERM, which takes effect at the next input row or the next record written; a second signal, or SIGKILL, exits at once and leaves a spillsort-* directory and an incomplete output file behind [default: $TMPDIR when set and not empty, else the system's temporary directory]
    #[arg(long, value_name = "DIR")]
    tmp_dir: Option<PathBuf>,

    /// Skip non-variant sites (where ref_allele == alt_allele)
    #[arg(long)]
    skip_non_variants: bool,

    /// Fold the partial rows VarDict writes for a multi-nucleotide variant (MNV), from reads that end, have a no-call (N) or a low-quality base inside it, back into the MNV; see the README
    #[arg(long)]
    fix_leaky_mnv: bool,

    /// Probability that a read of an MNV molecule counts toward any one partial row, at least 0.001 and below 0.5; with --fix-leaky-mnv, a partial folds when it has fewer reads (c) than its MNV (n) and P(Binomial(n + c, rate) >= c) >= 0.05. The rate varies with the library and the MNV's sequence context; see the README to estimate it
    #[arg(long, value_name = "RATE", default_value_t = LeakRate::default().value(), requires = "fix_leaky_mnv")]
    leaky_mnv_rate: f64,

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

    /// Label calls SAME_READ_POSITION when 2 or more ALT reads all place the variant at one read position, at an AF below this; not applied when not given
    #[arg(long, value_name = "MAX_AF")]
    filter_same_read_position: Option<f32>,

    /// Label calls STRAND_BIAS when the Fisher p-value of REF and ALT reads by strand is below this; not applied when not given
    #[arg(long, value_name = "MAX_P")]
    filter_strand_bias: Option<f32>,

    /// Only label STRAND_BIAS when the folded strand odds ratio is above this or the table has an empty cell; no odds-ratio limit when not given
    #[arg(long, value_name = "MIN_ODDS_RATIO", requires = "filter_strand_bias")]
    filter_strand_bias_min_odds_ratio: Option<f32>,

    /// Only label STRAND_BIAS when AF is below this; no AF limit when not given
    #[arg(long, value_name = "MAX_AF", requires = "filter_strand_bias")]
    filter_strand_bias_max_af: Option<f32>,

    /// Label calls LOW_AF when their AF is below this; not applied when not given
    #[arg(long, value_name = "MIN_AF")]
    filter_low_af: Option<f32>,

    /// Label calls LOW_QMEAN when the mean per-read variant base quality of their ALT reads is below this; not applied when not given
    #[arg(long, value_name = "MIN_QMEAN")]
    filter_low_qmean: Option<f32>,

    /// Label calls LOW_DP when their read depth is below this; not applied when not given
    #[arg(long, value_name = "MIN_DP")]
    filter_low_dp: Option<i32>,
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

    // The first SIGINT or SIGTERM stops the run cleanly, and a second exits at once.
    let stop = Arc::new(AtomicBool::new(false));
    let signal = Arc::new(AtomicUsize::new(0));
    for number in [SIGINT, SIGTERM] {
        flag::register_conditional_shutdown(number, 128 + number, Arc::clone(&stop))?;
        flag::register(number, Arc::clone(&stop))?;
        flag::register_usize(number, Arc::clone(&signal), number as usize)?;
    }

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

    let options = Options {
        output: opt.output,
        reference: opt.reference,
        sample: opt.sample.as_deref(),
        normal_sample: opt.normal_sample.as_deref(),
        skip_non_variants: opt.skip_non_variants,
        fix_leaky_mnv: opt.fix_leaky_mnv.then(|| {
            LeakRate::new(opt.leaky_mnv_rate).unwrap_or_else(|message| {
                error!("{message}");
                process::exit(1)
            })
        }),
        tmp_dir: opt.tmp_dir,
        signal: Some(&signal),
        filters: FilterThresholds {
            near_read_end: opt.filter_near_read_end,
            low_mean_mapq: opt.filter_low_mean_mapq,
            homopolymer_indel: opt.filter_homopolymer_indel,
            homopolymer_indel_max_af: opt.filter_homopolymer_indel_max_af,
            tandem_repeat_indel: opt.filter_tandem_repeat_indel,
            tandem_repeat_indel_max_af: opt.filter_tandem_repeat_indel_max_af,
            high_mean_mismatches: opt.filter_high_mean_mismatches,
            same_read_position: opt.filter_same_read_position,
            strand_bias: opt.filter_strand_bias,
            strand_bias_min_odds_ratio: opt.filter_strand_bias_min_odds_ratio,
            strand_bias_max_af: opt.filter_strand_bias_max_af,
            low_af: opt.filter_low_af,
            low_qmean: opt.filter_low_qmean,
            low_dp: opt.filter_low_dp,
        },
    };

    let to_file = options
        .output
        .as_ref()
        .is_some_and(|path| path.as_os_str() != "-");
    match vartovcf(input, &options) {
        Ok(exit_code) if exit_code > 128 => {
            let output = if to_file {
                "removed the output file"
            } else {
                "output already written may end in a partial record"
            };
            error!("Stopped by signal {}; {output}.", exit_code - 128);
            process::exit(exit_code)
        }
        Ok(exit_code) => process::exit(exit_code),
        Err(except) => {
            error!("{except}");
            process::exit(1)
        }
    }
}
