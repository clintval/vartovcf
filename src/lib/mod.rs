//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use csv::{Reader, ReaderBuilder, StringRecord};
use log::{info, warn};
use proglog::ProgLogBuilder;
use std::error;
use std::fmt::Debug;
use std::io::{Read, Write};
use std::ops::Deref;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::fai::{fasta_contigs_to_vcf_header, fasta_path_to_vcf_header};
use crate::fill::{Bed, BedLayout, Fill, parse_segment, read_bed};
use crate::filter::FilterThresholds;
use crate::mnv::{Fold, LeakRate, LeakyMnvFix, Merged, Rule, fold_partials};
use crate::normalize::{Normalized, normalize};
use crate::output::{Index, Output};
use crate::reads::ReadSources;
use crate::record::MIN_HOM_ALT_AF;
use crate::record::TumorNormalVariant;
use crate::record::TumorOnlyVariant;
use crate::record::{tumor_normal_fisher_p_value, tumor_normal_header, tumor_only_header};
use crate::reference::Reference;
use crate::site::{Counts, PerSample, SAMPLES, Walk};
use crate::sort::{Entry, Rank, RecordSorter, SortKey, resolve_tmp_dir};
use crate::vcf::{GenotypeAllele, Numeric, Record, RecordWriter, header_bytes};

pub mod fai;
pub mod fill;
pub mod filter;
pub mod fisher;
pub mod hicnt;
pub mod io;
pub mod mnv;
pub mod normalize;
pub mod output;
pub mod reads;
pub mod record;
pub mod reference;
pub mod site;
pub mod sort;
pub mod vcf;

/// The structural variant types VarDict writes, which are the values of the `SVTYPE` INFO field.
pub const VALID_SV_TYPES: &[&str] = &["DEL", "DUP", "INV"];

/// Namespace for path parts and extensions.
pub mod path {

    /// The Gzip extension.
    pub const GZIP_EXTENSION: &str = "gz";
}

/// The layouts of VarDictJava output rows `vartovcf` reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// One sample per row, from `vardict-java -b sample.bam`.
    TumorOnly,
    /// A tumor and its matched normal per row, from `vardict-java -b "tumor.bam|normal.bam"`.
    TumorNormal,
}

impl Layout {
    /// The 0-based column holding VarDict's `chr:start-end` segment in a row without `-fisher`.
    fn segment_column(&self) -> usize {
        match self {
            Layout::TumorOnly => 32,
            Layout::TumorNormal => 48,
        }
    }

    /// The number of samples in a row of this layout.
    fn samples(&self) -> usize {
        match self {
            Layout::TumorOnly => 1,
            Layout::TumorNormal => 2,
        }
    }

    /// The fewest columns a row of this layout has without `-fisher`.
    fn min_columns(&self) -> usize {
        match self {
            Layout::TumorOnly => 36,
            Layout::TumorNormal => 55,
        }
    }
}

/// A row's layout and whether VarDictJava ran with `-fisher`, whose strand-bias and tumor-normal
/// p-value and odds-ratio columns are dropped as each row is read.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Columns {
    layout: Layout,
    fisher: bool,
}

impl Columns {
    /// Every layout, without and with `-fisher`.
    const ALL: [Columns; 4] = [
        Columns::new(Layout::TumorOnly, false),
        Columns::new(Layout::TumorOnly, true),
        Columns::new(Layout::TumorNormal, false),
        Columns::new(Layout::TumorNormal, true),
    ];

    const fn new(layout: Layout, fisher: bool) -> Self {
        Columns { layout, fisher }
    }

    /// The 0-based columns `-fisher` adds, in order.
    fn fisher_columns(&self) -> &'static [usize] {
        match (self.layout, self.fisher) {
            (_, false) => &[],
            (Layout::TumorOnly, true) => &[20, 21],
            (Layout::TumorNormal, true) => &[25, 26, 45, 46, 59, 60],
        }
    }

    /// Maps a 0-based column of the row without `-fisher`'s columns to the column as read.
    fn input_column(&self, column: usize) -> usize {
        self.fisher_columns().iter().fold(column, |column, &added| {
            column + usize::from(column >= added)
        })
    }

    /// The 0-based column holding VarDict's `chr:start-end` segment.
    fn segment_column(&self) -> usize {
        self.input_column(self.layout.segment_column())
    }

    /// The fewest columns a row has.
    fn min_columns(&self) -> usize {
        self.layout.min_columns() + self.fisher_columns().len()
    }

    /// Copies a row and its position into `copy`, without the columns `-fisher` adds.
    fn copy_row(&self, record: &StringRecord, copy: &mut StringRecord) {
        let fisher_columns = self.fisher_columns();
        copy.clear();
        for (column, field) in record.iter().enumerate() {
            if !fisher_columns.contains(&column) {
                copy.push_field(field);
            }
        }
        copy.set_position(record.position().cloned());
    }

    /// Drops the columns `-fisher` adds, rebuilding the row in `scratch` to spare an allocation.
    fn drop_fisher_columns(&self, record: &mut StringRecord, scratch: &mut StringRecord) {
        if self.fisher {
            self.copy_row(record, scratch);
            std::mem::swap(record, scratch);
        }
    }
}

/// What `vartovcf` reads and writes, and how it converts the rows.
#[derive(Clone, Debug, Default)]
pub struct Options<'a> {
    /// The output VCF file, or standard output when `None`.
    pub output: Option<PathBuf>,
    /// The indexed FASTA reference VarDict was run with.
    pub reference: PathBuf,
    /// The tumor (or only) sample name, read from the input when not given.
    pub sample: Option<&'a str>,
    /// The matched normal sample name, required for tumor-normal input that does not name it.
    pub normal_sample: Option<&'a str>,
    /// Skip non-variant sites (where ref_allele == alt_allele).
    pub skip_non_variants: bool,
    /// Write INFO CONTEXT_3mer, the reference bases at POS-1, POS and POS+1, on every record.
    pub context_3mer: bool,
    /// Write a reference-only record at each position of VarDict's regions without one; see [`fill`].
    pub fill_gvcf: bool,
    /// The BED file VarDict was given, whose regions VarDict wrote no row for are filled too.
    pub fill_gvcf_bed: Option<PathBuf>,
    /// How VarDict turned the BED file's lines into its regions.
    pub fill_gvcf_bed_layout: BedLayout,
    /// How to fold the partial rows of each MNV back into it, by their counts at a leak rate or
    /// by their reads, or `None` to not fold; see [`mnv`].
    pub fix_leaky_mnv: Option<LeakyMnvFix>,
    /// The FILTER labels to apply.
    pub filters: FilterThresholds,
    /// The directory to spill records to while sorting, resolved by [`sort::resolve_tmp_dir`].
    pub tmp_dir: Option<PathBuf>,
    /// Set to a signal's number when one asks the run to stop.
    pub signal: Option<&'a AtomicUsize>,
}

/// Runs the tool `vartovcf` on an input VAR file or stream and writes the records to a VCF file.
///
/// Each VarDict region's rows are parsed and, when asked, their leaky MNV partials folded, both in
/// VarDict's own coordinates; each record is then normalized against the reference and sorted, and
/// each call is written once, after the input ends.
///
/// Returns the exit code: 0 for success, also when standard output's reader closes it early, or 128
/// plus the signal that stopped the run, after removing its temporary files, output file and index.
pub fn vartovcf<I: Read>(input: I, options: &Options) -> Result<i32, Box<dyn error::Error>> {
    let tmp_dir = resolve_tmp_dir(options.tmp_dir.as_deref())?;
    let mut output = Output::create(options.output.as_deref())?;
    let signal = || {
        options
            .signal
            .map_or(0, |signal| signal.load(Ordering::Relaxed))
    };
    match convert(input, options, &tmp_dir, &mut output, &signal) {
        Ok(Some(index)) => {
            if output.finish(index, &|| signal() != 0)? {
                Ok(0)
            } else {
                Ok(128 + signal() as i32)
            }
        }
        Ok(None) => {
            output.discard();
            Ok(128 + signal() as i32)
        }
        Err(_) if output.is_closed_pipe() => {
            output.discard();
            Ok(0)
        }
        Err(error) => {
            output.discard();
            Err(error)
        }
    }
}

/// Converts the input's rows into sorted VCF text in the output, returning the index a bgzipped
/// output of the reference needs, or `None` when a signal stopped the run first.
fn convert<I: Read>(
    input: I,
    options: &Options,
    tmp_dir: &Path,
    output: &mut Output,
    signal: &dyn Fn() -> usize,
) -> Result<Option<Index>, Box<dyn error::Error>> {
    if options.fill_gvcf_bed.is_some() && !options.fill_gvcf {
        return Err("--fill-gvcf-bed requires --fill-gvcf!".into());
    }
    if options.fill_gvcf && options.skip_non_variants {
        return Err("--fill-gvcf cannot be used with --skip-non-variants!".into());
    }
    let reference = Reference::from_path(&options.reference)?;
    let index = Index::for_longest_contig(reference.longest_contig_length());
    let bed = match &options.fill_gvcf_bed {
        Some(path) => {
            let lines = read_bed(path, |contig| reference.length(contig).ok())?;
            Some(Bed::new(path, options.fill_gvcf_bed_layout, lines))
        }
        None => None,
    };

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_reader(input);

    // Peek at the first row to learn the layout and sample before the header is made.
    let mut carry = StringRecord::new();
    let mut scratch = StringRecord::new();
    let mut peeked = reader.read_record(&mut carry)?;
    let columns = if peeked {
        let columns = detect_columns(&carry)?;
        check_width(&carry, columns)?;
        Some(columns)
    } else {
        None
    };
    let row_sample = peeked.then(|| carry.get(0).unwrap_or_default().to_string());
    let (sample, normal_sample) = resolve_samples(
        row_sample.as_deref(),
        columns.map(|columns| columns.layout),
        options.sample,
        options.normal_sample,
    )?;
    let columns = columns.unwrap_or(Columns::new(Layout::TumorOnly, false));
    let layout = columns.layout;
    if options.fill_gvcf && layout == Layout::TumorNormal {
        return Err("--fill-gvcf needs tumor-only input: VarDict writes no reference-only rows for a tumor-normal pair, so a position without a variant cannot be told from one without reads!".into());
    }

    let filters = &options.filters;
    let mut header = match &normal_sample {
        Some(normal_sample) => tumor_normal_header(
            &sample,
            normal_sample,
            filters,
            options.fix_leaky_mnv.as_ref(),
            options.context_3mer,
        ),
        None => tumor_only_header(
            &sample,
            filters,
            options.fix_leaky_mnv.as_ref(),
            options.context_3mer,
            options.fill_gvcf,
        ),
    };

    fasta_contigs_to_vcf_header(&options.reference, &mut header);
    fasta_path_to_vcf_header(&options.reference, &mut header)
        .expect("Adding FASTA path to header failed!");
    let header_text = header_bytes(&header)?;

    let folder = match &options.fix_leaky_mnv {
        None => None,
        Some(LeakyMnvFix::Counts(leak_rate)) => Some(Folder::Counts(*leak_rate)),
        Some(LeakyMnvFix::Reads { bams, options }) => {
            let sources = ReadSources::open(bams, *options)?;
            if peeked && sources.samples() != layout.samples() {
                return Err(format!(
                    "--reads names {} BAM(s) for {} sample(s): give one per sample, the tumor (or only) sample's first!",
                    sources.samples(),
                    layout.samples()
                )
                .into());
            }
            let names: Vec<&str> = std::iter::once(sample.as_str())
                .chain(normal_sample.as_deref())
                .collect();
            for warning in sources.sample_warnings(&names) {
                warn!("{warning}");
            }
            Some(Folder::Reads(sources))
        }
    };
    let mut sink = Sink {
        records: RecordWriter::new(header),
        reference,
        sorter: RecordSorter::new(tmp_dir)?,
        line: Vec::new(),
        seen: SeenRows::default(),
        filters,
        context_3mer: options.context_3mer,
        folder,
        folds: FoldCounts::default(),
        bed,
        counts: Vec::new(),
        progress: ProgLogBuilder::new()
            .name("main")
            .verb("Processed")
            .noun("variant records")
            .unit(100_000)
            .build(),
    };

    sink.enter_bed()?;

    let row_sample = row_sample.unwrap_or_default();
    let column = columns.segment_column();
    let mut segment = String::new();
    let mut ordinal: u64 = 0;
    let mut region: Vec<StringRecord> = Vec::new();

    while std::mem::take(&mut peeked) || read_row(&mut reader, &mut carry, columns)? {
        if signal() != 0 {
            return Ok(None);
        }
        let field = |index: usize| carry.get(index).unwrap_or_default();
        if field(column) != segment {
            sink.add_region(&region, columns, ordinal)?;
            region.clear();
            ordinal += 1;
            segment.clear();
            segment.push_str(field(column));
            sink.enter_region(&segment)?;
        }

        if field(5).is_empty() {
            continue; // A row without a REF allele is one not to deserialize.
        }
        if options.skip_non_variants && field(5).eq_ignore_ascii_case(field(6)) {
            columns.drop_fisher_columns(&mut carry, &mut scratch);
            let (calls, _) = parse_row(&carry, columns)?;
            sink.add_counts(&calls[..layout.samples()], ordinal)?;
            continue;
        }

        if field(0) != row_sample {
            let message = format!("Expected sample '{}' found '{}'!", row_sample, field(0));
            return Err(message.into());
        };
        if options.fix_leaky_mnv.is_some() {
            region.push(compact(&carry, columns));
        } else {
            columns.drop_fisher_columns(&mut carry, &mut scratch);
            let (calls, status) = parse_row(&carry, columns)?;
            sink.add(&calls[..layout.samples()], status, ordinal)?;
        }
    }
    sink.add_region(&region, columns, ordinal)?;
    if let Some(Folder::Reads(sources)) = &sink.folder {
        let summary = sources.summary;
        info!(
            "Reads decided {} partial row(s) and sample(s) under {} MNV pileup(s); {} reconciled with VarDict's ALT reads, and {} of the {} whose reads would have moved did not.",
            summary.decisions,
            summary.pileups,
            summary.reconciled,
            summary.decisive_unreconciled,
            summary.decisive
        );
        if summary.empty * 2 > summary.decisions {
            warn!(
                "No read in the reads fed {} of the {} partial row(s) and sample(s) VarDict counted reads toward. Check that --reads names the BAMs VarDict was run on, in sample order, and that --reads-base-quality, --reads-min-mapq, --reads-exclude-flags, --reads-max-mismatches and --reads-mnv-gap match the run.",
                summary.empty, summary.decisions
            );
        } else if summary.decisive >= 2 && summary.decisive_unreconciled * 2 > summary.decisive {
            warn!(
                "{} of the {} partial row(s) and sample(s) whose reads would have moved did not find the reads VarDict counted, so they kept every read. Check that --reads names the BAMs VarDict was run on, in sample order, and that --reads-base-quality, --reads-min-mapq, --reads-exclude-flags, --reads-max-mismatches and --reads-mnv-gap match the run.",
                summary.decisive_unreconciled, summary.decisive
            );
        }
    }
    if sink.folder.is_some() {
        let folds = sink.folds;
        info!(
            "Folded {} partial row(s) whole and {} in part into {} MNV(s), moving {} ALT read(s).",
            folds.whole, folds.split, folds.mnvs, folds.reads
        );
    }
    output
        .write_all(&header_text)
        .map_err(|error| format!("Could not write the output: {error}!"))?;
    let Some(fills) = sink.write(output, options.fill_gvcf, &|| signal() != 0)? else {
        return Ok(None);
    };
    if options.fill_gvcf {
        info!(
            "Filled {} position(s) without a reference-only record: {} from VarDict's rows, {} inside a deletion, MNV or Complex call and {} with no reads.",
            fills.rows + fills.spanned + fills.no_reads,
            fills.rows,
            fills.spanned,
            fills.no_reads
        );
    }
    Ok(Some(index))
}

/// Parses a row into its calls, one per sample with the tumor (or only) sample first and a default
/// call filling the second slot of a tumor-only row, and its tumor-normal VARDICT_STATUS.
fn parse_row(
    record: &StringRecord,
    columns: Columns,
) -> Result<([TumorOnlyVariant<'_>; 2], Option<&str>), String> {
    let describe = |error| describe_parse_error(error, columns);
    Ok(match columns.layout {
        Layout::TumorOnly => {
            let call: TumorOnlyVariant = record.deserialize(None).map_err(describe)?;
            ([call, TumorOnlyVariant::default()], None)
        }
        Layout::TumorNormal => {
            let row: TumorNormalVariant = record.deserialize(None).map_err(describe)?;
            ([row.tumor_call(), row.normal_call()], Some(row.status))
        }
    })
}

/// Normalizes records and adds them to the sort, which writes each call once.
struct Sink<'r> {
    records: RecordWriter,
    reference: Reference,
    sorter: RecordSorter,
    line: Vec<u8>,
    seen: SeenRows,
    filters: &'r FilterThresholds,
    context_3mer: bool,
    folder: Option<Folder>,
    folds: FoldCounts,
    bed: Option<Bed>,
    counts: Vec<u8>,
    progress: proglog::ProgLog,
}

/// How many positions were filled, by how their depth was found.
#[derive(Clone, Copy, Debug, Default)]
struct FillCounts {
    rows: usize,
    spanned: usize,
    no_reads: usize,
}

/// How many rows the folds touched.
#[derive(Clone, Copy, Debug, Default)]
struct FoldCounts {
    mnvs: usize,
    whole: usize,
    split: usize,
    reads: i64,
}

/// How a sink folds leaky MNV partials.
enum Folder {
    Counts(LeakRate),
    Reads(ReadSources),
}

impl Sink<'_> {
    /// Adds a record already normalized, from the region with this ordinal, to the sort, ranked
    /// against other records of its call by its calls and by the row it was written for: the first
    /// region's records first, then the best supported, by [`Rank`], so the sort writes the one each
    /// call keeps. Returns whether the writer kept it.
    fn add_record(
        &mut self,
        calls: &[TumorOnlyVariant],
        normalized: &Normalized,
        status: Option<&str>,
        extras: Extras,
        region: u64,
        row: &TumorOnlyVariant,
    ) -> Result<bool, Box<dyn error::Error>> {
        let written: Vec<TumorOnlyVariant> =
            calls.iter().map(|call| normalized.apply(call)).collect();
        let context = self.context(&written[0])?;
        let extras = Extras {
            context: context.as_deref(),
            ..extras
        };
        let Some(record) = build_record(&mut self.records, &written, status, self.filters, extras)?
        else {
            return Ok(false);
        };
        let key = sort_key(record.rid(), normalized);
        self.line.clear();
        record.write(&mut self.line)?;
        let rank = Rank::new(region, calls, row);
        self.sorter.push(key, rank, &self.line)?;
        self.progress.record();
        Ok(true)
    }

    /// The reference's 3-mer around a normalized call's POS, when asked for.
    fn context(&mut self, call: &TumorOnlyVariant) -> Result<Option<String>, String> {
        self.context_3mer
            .then(|| self.reference.context(call.contig, call.start, 1))
            .transpose()
    }

    /// Adds the rows of one VarDict region, a block it writes together for one target, with the
    /// partial rows of each MNV folded into it; VarDict orders a tumor-only region's positions by
    /// hash, so a partial row can come before or after its MNV.
    ///
    /// A substitution an earlier region added, or folded as a partial, is not added again, and is
    /// not folded into an MNV here, so its reads are never counted twice; a partial an MNV an earlier
    /// region already wrote takes here is dropped, as that MNV is not written again. The set of
    /// substitutions seen is kept for this, and only for substitutions, so it grows with the run's
    /// substitution rows, not its rows.
    fn add_region(
        &mut self,
        region: &[StringRecord],
        columns: Columns,
        ordinal: u64,
    ) -> Result<(), Box<dyn error::Error>> {
        let samples = columns.layout.samples();
        let mut parsed: Vec<TumorOnlyVariant> = Vec::with_capacity(region.len() * samples);
        let mut statuses: Vec<Option<&str>> = Vec::with_capacity(region.len());
        for record in region {
            let (row_calls, status) = parse_row(record, columns)?;
            parsed.extend(row_calls.into_iter().take(samples));
            statuses.push(status);
        }
        let calls: Vec<&[TumorOnlyVariant]> = parsed.chunks(samples).collect();
        let normalized = calls
            .iter()
            .map(|row| normalize(&mut self.reference, &row[0]))
            .collect::<Result<Vec<_>, _>>()?;
        let written: Vec<bool> = calls
            .iter()
            .zip(&normalized)
            .map(|(row, call)| self.seen.contains(row[0].contig, call))
            .collect();
        let folds = match &mut self.folder {
            None => vec![Fold::Keep; calls.len()],
            Some(Folder::Counts(leak_rate)) => {
                fold_partials(&calls, &written, Rule::Counts(*leak_rate))
            }
            Some(Folder::Reads(sources)) => {
                let fates = sources.fates(&calls, &written, &mut self.reference)?;
                fold_partials(&calls, &written, Rule::Reads(&fates))
            }
        };
        let missing = [i32::missing(); SAMPLES];
        let extras = Extras {
            site_dp: &missing[..samples],
            ..Extras::default()
        };
        for (row, fold) in folds.into_iter().enumerate() {
            match &fold {
                Fold::Keep => {}
                Fold::Into(_) => self.folds.whole += 1,
                Fold::Split { .. } => self.folds.split += 1,
                Fold::Merged(merged) => {
                    self.folds.mnvs += 1;
                    self.folds.reads += merged
                        .partial_reads
                        .iter()
                        .map(|&r| i64::from(r))
                        .sum::<i64>();
                }
            }
            let call = &normalized[row];
            match fold {
                _ if written[row] => {}
                Fold::Keep => {
                    self.add_record(
                        calls[row],
                        call,
                        statuses[row],
                        extras,
                        ordinal,
                        &calls[row][0],
                    )?;
                }
                Fold::Merged(merged) => {
                    self.add_record(
                        &merged.calls,
                        call,
                        statuses[row],
                        Extras {
                            merged: Some(&merged),
                            ..extras
                        },
                        ordinal,
                        &calls[row][0],
                    )?;
                }
                Fold::Split {
                    calls: kept, moved, ..
                } => {
                    self.add_record(
                        &kept,
                        call,
                        statuses[row],
                        Extras {
                            moved: Some(&moved),
                            ..extras
                        },
                        ordinal,
                        &kept[0],
                    )?;
                }
                Fold::Into(_) => {}
            }
        }
        for (row, call) in calls.iter().zip(&normalized) {
            if is_substitution(call) {
                self.seen.insert(row[0].contig, call)?;
            }
        }
        for row in &calls {
            self.add_counts(row, ordinal)?;
        }
        Ok(())
    }

    /// Normalizes a row's record, from the region with this ordinal, and adds it to the sort with
    /// what VarDict counted in the row.
    fn add(
        &mut self,
        calls: &[TumorOnlyVariant],
        status: Option<&str>,
        region: u64,
    ) -> Result<(), Box<dyn error::Error>> {
        let normalized = normalize(&mut self.reference, &calls[0])?;
        let missing = [i32::missing(); SAMPLES];
        let extras = Extras {
            site_dp: &missing[..calls.len()],
            ..Extras::default()
        };
        self.add_record(calls, &normalized, status, extras, region, &calls[0])?;
        self.add_counts(calls, region)
    }

    /// Adds what VarDict counted in a row, its calls one per sample, of the region with this
    /// ordinal to the sort, under the position VarDict counted it at.
    fn add_counts(
        &mut self,
        row: &[TumorOnlyVariant],
        region: u64,
    ) -> Result<(), Box<dyn error::Error>> {
        let (position, counts) = Counts::of(row);
        if position == 0 {
            return Ok(());
        }
        let rid = self.records.rid(row[0].contig)?;
        self.counts.clear();
        counts.encode(&mut self.counts);
        self.sorter
            .push_counts(rid, position, region, &self.counts)?;
        Ok(())
    }

    /// Adds a VarDict region, from its `chr:start-end` column, to the sort, clipped to its contig,
    /// after checking it against the BED file of `--fill-gvcf-bed`.
    fn enter_region(&mut self, segment: &str) -> Result<(), Box<dyn error::Error>> {
        let (contig, start, end) = parse_segment(segment).ok_or_else(|| {
            format!("Expected VarDict's region as chr:start-end, found '{segment}'!")
        })?;
        if let Some(bed) = &self.bed {
            bed.check(contig, start, end)?;
        }
        let length = self.reference.length(contig)?;
        let rid = self.records.rid(contig)?;
        let (start, end) = (start.max(1), end.min(length));
        if start <= end {
            self.sorter.push_region(rid, start, end)?;
        }
        Ok(())
    }

    /// Adds the region of each line of the BED file of `--fill-gvcf-bed` to the sort, clipped to
    /// its contig.
    fn enter_bed(&mut self) -> Result<(), Box<dyn error::Error>> {
        let Some(bed) = &self.bed else {
            return Ok(());
        };
        for (contig, start, end) in bed.regions() {
            let rid = self.records.rid(contig)?;
            let length = self.reference.length(contig)? as i64;
            let (start, end) = (start.max(1), end.min(length));
            if start <= end {
                self.sorter.push_region(rid, start as u64, end as u64)?;
            }
        }
        Ok(())
    }

    /// Writes the sorted records to the output, returning how many positions were filled, or `None`
    /// when `stop` asked the run to stop first.
    fn write(
        self,
        output: &mut Output,
        fill_gvcf: bool,
        stop: &dyn Fn() -> bool,
    ) -> Result<Option<FillCounts>, Box<dyn error::Error>> {
        let mut entries = self.sorter.entries()?;
        let mut writer = Writer {
            records: self.records,
            reference: self.reference,
            filters: self.filters,
            context_3mer: self.context_3mer,
            fill_gvcf,
            output,
            stop,
            line: Vec::new(),
            fills: FillCounts::default(),
            contig: String::new(),
            rid: None,
            walk: Walk::default(),
            end: 0,
            position: 0,
            counts_region: None,
            wrote_record: false,
            read_record: false,
            read_reference: false,
            over_depth: 0,
            first_over_depth: None,
        };
        while let Some(entry) = entries.next_entry()? {
            if stop() || !writer.entry(entry)? {
                return Ok(None);
            }
        }
        if !writer.leave_contig()? {
            return Ok(None);
        }
        if let Some(first) = &writer.first_over_depth {
            warn!(
                "The reads VarDict counted at {} position(s) add up to more than VarDict's DP there without a row labeled Complex to explain it, first at {first}; SITE_DP keeps their sum.",
                writer.over_depth
            );
        }
        if fill_gvcf && writer.read_record && !writer.read_reference {
            warn!(
                "No reference-only row was read: --fill-gvcf expects VarDict's pileup (-p) output, which writes one at each position without an allele, and filled every other position as if VarDict had counted no reads there."
            );
        }
        Ok(Some(writer.fills))
    }
}

/// Writes the sorted entries as VCF records in one pass, holding only what VarDict counted at the
/// position it is at and the deletions and MNVs spanning it. Each record gets FORMAT SITE_DP at its
/// POS, missing outside every VarDict region, and with `--fill-gvcf` each position of a region
/// without a reference-only record gets one. Where regions overlap, a position's counts are taken
/// from the first region that counted a row there.
struct Writer<'w> {
    records: RecordWriter,
    reference: Reference,
    filters: &'w FilterThresholds,
    context_3mer: bool,
    fill_gvcf: bool,
    output: &'w mut Output,
    stop: &'w dyn Fn() -> bool,
    line: Vec<u8>,
    fills: FillCounts,
    contig: String,
    rid: Option<u32>,
    walk: Walk,
    over_depth: usize,
    first_over_depth: Option<String>,
    end: u64,
    position: u64,
    counts_region: Option<u64>,
    wrote_record: bool,
    read_record: bool,
    read_reference: bool,
}

impl Writer<'_> {
    /// Writes what an entry adds, first leaving the positions and contig before it; returns false
    /// when `stop` asked the run to stop.
    fn entry(&mut self, entry: Entry) -> Result<bool, Box<dyn error::Error>> {
        let (rid, position) = entry.locus();
        if self.rid != Some(rid) {
            if !self.leave_contig()? {
                return Ok(false);
            }
            self.enter_contig(rid)?;
        }
        if position > self.position {
            self.leave_position()?;
            if !self.fill_gap(self.position + 1, position - 1)? {
                return Ok(false);
            }
            self.walk.at(position);
            self.position = position;
            self.counts_region = None;
            self.wrote_record = false;
        }
        match entry {
            Entry::Region { end, .. } => self.end = self.end.max(end),
            Entry::Counts { region, counts, .. } => {
                if self.counts_region.is_none_or(|first| first == region) {
                    self.counts_region = Some(region);
                    self.walk.add(&Counts::decode(counts));
                }
            }
            Entry::Record {
                reference_only,
                line,
                ..
            } => {
                let site = self.walk.site();
                let counted = site.fill == Fill::Rows;
                if !self.wrote_record && !reference_only && counted {
                    self.fill(true)?;
                }
                self.wrote_record = true;
                self.read_reference |= reference_only;
                self.read_record = true;
                let site_dp: PerSample = if counted {
                    site.site_dp
                } else {
                    [i32::missing(); SAMPLES]
                };
                stamp(line, &site_dp, &mut self.line);
                self.write_line()?;
            }
        }
        Ok(true)
    }

    /// Starts the contig with index `rid`.
    fn enter_contig(&mut self, rid: u32) -> Result<(), Box<dyn error::Error>> {
        let (name, _) = self
            .records
            .header()
            .contigs()
            .get_index(rid as usize)
            .ok_or_else(|| format!("No contig has index {rid} in the VCF header!"))?;
        self.contig.clear();
        self.contig.push_str(name);
        self.rid = Some(rid);
        self.walk = Walk::default();
        self.end = 0;
        self.position = 0;
        Ok(())
    }

    /// Finishes the contig: the position it is at, then every position of its regions after it;
    /// returns false when `stop` asked the run to stop.
    fn leave_contig(&mut self) -> Result<bool, Box<dyn error::Error>> {
        if self.rid.is_none() {
            return Ok(true);
        }
        self.leave_position()?;
        self.fill_gap(self.position + 1, self.end)
    }

    /// Finishes the position it is at, filling it when no record was written there.
    fn leave_position(&mut self) -> Result<(), Box<dyn error::Error>> {
        if self.position == 0 {
            return Ok(());
        }
        if self.walk.site().exceeds_depth {
            self.over_depth += 1;
            if self.first_over_depth.is_none() {
                self.first_over_depth = Some(format!("{}:{}", self.contig, self.position));
            }
        }
        if !self.wrote_record {
            self.fill(false)?;
        }
        Ok(())
    }

    /// Fills each position from `from` to `to`, none of which has an entry, inside the regions
    /// entered; returns false when `stop` asked the run to stop.
    fn fill_gap(&mut self, from: u64, to: u64) -> Result<bool, Box<dyn error::Error>> {
        if !self.fill_gvcf {
            return Ok(true);
        }
        for position in from..=to.min(self.end) {
            if (self.stop)() {
                return Ok(false);
            }
            self.walk.at(position);
            self.position = position;
            self.fill(false)?;
        }
        Ok(true)
    }

    /// Writes a reference-only record at the position it is at, from what VarDict counted there,
    /// when filling and the position lies in a region entered; a record before ALT records at the
    /// position carries depth only, so its genotype is missing.
    fn fill(&mut self, before_alt: bool) -> Result<(), Box<dyn error::Error>> {
        if !self.fill_gvcf || self.position > self.end {
            return Ok(());
        }
        let site = self.walk.site();
        let position = self.position;
        let base = [self.reference.base(&self.contig, position - 1)?];
        let base = std::str::from_utf8(&base).map_err(|_| {
            format!(
                "Expected a reference base at {}:{position}, found byte {}!",
                self.contig, base[0]
            )
        })?;
        let context = if self.context_3mer {
            Some(self.reference.context(&self.contig, position, 1)?)
        } else {
            None
        };
        let call = TumorOnlyVariant {
            contig: &self.contig,
            start: position,
            end: position,
            ref_allele: base,
            alt_allele: base,
            depth: site.depth[0],
            ref_forward: site.ref_forward[0],
            ref_reverse: site.ref_reverse[0],
            ..TumorOnlyVariant::default()
        };
        let extras = Extras {
            context: context.as_deref(),
            site_dp: &site.site_dp[..1],
            unknown_genotype: before_alt,
            ..Extras::default()
        };
        let calls = std::slice::from_ref(&call);
        let Some(mut record) = build_record(&mut self.records, calls, None, self.filters, extras)?
        else {
            return Ok(());
        };
        record.push_info_string("GVCF_FILL", site.fill.as_str())?;
        self.line.clear();
        record.write(&mut self.line)?;
        self.write_line()?;
        match site.fill {
            Fill::Rows => self.fills.rows += 1,
            Fill::Spanned => self.fills.spanned += 1,
            Fill::NoReads => self.fills.no_reads += 1,
        }
        Ok(())
    }

    /// Writes the line held to the output.
    fn write_line(&mut self) -> Result<(), String> {
        self.output
            .write_all(&self.line)
            .map_err(|error| format!("Could not write the output: {error}!"))
    }
}

/// Copies a record's VCF line into `out` with each sample's FORMAT SITE_DP value replaced by the
/// sample's of `site_dp`, written `.` when missing.
fn stamp(line: &[u8], site_dp: &[i32], out: &mut Vec<u8>) {
    out.clear();
    let body = line.strip_suffix(b"\n").unwrap_or(line);
    let mut index = None;
    for (column, field) in body.split(|&byte| byte == b'\t').enumerate() {
        if column > 0 {
            out.push(b'\t');
        }
        match (column, index) {
            (8, _) => {
                index = field
                    .split(|&byte| byte == b':')
                    .position(|key| key == b"SITE_DP");
                out.extend_from_slice(field);
            }
            (9.., Some(index)) => {
                for (place, value) in field.split(|&byte| byte == b':').enumerate() {
                    if place > 0 {
                        out.push(b':');
                    }
                    if place != index {
                        out.extend_from_slice(value);
                        continue;
                    }
                    match site_dp.get(column - 9).filter(|depth| !depth.is_missing()) {
                        Some(depth) => out.extend_from_slice(depth.to_string().as_bytes()),
                        None => out.push(b'.'),
                    }
                }
            }
            _ => out.extend_from_slice(field),
        }
    }
    out.push(b'\n');
}

/// The sort key of a record on the contig with this index, which is also what makes two records
/// the same call.
fn sort_key(rid: u32, call: &Normalized) -> SortKey<'_> {
    let alt_allele = if call.ref_allele == call.alt_allele {
        &[][..]
    } else {
        call.alt_allele.as_bytes()
    };
    SortKey {
        rid,
        pos: call.start - 1,
        ref_allele: call.ref_allele.as_bytes(),
        alt_allele,
        end: call.end,
    }
}

/// Whether a normalized call substitutes bases, as only a substitution can be a leaky MNV's partial.
fn is_substitution(call: &Normalized) -> bool {
    call.ref_allele.len() == call.alt_allele.len()
        && call.ref_allele != call.alt_allele
        && call
            .alt_allele
            .bytes()
            .all(|base| base.is_ascii_alphabetic())
}

/// Copies a row into a record holding only its fields, without `-fisher`'s, and its position, as a
/// region's rows are held until the region ends.
fn compact(record: &StringRecord, columns: Columns) -> StringRecord {
    let mut copy = StringRecord::with_capacity(record.as_slice().len(), record.len());
    columns.copy_row(record, &mut copy);
    copy
}

/// What a record carries besides its calls: the merge of an MNV that absorbed partials, the reads
/// a partial lost to its MNV per sample, its reference 3-mer, its depth at POS per sample, FORMAT
/// SITE_DP, left out when empty, and whether its genotype is missing whatever its reads.
#[derive(Clone, Copy, Default)]
struct Extras<'e> {
    merged: Option<&'e Merged<'e>>,
    moved: Option<&'e [i32]>,
    context: Option<&'e str>,
    site_dp: &'e [i32],
    unknown_genotype: bool,
}

/// Builds one VCF record with a sample column per call, the tumor (or only) sample first; FILTER and
/// the site's INFO come from the first call, and a tumor-normal record carries its VARDICT_STATUS.
/// The record also carries its [`Extras`]. Returns `None` for a record the writer skips.
fn build_record<'w>(
    records: &'w mut RecordWriter,
    calls: &[TumorOnlyVariant],
    status: Option<&str>,
    filters: &FilterThresholds,
    extras: Extras,
) -> Result<Option<Record<'w>>, Box<dyn error::Error>> {
    let var = &calls[0];
    let Some(mut variant) = records.record(
        var.contig,
        var.start,
        var.ref_allele,
        &var.alt_allele_for_vcf(),
    )?
    else {
        return Ok(None);
    };

    if filters.any() && var.has_alt_reads() {
        let labels = filters.labels(var);
        if labels.is_empty() {
            variant.push_filter("PASS")?;
        }
        for label in labels {
            variant.push_filter(label)?;
        }
    }

    if let Some(class) = var.variant_class() {
        variant.push_info_string("TYPE", class)?;
    }

    if let Some(shift) = var.indel_3p_shift_value() {
        variant.push_info_integer("INDEL_3P_SHIFT", &[shift])?;
    }

    if let (Some(copies), Some(unit_length)) =
        (var.repeat_unit_copies_value(), var.repeat_unit_len_value())
    {
        variant.push_info_float("REPEAT_UNIT_COPIES", &[copies])?;
        variant.push_info_integer("REPEAT_UNIT_LEN", &[unit_length])?;
    }

    if VALID_SV_TYPES.contains(&var.variant_type) {
        if var.end < var.start {
            warn!(
                "The {} at {}:{} ends at {}, before its POS; its INFO END is written as VarDict gives it.",
                var.variant_type, var.contig, var.start, var.end
            );
        }
        variant.push_info_integer("END", &[var.end as i32])?;
        let sv_call = calls
            .iter()
            .find(|call| call.has_alt_reads())
            .unwrap_or(var);
        let sv_length = sv_call.sv_length().ok_or_else(|| {
            format!(
                "{}:{}: cannot read the {} length from the genotype column '{}'",
                var.contig, var.start, var.variant_type, sv_call.gt
            )
        })?;
        variant.push_info_integer("SVLEN", &[sv_length])?;
        variant.push_info_string("SVTYPE", var.variant_type)?;
    }

    if let (Some(status), [tumor, normal]) = (status, calls) {
        variant.push_info_string("VARDICT_STATUS", status)?;
        variant.push_info_float(
            "TUMOR_NORMAL_FISHER_P",
            &[tumor_normal_fisher_p_value(tumor, normal)],
        )?;
    }

    if let Some(merged) = extras.merged {
        variant.push_info_integer("MNV_PARTIALS", &[merged.partials as i32])?;
    }

    if let Some(context) = extras.context {
        variant.push_info_string("CONTEXT_3mer", context)?;
    }

    let genotypes = Values::<_, 4>::collect(
        GenotypeAllele::UnphasedMissing,
        calls.iter().flat_map(|call| {
            let missing = &[GenotypeAllele::UnphasedMissing; 2][..];
            let alleles = if extras.unknown_genotype {
                missing
            } else {
                call.gt_value(MIN_HOM_ALT_AF)
            };
            alleles.iter().copied()
        }),
    );
    variant.push_genotypes(&genotypes)?;
    variant.push_format_integer("AD", &each_allele(calls, TumorOnlyVariant::ad_value))?;
    variant.push_format_integer("ADF", &each_allele(calls, TumorOnlyVariant::adf_value))?;
    variant.push_format_integer("ADR", &each_allele(calls, TumorOnlyVariant::adr_value))?;
    variant.push_format_float(
        "STRAND_BIAS_FISHER_P",
        &each(calls, TumorOnlyVariant::strand_bias_fisher_p_value),
    )?;
    variant.push_format_integer("DP", &each(calls, |call| call.depth))?;
    variant.push_format_float("AF", &each(calls, TumorOnlyVariant::af_value))?;
    variant.push_format_integer("HICNT", &each(calls, TumorOnlyVariant::hicnt_value))?;
    variant.push_format_float(
        "REALIGNED_FRAC_OF_DP",
        &each(calls, TumorOnlyVariant::realigned_frac_of_dp_value),
    )?;
    variant.push_format_float(
        "MEAN_DIST_TO_READ_END",
        &each(calls, TumorOnlyVariant::mean_dist_to_read_end_value),
    )?;
    variant.push_format_integer(
        "ALT_READ_POS_VARIES",
        &each(calls, TumorOnlyVariant::alt_read_pos_varies_value),
    )?;
    variant.push_format_float("QMEAN", &each(calls, TumorOnlyVariant::qmean_value))?;
    variant.push_format_float("MEAN_MAPQ", &each(calls, TumorOnlyVariant::mean_mapq_value))?;
    variant.push_format_float(
        "MEAN_MISMATCHES",
        &each(calls, TumorOnlyVariant::mean_mismatches_value),
    )?;

    if let Some(merged) = extras.merged {
        variant.push_format_integer("MNV_PARTIAL_READS", &merged.partial_reads)?;
    }
    if let Some(moved) = extras.moved {
        variant.push_format_integer("MNV_MOVED_READS", moved)?;
    }

    if let Some(counts) = calls
        .iter()
        .map(TumorOnlyVariant::sv_read_counts)
        .collect::<Option<Vec<_>>>()
    {
        variant.push_format_integer("SV_SOFTCLIP_READS", &each(&counts, |count| count.0))?;
        variant.push_format_integer("SV_DISCORDANT_READS", &each(&counts, |count| count.1))?;
    }
    variant.push_format_integer("SITE_DP", extras.site_dp)?;

    Ok(Some(variant))
}

/// Collects one value per sample.
fn each<T, V: Copy + Default>(items: &[T], value: impl Fn(&T) -> V) -> Values<V, 2> {
    Values::collect(V::default(), items.iter().map(value))
}

/// Collects one value per allele for each sample, sample after sample.
fn each_allele<T>(items: &[T], value: impl Fn(&T) -> Vec<i32>) -> Values<i32, 4> {
    Values::collect(0, items.iter().flat_map(value))
}

/// Up to `N` values held on the stack, sparing a vector for each FORMAT field of each record.
struct Values<V, const N: usize> {
    values: [V; N],
    len: usize,
}

impl<V: Copy, const N: usize> Values<V, N> {
    /// Collects the values after a placeholder fills the unused slots.
    fn collect(fill: V, values: impl IntoIterator<Item = V>) -> Self {
        let mut collected = Values {
            values: [fill; N],
            len: 0,
        };
        for value in values {
            collected.values[collected.len] = value;
            collected.len += 1;
        }
        collected
    }
}

impl<V, const N: usize> Deref for Values<V, N> {
    type Target = [V];

    fn deref(&self) -> &[V] {
        &self.values[..self.len]
    }
}

/// The substitutions added so far, each held as a fixed-size key, for skipping one added before.
#[derive(Default)]
struct SeenRows {
    keys: AHashSet<RowKey>,
    strings: AHashMap<Box<str>, u32>,
}

/// A call's contig, start, END, REF and ALT, with each string held as an id.
#[derive(Hash, PartialEq, Eq)]
struct RowKey {
    contig: u32,
    start: u32,
    end: u32,
    ref_allele: u32,
    alt_allele: u32,
}

impl SeenRows {
    /// Records a normalized call, returning whether it is new; a position needs to fit in 32 bits.
    fn insert(&mut self, contig: &str, call: &Normalized) -> Result<bool, String> {
        let (Ok(start), Ok(end)) = (call.start.try_into(), call.end.try_into()) else {
            return Err(format!(
                "Expected a start and end below 4294967296, found {contig}:{}-{}!",
                call.start, call.end
            ));
        };
        let key = RowKey {
            start,
            end,
            contig: self.id(contig),
            ref_allele: self.id(&call.ref_allele),
            alt_allele: self.id(&call.alt_allele),
        };
        Ok(self.keys.insert(key))
    }

    /// Whether a normalized call was recorded.
    fn contains(&self, contig: &str, call: &Normalized) -> bool {
        let (Ok(start), Ok(end)) = (call.start.try_into(), call.end.try_into()) else {
            return false;
        };
        let id = |string: &str| self.strings.get(string).copied();
        match (id(contig), id(&call.ref_allele), id(&call.alt_allele)) {
            (Some(contig), Some(ref_allele), Some(alt_allele)) => self.keys.contains(&RowKey {
                contig,
                start,
                end,
                ref_allele,
                alt_allele,
            }),
            _ => false,
        }
    }

    /// Returns the id of a string, giving it the next id the first time it is seen.
    fn id(&mut self, string: &str) -> u32 {
        if let Some(&id) = self.strings.get(string) {
            return id;
        }
        let id = u32::try_from(self.strings.len()).expect("Over 2^32 distinct strings!");
        self.strings.insert(string.into(), id);
        id
    }
}

/// Reads the next row, refusing one whose columns differ from the first row's.
fn read_row<R: Read>(
    reader: &mut Reader<R>,
    record: &mut StringRecord,
    columns: Columns,
) -> Result<bool, Box<dyn error::Error>> {
    if !reader.read_record(record)? {
        return Ok(false);
    }
    if !is_segment(record, columns.segment_column()) {
        return Err(layout_error(record, Some(columns)).into());
    }
    check_width(record, columns)?;
    Ok(true)
}

/// Refuses a row with fewer columns than its layout has.
fn check_width(record: &StringRecord, columns: Columns) -> Result<(), String> {
    if record.len() >= columns.min_columns() {
        return Ok(());
    }
    let line = record.position().map_or(0, |position| position.line());
    Err(format!(
        "Expected at least {} columns on line {line}, found {}!",
        columns.min_columns(),
        record.len()
    ))
}

/// Detects the columns of a row from the column holding its segment.
fn detect_columns(record: &StringRecord) -> Result<Columns, String> {
    Columns::ALL
        .into_iter()
        .find(|columns| is_segment(record, columns.segment_column()))
        .ok_or_else(|| layout_error(record, None))
}

/// Whether a row holds a `chr:start-end` segment on its own contig at a 0-based column.
fn is_segment(record: &StringRecord, index: usize) -> bool {
    let contig = record.get(2).unwrap_or_default();
    record
        .get(index)
        .and_then(|field| {
            field
                .strip_prefix(contig)?
                .strip_prefix(':')?
                .split_once('-')
        })
        .is_some_and(|(start, end)| start.parse::<u64>().is_ok() && end.parse::<u64>().is_ok())
}

/// Describes why a row matches no layout, or not the first row's columns.
fn layout_error(record: &StringRecord, first: Option<Columns>) -> String {
    let line = record.position().map_or(0, |position| position.line());
    let contig = record.get(2).unwrap_or_default();
    match first {
        Some(columns) => format!(
            "Expected a {contig}:start-end segment in column {} on line {line}, as on the first row!",
            columns.segment_column() + 1
        ),
        None => format!(
            "Expected a {contig}:start-end segment in column 33 or 35 (tumor-only) or 49 or 53 (tumor-normal) on line {line}!"
        ),
    }
}

/// Resolves the tumor and normal sample names from the first row's sample, the layout and the
/// arguments. A tumor-normal row names both samples as `tumor|normal` when VarDict was given both
/// names with a BED file of regions, and only the tumor otherwise. Given names, also accepted as
/// `tumor|normal`, take the place of a tumor-normal input's and must match a tumor-only input's.
fn resolve_samples(
    row_sample: Option<&str>,
    layout: Option<Layout>,
    sample: Option<&str>,
    normal_sample: Option<&str>,
) -> Result<(String, Option<String>), String> {
    let (sample, normal_sample) = match (sample.and_then(|name| name.split_once('|')), layout) {
        (Some(_), Some(Layout::TumorNormal) | None) if normal_sample.is_some() => {
            return Err("Give the normal sample's name in --sample \"tumor|normal\" or in --normal-sample, not both!".to_string());
        }
        (Some((tumor, normal)), Some(Layout::TumorNormal) | None) => (Some(tumor), Some(normal)),
        _ => (sample, normal_sample),
    };
    let (row_tumor, row_normal) = match (layout, row_sample.map(|row| row.split_once('|'))) {
        (Some(Layout::TumorNormal), Some(Some((tumor, normal)))) => (Some(tumor), Some(normal)),
        _ => (row_sample, None),
    };
    if let (Some(Layout::TumorOnly), Some(given), Some(row)) = (layout, sample, row_tumor)
        && given != row
    {
        return Err(format!("Expected sample '{given}' found '{row}'!"));
    }
    let tumor = given_or_input(sample, row_tumor).ok_or("The input has no rows to take a sample name from: give --sample, and --normal-sample for tumor-normal output!")?;
    let normal = match layout {
        Some(Layout::TumorOnly) if normal_sample.is_some() => {
            return Err("--normal-sample was given but the input is tumor-only!".to_string());
        }
        Some(Layout::TumorNormal) => Some(
            given_or_input(normal_sample, row_normal)
                .ok_or("The input is tumor-normal but no --normal-sample was given!")?,
        ),
        _ => normal_sample,
    };
    for name in [Some(tumor), normal].into_iter().flatten() {
        if name.trim().is_empty() {
            return Err("Sample names must not be empty or blank!".to_string());
        }
        if name.contains(['\t', '\n', '\r']) {
            return Err("Sample names must not contain tabs or line breaks!".to_string());
        }
    }
    if normal == Some(tumor) {
        return Err(format!(
            "The tumor and normal samples need different names, but both are '{tumor}'!"
        ));
    }
    Ok((tumor.to_string(), normal.map(str::to_string)))
}

/// Returns the given sample name, warning when the input names the sample differently, or else
/// the input's name.
fn given_or_input<'a>(given: Option<&'a str>, input: Option<&'a str>) -> Option<&'a str> {
    if let (Some(given), Some(input)) = (given, input)
        && given != input
    {
        warn!("The input names a sample '{input}'; writing it as '{given}' as given.");
    }
    given.or(input)
}

/// Names the 1-based line and input column of a record that could not be deserialized.
fn describe_parse_error(error: csv::Error, columns: Columns) -> String {
    match error.kind() {
        csv::ErrorKind::Deserialize {
            pos: Some(pos),
            err,
        } => match err.field() {
            Some(index) => format!(
                "Could not parse column {} on line {}: {}!",
                columns.input_column(index as usize) + 1,
                pos.line(),
                err.kind()
            ),
            None => format!("Could not parse line {}: {}!", pos.line(), err.kind()),
        },
        _ => error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use anyhow::Result;
    use file_diff::diff;
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_vartovcf_run() -> Result<(), Box<dyn std::error::Error>> {
        let sample = "dna00001";
        let input = BufReader::new(File::open("tests/calls.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(
            input,
            &Options {
                output: Some(output.path().into()),
                reference: reference.clone(),
                sample: Some(sample),
                ..Default::default()
            },
        )?;
        assert_eq!(exit, 0);
        assert!(diff(output.path().to_str().unwrap(), "tests/calls.vcf"));
        Ok(())
    }

    #[test]
    fn test_writing_stops_once_a_signal_arrives() {
        let rows = std::fs::read_to_string("tests/calls.g.var").unwrap();
        let calls = std::cell::Cell::new(0);
        let signal = || {
            calls.set(calls.get() + 1);
            if calls.get() > rows.lines().count() {
                2
            } else {
                0
            }
        };
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("out.vcf");
        let mut output = Output::create(Some(&path)).unwrap();
        let options = Options {
            reference: PathBuf::from("tests/reference.fa"),
            sample: Some("dna00001"),
            ..Default::default()
        };
        let index =
            crate::convert(rows.as_bytes(), &options, dir.path(), &mut output, &signal).unwrap();
        assert!(index.is_none());
        assert_eq!(calls.get(), rows.lines().count() + 1);
    }

    #[test]
    fn test_filling_stops_inside_a_gap_once_a_signal_arrives() {
        let rows = std::fs::read_to_string("tests/calls.contig-ends.var").unwrap();
        let stop_at = rows.lines().count() + 20;
        let calls = std::cell::Cell::new(0);
        let signal = || {
            calls.set(calls.get() + 1);
            if calls.get() > stop_at { 2 } else { 0 }
        };
        let dir = tempfile::TempDir::new().unwrap();
        let mut output = Output::create(Some(&dir.path().join("out.vcf"))).unwrap();
        let options = Options {
            reference: PathBuf::from("tests/reference.fa"),
            sample: Some("dna00001"),
            fill_gvcf: true,
            ..Default::default()
        };
        let index =
            crate::convert(rows.as_bytes(), &options, dir.path(), &mut output, &signal).unwrap();
        assert!(index.is_none());
        assert_eq!(calls.get(), stop_at + 1);
    }

    #[test]
    fn test_stamp_writes_each_samples_site_dp_in_place() {
        let mut out = Vec::new();
        stamp(
            b"c\t5\t.\tA\tC\t.\t.\t.\tGT:DP:SITE_DP:AF\t0/1:30:.:0.5\t0/0:20:.:0\n",
            &[25, i32::missing()],
            &mut out,
        );
        assert_eq!(
            out,
            b"c\t5\t.\tA\tC\t.\t.\t.\tGT:DP:SITE_DP:AF\t0/1:30:25:0.5\t0/0:20:.:0\n"
        );
        stamp(b"c\t5\t.\tA\tC\t.\t.\t.\n", &[25], &mut out);
        assert_eq!(out, b"c\t5\t.\tA\tC\t.\t.\t.\n");
    }

    #[test]
    fn test_vartovcf_refuses_a_fill_bed_without_fill_gvcf() {
        let input = BufReader::new(File::open("tests/calls.fill-gvcf.var").unwrap());
        let error = vartovcf(
            input,
            &Options {
                output: Some("/dev/null".into()),
                reference: PathBuf::from("tests/reference.fa"),
                fill_gvcf_bed: Some(PathBuf::from("tests/calls.fill-gvcf.bed")),
                ..Default::default()
            },
        )
        .unwrap_err();
        assert_eq!(error.to_string(), "--fill-gvcf-bed requires --fill-gvcf!");
    }

    #[test]
    fn test_vartovcf_refuses_a_fill_that_skips_non_variant_rows() {
        let input = BufReader::new(File::open("tests/calls.fill-gvcf.var").unwrap());
        let error = vartovcf(
            input,
            &Options {
                output: Some("/dev/null".into()),
                reference: PathBuf::from("tests/reference.fa"),
                fill_gvcf: true,
                skip_non_variants: true,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert_eq!(
            error.to_string(),
            "--fill-gvcf cannot be used with --skip-non-variants!"
        );
    }

    #[test]
    fn test_a_bare_zero_strand_bias_column_is_accepted() -> Result<(), Box<dyn std::error::Error>> {
        let rows = std::fs::read_to_string("tests/calls.var")?;
        let mut fields: Vec<&str> = rows.lines().nth(2).unwrap().split('\t').collect();
        fields[15] = "0";
        let input = format!("{}\n", fields.join("\t"));
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let exit = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/reference.fa"),
                sample: Some("dna00001"),
                ..Default::default()
            },
        )?;
        assert_eq!(exit, 0);
        Ok(())
    }

    #[rstest]
    #[case(
        "tests/calls.var",
        36,
        "Expected at least 38 columns on line 1, found 36!"
    )]
    #[case(
        "tests/calls.no-fisher.var",
        34,
        "Expected at least 36 columns on line 1, found 34!"
    )]
    #[case(
        "tests/calls.tumor-normal.var",
        57,
        "Expected at least 61 columns on line 1, found 57!"
    )]
    #[case(
        "tests/calls.tumor-normal.no-fisher.var",
        50,
        "Expected at least 55 columns on line 1, found 50!"
    )]
    fn test_refuses_a_truncated_row(
        #[case] path: &str,
        #[case] columns: usize,
        #[case] message: &str,
    ) {
        let rows = std::fs::read_to_string(path).unwrap();
        let fields: Vec<&str> = rows.lines().next().unwrap().split('\t').collect();
        let input = format!("{}\n", fields[..columns].join("\t"));
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let result = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/reference.fa"),
                normal_sample: path.contains("tumor-normal").then_some("N"),
                ..Default::default()
            },
        );
        assert_eq!(result.unwrap_err().to_string(), message);
    }

    #[test]
    fn test_refuses_a_truncated_row_after_the_first() {
        let rows = std::fs::read_to_string("tests/calls.var").unwrap();
        let mut lines = rows.lines();
        let first = lines.next().unwrap();
        let second: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let input = format!("{first}\n{}\n", second[..36].join("\t"));
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let result = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/reference.fa"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Expected at least 38 columns on line 2, found 36!"
        );
    }

    #[test]
    fn test_an_unreadable_sv_length_names_the_genotype_column_it_read() {
        let rows = std::fs::read_to_string("tests/calls.tumor-normal.edge.var").unwrap();
        let input = rows.replace("A/-199", "A/A");
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let result = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/tumor-normal.fa"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "cK:160: cannot read the DEL length from the genotype column 'A/A'"
        );
    }

    #[test]
    fn test_refuses_a_mismatched_tumor_only_sample() {
        let input = BufReader::new(File::open("tests/calls.var").unwrap());
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            &Options {
                output: Some(output.path().into()),
                reference: reference.clone(),
                sample: Some("XXXXXXXX"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Expected sample 'XXXXXXXX' found 'dna00001'!"
        );
    }

    /// Drops the columns VarDictJava's `-fisher` adds from rows of either layout.
    fn without_fisher_columns(rows: &str) -> String {
        rows.lines()
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                let added: &[usize] = match fields.len() {
                    38 | 39 => &[20, 21],
                    61 => &[25, 26, 45, 46, 59, 60],
                    width => panic!("Unexpected row of {width} columns!"),
                };
                let kept: Vec<&str> = fields
                    .iter()
                    .enumerate()
                    .filter(|(column, _)| !added.contains(column))
                    .map(|(_, field)| *field)
                    .collect();
                format!("{}\n", kept.join("\t"))
            })
            .collect()
    }

    fn convert(input: &str, options: &Options) -> Result<Vec<u8>, Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let options = Options {
            output: Some(output.path().into()),
            ..options.clone()
        };
        assert_eq!(vartovcf(input.as_bytes(), &options)?, 0);
        Ok(std::fs::read(output.path())?)
    }

    #[rstest]
    #[case("tests/calls.var", "tests/reference.fa", false)]
    #[case("tests/calls.filters.var", "tests/reference.fa", false)]
    #[case("tests/calls.g.var", "tests/reference.fa", false)]
    #[case("tests/calls.complex.var", "tests/reference.fa", false)]
    #[case("tests/calls.normalize.var", "tests/reference.fa", false)]
    #[case("tests/calls.leaky-mnv.var", "tests/reference.fa", true)]
    #[case("tests/calls.leaky-mnv.overlapping.var", "tests/reference.fa", true)]
    #[case("tests/calls.tumor-normal.var", "tests/tumor-normal.fa", false)]
    #[case("tests/calls.tumor-normal.edge.var", "tests/tumor-normal.fa", false)]
    #[case("tests/calls.tumor-normal.hicnt.var", "tests/tumor-normal.fa", false)]
    #[case(
        "tests/calls.tumor-normal.normalize.var",
        "tests/tumor-normal.fa",
        false
    )]
    #[case(
        "tests/calls.tumor-normal.leaky-mnv.var",
        "tests/tumor-normal.fa",
        true
    )]
    fn test_rows_without_fisher_columns_convert_the_same(
        #[case] path: &str,
        #[case] reference: &str,
        #[case] fix_leaky_mnv: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let rows = std::fs::read_to_string(path)?;
        let options = Options {
            reference: PathBuf::from(reference),
            normal_sample: path.contains("tumor-normal").then_some("N"),
            fix_leaky_mnv: fix_leaky_mnv.then_some(LeakyMnvFix::Counts(LeakRate::default())),
            filters: FilterThresholds {
                near_read_end: Some(8.0),
                low_mean_mapq: Some(10.0),
                high_mean_mismatches: Some(5.25),
                same_read_position: Some(0.35),
                strand_bias: Some(0.01),
                strand_bias_min_odds_ratio: Some(5.0),
                strand_bias_max_af: Some(0.25),
                low_af: Some(0.0001),
                low_qmean: Some(30.0),
                low_dp: Some(3),
                ..Default::default()
            },
            ..Default::default()
        };
        let expected = convert(&rows, &options)?;
        assert_eq!(convert(&without_fisher_columns(&rows), &options)?, expected);
        Ok(())
    }

    #[rstest]
    #[case("tests/calls.var")]
    #[case("tests/calls.tumor-normal.var")]
    #[case("tests/calls.tumor-normal.edge.var")]
    fn test_strand_bias_p_values_round_to_the_ones_vardictjava_prints(#[case] path: &str) {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(path)
            .unwrap();
        let mut checked = 0;
        for row in reader.records() {
            let row = row.unwrap();
            let columns = detect_columns(&row).unwrap();
            let mut record = row.clone();
            columns.drop_fisher_columns(&mut record, &mut StringRecord::new());
            let (calls, _) = parse_row(&record, columns).unwrap();
            let printed = columns.fisher_columns().iter().step_by(2);
            for (call, &column) in calls.iter().zip(printed).take(columns.layout.samples()) {
                if call.has_alt_reads() {
                    let printed: f64 = row[column].parse().unwrap();
                    let p = call.strand_bias_fisher_p_value() as f64;
                    assert!((p - printed).abs() <= 5e-6, "{p} printed as {printed}");
                    checked += 1;
                }
            }
        }
        assert!(checked > 0);
    }

    #[test]
    fn test_refuses_a_row_without_the_first_rows_fisher_columns() {
        let with_fisher = std::fs::read_to_string("tests/calls.var").unwrap();
        let without_fisher = std::fs::read_to_string("tests/calls.no-fisher.var").unwrap();
        let input = format!(
            "{}\n{}\n",
            with_fisher.lines().nth(2).unwrap(),
            without_fisher.lines().nth(2).unwrap()
        );
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let result = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/reference.fa"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Expected a chr1:start-end segment in column 35 on line 2, as on the first row!"
        );
    }

    #[rstest]
    #[case(
        "tests/calls.var",
        2,
        22,
        "Could not parse column 23 on line 1: invalid float literal!"
    )]
    #[case(
        "tests/calls.no-fisher.var",
        2,
        20,
        "Could not parse column 21 on line 1: invalid float literal!"
    )]
    #[case(
        "tests/calls.tumor-normal.var",
        0,
        47,
        "Could not parse column 48 on line 1: invalid digit found in string!"
    )]
    #[case(
        "tests/calls.tumor-normal.no-fisher.var",
        0,
        43,
        "Could not parse column 44 on line 1: invalid digit found in string!"
    )]
    fn test_a_parse_error_names_the_column_as_read(
        #[case] path: &str,
        #[case] row: usize,
        #[case] column: usize,
        #[case] message: &str,
    ) {
        let rows = std::fs::read_to_string(path).unwrap();
        let mut fields: Vec<&str> = rows.lines().nth(row).unwrap().split('\t').collect();
        fields[column] = "x";
        let input = format!("{}\n", fields.join("\t"));
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let result = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from(if path.contains("tumor-normal") {
                    "tests/tumor-normal.fa"
                } else {
                    "tests/reference.fa"
                }),
                normal_sample: path.contains("tumor-normal").then_some("N"),
                ..Default::default()
            },
        );
        assert_eq!(result.unwrap_err().to_string(), message);
    }

    #[rstest]
    #[case("tests/calls.var", Columns::new(Layout::TumorOnly, true))]
    #[case("tests/calls.no-fisher.var", Columns::new(Layout::TumorOnly, false))]
    fn test_reads_the_crispr_column_with_or_without_fisher(
        #[case] path: &str,
        #[case] columns: Columns,
    ) {
        let rows = std::fs::read_to_string(path).unwrap();
        let mut record = first_record(&format!("{}\t-3\n", rows.lines().nth(2).unwrap()));
        assert_eq!(detect_columns(&record), Ok(columns));
        columns.drop_fisher_columns(&mut record, &mut StringRecord::new());
        let (calls, _) = parse_row(&record, columns).unwrap();
        assert_eq!(calls[0].distance_to_crispr_site, Some(-3));
        assert_eq!(calls[0].mean_mapping_quality, 60.0);
    }

    #[test]
    fn test_malformed_row_is_an_error() -> Result<(), Box<dyn std::error::Error>> {
        let input = BufReader::new(File::open("tests/calls.malformed.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            &Options {
                output: Some(output.path().into()),
                reference: reference.clone(),
                sample: Some("dna00001"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Could not parse column 8 on line 2: invalid digit found in string!"
        );
        Ok(())
    }

    #[test]
    fn test_refuses_rows_without_a_segment() {
        let input = "dna00001\tBRINP3\tchr1\t190098265\t190098265\tA\tA\n".as_bytes();
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let result = vartovcf(
            input,
            &Options {
                output: Some(output.path().into()),
                reference: reference.clone(),
                sample: Some("dna00001"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Expected a chr1:start-end segment in column 33 or 35 (tumor-only) or 49 or 53 (tumor-normal) on line 1!"
        );
    }

    fn first_record(input: &str) -> StringRecord {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(input.as_bytes());
        let mut record = StringRecord::new();
        assert!(reader.read_record(&mut record).unwrap());
        record
    }

    #[rstest]
    #[case("tests/calls.var", Layout::TumorOnly, true, 38)]
    #[case("tests/calls.no-fisher.var", Layout::TumorOnly, false, 36)]
    #[case("tests/calls.tumor-normal.var", Layout::TumorNormal, true, 61)]
    #[case(
        "tests/calls.tumor-normal.no-fisher.var",
        Layout::TumorNormal,
        false,
        55
    )]
    fn test_detect_columns(
        #[case] path: &str,
        #[case] layout: Layout,
        #[case] fisher: bool,
        #[case] width: usize,
    ) {
        let mut record = first_record(&std::fs::read_to_string(path).unwrap());
        let columns = detect_columns(&record);
        assert_eq!(columns, Ok(Columns::new(layout, fisher)));
        assert_eq!(record.len(), width);
        columns
            .unwrap()
            .drop_fisher_columns(&mut record, &mut StringRecord::new());
        assert_eq!(record.len(), layout.min_columns());
        assert!(is_segment(&record, layout.segment_column()));
    }

    #[test]
    fn test_detect_columns_refuses_rows_without_a_segment() {
        let input = "dna00001\tBRINP3\tchr1\t190098265\t190098265\tA\tA\n";
        assert_eq!(
            detect_columns(&first_record(input)),
            Err("Expected a chr1:start-end segment in column 33 or 35 (tumor-only) or 49 or 53 (tumor-normal) on line 1!".to_string())
        );
    }

    #[rstest]
    #[case(Some("dna00001"), Some(Layout::TumorOnly), None, None, Ok(("dna00001", None)))]
    #[case(Some("dna00001"), Some(Layout::TumorOnly), Some("dna00001"), None, Ok(("dna00001", None)))]
    #[case(
        Some("dna00001"),
        Some(Layout::TumorOnly),
        Some("XXXX"),
        None,
        Err("Expected sample 'XXXX' found 'dna00001'!")
    )]
    #[case(
        Some("dna00001"),
        Some(Layout::TumorOnly),
        Some("a|b"),
        None,
        Err("Expected sample 'a|b' found 'dna00001'!")
    )]
    #[case(
        Some("dna00001"),
        Some(Layout::TumorOnly),
        None,
        Some("normal"),
        Err("--normal-sample was given but the input is tumor-only!")
    )]
    #[case(
        Some("tumor"),
        Some(Layout::TumorNormal),
        None,
        None,
        Err("The input is tumor-normal but no --normal-sample was given!")
    )]
    #[case(Some("tumor"), Some(Layout::TumorNormal), None, Some("normal"), Ok(("tumor", Some("normal"))))]
    #[case(
        None,
        None,
        None,
        None,
        Err(
            "The input has no rows to take a sample name from: give --sample, and --normal-sample for tumor-normal output!"
        )
    )]
    #[case(Some("T|N"), Some(Layout::TumorNormal), None, None, Ok(("T", Some("N"))))]
    #[case(Some("T|N"), Some(Layout::TumorNormal), Some("T"), Some("N"), Ok(("T", Some("N"))))]
    #[case(Some("T|N"), Some(Layout::TumorNormal), None, Some("X"), Ok(("T", Some("X"))))]
    #[case(Some("T|N"), Some(Layout::TumorNormal), Some("X"), None, Ok(("X", Some("N"))))]
    #[case(Some("T.bam|N"), Some(Layout::TumorNormal), Some("T"), Some("N"), Ok(("T", Some("N"))))]
    #[case(Some("T|N"), Some(Layout::TumorNormal), Some("A|B"), None, Ok(("A", Some("B"))))]
    #[case(Some("T"), Some(Layout::TumorNormal), Some("A|B"), None, Ok(("A", Some("B"))))]
    #[case(None, None, Some("A|B"), None, Ok(("A", Some("B"))))]
    #[case(
        None,
        None,
        Some(""),
        None,
        Err("Sample names must not be empty or blank!")
    )]
    #[case(
        None,
        None,
        Some(" "),
        None,
        Err("Sample names must not be empty or blank!")
    )]
    #[case(
        None,
        None,
        Some("T\tX"),
        None,
        Err("Sample names must not contain tabs or line breaks!")
    )]
    #[case(
        Some("T|N"),
        Some(Layout::TumorNormal),
        Some("T|N"),
        Some("N"),
        Err(
            "Give the normal sample's name in --sample \"tumor|normal\" or in --normal-sample, not both!"
        )
    )]
    #[case(
        Some("T"),
        Some(Layout::TumorNormal),
        None,
        Some(""),
        Err("Sample names must not be empty or blank!")
    )]
    #[case(
        Some("T|"),
        Some(Layout::TumorNormal),
        None,
        None,
        Err("Sample names must not be empty or blank!")
    )]
    #[case(
        Some("T"),
        Some(Layout::TumorNormal),
        None,
        Some("T"),
        Err("The tumor and normal samples need different names, but both are 'T'!")
    )]
    #[case(
        Some("T|T"),
        Some(Layout::TumorNormal),
        None,
        None,
        Err("The tumor and normal samples need different names, but both are 'T'!")
    )]
    #[case(None, None, Some("tumor"), None, Ok(("tumor", None)))]
    #[case(None, None, Some("tumor"), Some("normal"), Ok(("tumor", Some("normal"))))]
    fn test_resolve_samples(
        #[case] row_sample: Option<&str>,
        #[case] layout: Option<Layout>,
        #[case] sample: Option<&str>,
        #[case] normal_sample: Option<&str>,
        #[case] expected: Result<(&str, Option<&str>), &str>,
    ) {
        let expected = expected
            .map(|(tumor, normal)| (tumor.to_string(), normal.map(str::to_string)))
            .map_err(str::to_string);
        assert_eq!(
            resolve_samples(row_sample, layout, sample, normal_sample),
            expected
        );
    }

    #[rstest]
    #[case(None, vec!["tumor"])]
    #[case(Some("normal"), vec!["tumor", "normal"])]
    fn test_empty_input_writes_a_header_for_the_named_samples(
        #[case] normal_sample: Option<&str>,
        #[case] expected: Vec<&str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(
            "".as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: reference.clone(),
                sample: Some("tumor"),
                normal_sample,
                ..Default::default()
            },
        )?;
        assert_eq!(exit, 0);
        let mut reader = noodles::vcf::io::reader::Builder::default()
            .build_from_path(output.path())
            .expect("Error opening output file!");
        let samples: Vec<String> = reader
            .read_header()?
            .sample_names()
            .iter()
            .map(|sample| sample.to_string())
            .collect();
        assert_eq!(samples, expected);
        Ok(())
    }

    #[rstest]
    #[case(false, 11)]
    #[case(true, 7)]
    fn test_only_the_flag_folds_leaky_mnv_partials(
        #[case] fix_leaky_mnv: bool,
        #[case] expected: usize,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let input = BufReader::new(File::open("tests/calls.leaky-mnv.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        vartovcf(
            input,
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/reference.fa"),
                sample: Some("dna00001"),
                fix_leaky_mnv: fix_leaky_mnv.then_some(LeakyMnvFix::Counts(LeakRate::default())),
                ..Default::default()
            },
        )?;
        let mut reader = noodles::vcf::io::reader::Builder::default()
            .build_from_path(output.path())
            .expect("Error opening output file!");
        let header = reader.read_header()?;
        assert_eq!(header.infos().contains_key("MNV_PARTIALS"), fix_leaky_mnv);
        assert_eq!(
            header.formats().contains_key("MNV_PARTIAL_READS"),
            fix_leaky_mnv
        );
        assert_eq!(reader.records().count(), expected);
        Ok(())
    }

    fn row(
        contig: &'static str,
        start: u64,
        end: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> (&'static str, Normalized) {
        let call = Normalized {
            start,
            end,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
            indel_3p_shift: None,
        };
        (contig, call)
    }

    #[test]
    fn test_seen_rows_drops_the_same_call_twice() {
        let mut seen = SeenRows::default();
        let (contig, call) = row("chr1", 10, 12, "ACG", "A");
        assert_eq!(seen.insert(contig, &call), Ok(true));
        assert_eq!(seen.insert(contig, &call), Ok(false));
    }

    #[rstest]
    #[case::contig(row("chr2", 10, 12, "ACG", "A"))]
    #[case::start(row("chr1", 11, 12, "ACG", "A"))]
    #[case::end(row("chr1", 10, 13, "ACG", "A"))]
    #[case::ref_allele(row("chr1", 10, 12, "ACT", "A"))]
    #[case::alt_allele(row("chr1", 10, 12, "ACG", "<DEL>"))]
    #[case::swapped_alleles(row("chr1", 10, 12, "A", "ACG"))]
    fn test_seen_rows_keeps_a_call_differing_in_one_field(
        #[case] other: (&'static str, Normalized),
    ) {
        let mut seen = SeenRows::default();
        let (contig, call) = row("chr1", 10, 12, "ACG", "A");
        assert_eq!(seen.insert(contig, &call), Ok(true));
        assert_eq!(seen.insert(other.0, &other.1), Ok(true));
        assert_eq!(seen.insert(other.0, &other.1), Ok(false));
    }

    #[test]
    fn test_seen_rows_holds_each_distinct_string_once() {
        let mut seen = SeenRows::default();
        for start in 1..=100 {
            for alt_allele in ["A", "T"] {
                let (contig, call) = row("chr1", start, start, "A", alt_allele);
                seen.insert(contig, &call).unwrap();
            }
        }
        assert_eq!(seen.keys.len(), 200);
        assert_eq!(seen.strings.len(), 3);
    }

    #[rstest]
    #[case(1 << 32, 1 << 32)]
    #[case(1, 1 << 32)]
    fn test_seen_rows_refuses_a_position_beyond_32_bits(#[case] start: u64, #[case] end: u64) {
        let mut seen = SeenRows::default();
        let (contig, call) = row("chr1", start, end, "A", "T");
        assert_eq!(
            seen.insert(contig, &call),
            Err(format!(
                "Expected a start and end below 4294967296, found chr1:{start}-{end}!"
            ))
        );
    }

    #[test]
    fn test_refuses_a_position_beyond_the_reference() {
        let rows = std::fs::read_to_string("tests/calls.var").unwrap();
        let mut fields: Vec<&str> = rows.lines().nth(2).unwrap().split('\t').collect();
        fields[3] = "4294967296";
        fields[4] = "4294967296";
        let input = format!("{}\n", fields.join("\t"));
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let result = vartovcf(
            input.as_bytes(),
            &Options {
                output: Some(output.path().into()),
                reference: PathBuf::from("tests/reference.fa"),
                ..Default::default()
            },
        );
        assert_eq!(
            result.unwrap_err().to_string(),
            "Position chr1:4294967296 is beyond the end of chr1, which has 4000 bases in the reference: give --reference the FASTA VarDict was run with!"
        );
    }

    #[test]
    fn test_skip_non_variants_flag() -> Result<(), Box<dyn std::error::Error>> {
        let sample = "dna00001";
        let input = BufReader::new(File::open("tests/calls.g.var")?);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let reference = PathBuf::from("tests/reference.fa");
        let exit = vartovcf(
            input,
            &Options {
                output: Some(output.path().into()),
                reference: reference.clone(),
                sample: Some(sample),
                skip_non_variants: true,
                ..Default::default()
            },
        )?;
        assert_eq!(exit, 0);

        // Read the output and verify no non-variant sites exist
        use noodles::vcf::variant::record::AlternateBases;
        let mut reader = noodles::vcf::io::reader::Builder::default()
            .build_from_path(output.path())
            .expect("Error opening output file!");
        reader.read_header()?;

        let mut record_count = 0;
        for record_result in reader.records() {
            let record = record_result?;
            let alternate_bases = record.alternate_bases();
            let alleles: Vec<&[u8]> = std::iter::once(Ok(record.reference_bases()))
                .chain(alternate_bases.iter())
                .map(|allele| allele.map(str::as_bytes))
                .collect::<std::io::Result<_>>()?;

            // Verify that REF != ALT (i.e., no non-variant sites)
            // ALT should not be "." which represents non-variants
            assert_ne!(alleles.len(), 0, "Record should have alleles");
            if alleles.len() >= 2 {
                let alt = alleles[1];
                assert_ne!(
                    alt, b".",
                    "Non-variant site found when --skip-non-variants is true"
                );
            }
            record_count += 1;
        }

        // Verify we got some records (not all were filtered)
        assert!(record_count > 0, "No records found in output");
        Ok(())
    }
}
