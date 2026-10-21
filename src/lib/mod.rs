//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]

use ahash::{AHashMap, AHashSet};
use anyhow::Result;
use csv::{Reader, ReaderBuilder, StringRecord};
use log::warn;
use proglog::ProgLogBuilder;
use rust_htslib::bcf::Format;
use rust_htslib::bcf::Writer as VcfWriter;
use rust_htslib::bcf::record::{GenotypeAllele, Numeric};
use std::cmp::Reverse;
use std::error;
use std::fmt::Debug;
use std::io::Read;
use std::ops::Deref;
use std::path::PathBuf;

use crate::fai::{fasta_contigs_to_vcf_header, fasta_path_to_vcf_header};
use crate::filter::FilterThresholds;
use crate::io::has_gzip_ext;
use crate::mnv::{Fold, Merged, fold_partials};
use crate::normalize::{Normalized, normalize};
use crate::record::MIN_HOM_ALT_AF;
use crate::record::TumorNormalVariant;
use crate::record::TumorOnlyVariant;
use crate::record::{tumor_normal_fisher_p_value, tumor_normal_header, tumor_only_header};
use crate::reference::Reference;

pub mod fai;
pub mod filter;
pub mod fisher;
pub mod hicnt;
pub mod io;
pub mod mnv;
pub mod normalize;
pub mod record;
pub mod reference;

/// The structural variant types VarDict writes, which are the values of the `SVTYPE` INFO field.
pub const VALID_SV_TYPES: &[&str] = &["DEL", "DUP", "INV"];

/// Namespace for path parts and extensions.
pub mod path {

    /// The Gzip extension.
    pub const GZIP_EXTENSION: &str = "gz";
}

/// The layouts of VarDictJava output rows `vartovcf` reads, each from a run with `--fisher`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// One sample per row, from `vardict-java -b sample.bam`.
    TumorOnly,
    /// A tumor and its matched normal per row, from `vardict-java -b "tumor.bam|normal.bam"`.
    TumorNormal,
}

impl Layout {
    /// The 0-based column holding VarDict's `chr:start-end` segment.
    fn segment_column(&self) -> usize {
        match self {
            Layout::TumorOnly => 34,
            Layout::TumorNormal => 52,
        }
    }

    /// The number of samples in a row of this layout.
    fn samples(&self) -> usize {
        match self {
            Layout::TumorOnly => 1,
            Layout::TumorNormal => 2,
        }
    }

    /// The fewest columns a row of this layout has.
    fn min_columns(&self) -> usize {
        match self {
            Layout::TumorOnly => 38,
            Layout::TumorNormal => 61,
        }
    }
}

/// The 0-based segment columns of the tumor-only and tumor-normal layouts without `--fisher`.
const NO_FISHER_SEGMENT_COLUMNS: [usize; 2] = [32, 48];

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
    /// Fold the partial rows of each MNV back into it, see [`mnv`].
    pub fix_leaky_mnv: bool,
    /// The FILTER labels to apply.
    pub filters: FilterThresholds,
}

/// Runs the tool `vartovcf` on an input VAR file or stream and writes the records to a VCF file.
///
/// Each VarDict region's rows are parsed and, when asked, their leaky MNV partials folded, both in
/// VarDict's own coordinates; each record is then normalized against the reference and written
/// unless an earlier record normalized to the same call.
///
/// Returns the exit code for success (0).
pub fn vartovcf<I: Read>(input: I, options: &Options) -> Result<i32, Box<dyn error::Error>> {
    let reference = Reference::from_path(&options.reference)?;

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .flexible(true)
        .from_reader(input);

    // Peek at the first row to learn the layout and sample before the writer emits the header.
    let mut carry = StringRecord::new();
    let mut peeked = reader.read_record(&mut carry)?;
    let layout = if peeked {
        let layout = detect_layout(&carry)?;
        check_width(&carry, layout)?;
        Some(layout)
    } else {
        None
    };
    let row_sample = peeked.then(|| carry.get(0).unwrap_or_default().to_string());
    let (sample, normal_sample) = resolve_samples(
        row_sample.as_deref(),
        layout,
        options.sample,
        options.normal_sample,
    )?;
    let layout = layout.unwrap_or(Layout::TumorOnly);

    let filters = &options.filters;
    let mut header = match &normal_sample {
        Some(normal_sample) => {
            tumor_normal_header(&sample, normal_sample, filters, options.fix_leaky_mnv)
        }
        None => tumor_only_header(&sample, filters, options.fix_leaky_mnv),
    };

    fasta_contigs_to_vcf_header(&options.reference, &mut header);
    fasta_path_to_vcf_header(&options.reference, &mut header)
        .expect("Adding FASTA path to header failed!");

    let writer = match &options.output {
        Some(output) => VcfWriter::from_path(output, &header, !has_gzip_ext(output), Format::Vcf),
        None => VcfWriter::from_stdout(&header, true, Format::Vcf),
    }
    .expect("Could not build a VCF writer!");

    let mut sink = Sink {
        writer,
        reference,
        seen: SeenRows::default(),
        filters,
        progress: ProgLogBuilder::new()
            .name("main")
            .verb("Processed")
            .noun("variant records")
            .unit(100_000)
            .build(),
    };

    let row_sample = row_sample.unwrap_or_default();
    let mut region: Vec<StringRecord> = Vec::new();

    while std::mem::take(&mut peeked) || read_row(&mut reader, &mut carry, layout)? {
        if carry.get(5).is_none_or(|f| f.is_empty()) {
            continue; // If the 5th field is empty, it's a record we need to avoid deserializing.
        }

        let field = |index: usize| carry.get(index).unwrap_or_default();
        if options.skip_non_variants && field(5).eq_ignore_ascii_case(field(6)) {
            continue;
        }

        if field(0) != row_sample {
            let message = format!("Expected sample '{}' found '{}'!", row_sample, field(0));
            return Err(message.into());
        };

        let column = layout.segment_column();
        if region
            .last()
            .is_some_and(|last| last.get(column) != carry.get(column))
        {
            sink.write_region(&region, layout, options.fix_leaky_mnv)?;
            region.clear();
        }
        region.push(compact(&carry));
    }
    sink.write_region(&region, layout, options.fix_leaky_mnv)?;

    Ok(0)
}

/// Copies a row into a record holding only its fields and position, as a region's rows are held
/// until the region ends.
fn compact(record: &StringRecord) -> StringRecord {
    let mut copy = StringRecord::with_capacity(record.as_slice().len(), record.len());
    copy.extend(record.iter());
    copy.set_position(record.position().cloned());
    copy
}

/// Parses a row into its calls, one per sample with the tumor (or only) sample first and a default
/// call filling the second slot of a tumor-only row, and its tumor-normal VARDICT_STATUS.
fn parse_row(
    record: &StringRecord,
    layout: Layout,
) -> Result<([TumorOnlyVariant<'_>; 2], Option<&str>), String> {
    Ok(match layout {
        Layout::TumorOnly => {
            let call: TumorOnlyVariant = record.deserialize(None).map_err(describe_parse_error)?;
            ([call, TumorOnlyVariant::default()], None)
        }
        Layout::TumorNormal => {
            let row: TumorNormalVariant = record.deserialize(None).map_err(describe_parse_error)?;
            ([row.tumor_call(), row.normal_call()], Some(row.status))
        }
    })
}

/// Normalizes records and writes each call once.
struct Sink<'r> {
    writer: VcfWriter,
    reference: Reference,
    seen: SeenRows,
    filters: &'r FilterThresholds,
    progress: proglog::ProgLog,
}

impl Sink<'_> {
    /// Writes the rows of one VarDict region, a block it writes together for one target, each call
    /// once; VarDict orders a tumor-only region's positions by hash, so rows can come in any order.
    ///
    /// When asked, the partial rows of each MNV are folded into it first. Rows of the region that
    /// normalize to one call hold different reads of it, and the best supported is written, by
    /// [`Support`]. A call an earlier region wrote, or folded as a partial, is not written again,
    /// and is not folded into an MNV here, so its reads are never counted twice.
    fn write_region(
        &mut self,
        region: &[StringRecord],
        layout: Layout,
        fix_leaky_mnv: bool,
    ) -> Result<(), Box<dyn error::Error>> {
        let samples = layout.samples();
        let mut parsed: Vec<TumorOnlyVariant> = Vec::with_capacity(region.len() * samples);
        let mut statuses: Vec<Option<&str>> = Vec::with_capacity(region.len());
        for record in region {
            let (row_calls, status) = parse_row(record, layout)?;
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
        let folds = if fix_leaky_mnv {
            fold_partials(&calls, &written)
        } else {
            vec![Fold::Keep; calls.len()]
        };

        let mut best: AHashMap<CallKey, (usize, Support)> = AHashMap::new();
        for (row, fold) in folds.iter().enumerate() {
            let support = match fold {
                _ if written[row] => continue,
                Fold::Keep => Support::new(calls[row], &calls[row][0]),
                Fold::Merged(merged) => Support::new(&merged.calls, &calls[row][0]),
                Fold::Into(_) => continue,
            };
            let key = CallKey::new(calls[row][0].contig, &normalized[row]);
            let entry = best.entry(key).or_insert((row, support));
            if support < entry.1 {
                *entry = (row, support);
            }
        }

        for (row, fold) in folds.iter().enumerate() {
            let key = CallKey::new(calls[row][0].contig, &normalized[row]);
            if best.get(&key).is_none_or(|&(chosen, _)| chosen != row) {
                continue;
            }
            let (row_calls, merged) = match fold {
                Fold::Merged(merged) => (&merged.calls[..], Some(merged)),
                _ => (calls[row], None),
            };
            let row_calls: Vec<TumorOnlyVariant> = row_calls
                .iter()
                .map(|call| normalized[row].apply(call))
                .collect();
            write_record(
                &mut self.writer,
                &row_calls,
                statuses[row],
                self.filters,
                merged,
            )?;
            self.progress.record();
        }
        for (row, call) in calls.iter().zip(&normalized) {
            self.seen.insert(row[0].contig, call)?;
        }
        Ok(())
    }
}

/// How well a record supports its call, ordered best first: the most ALT reads in the tumor (or
/// only) sample, then in all samples, then the most depth in the tumor, then VarDict's lowest
/// position, REF and ALT for the row, so the choice never depends on the order of the rows.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Support<'a> {
    tumor_reads: Reverse<i64>,
    reads: Reverse<i64>,
    tumor_depth: Reverse<i64>,
    start: u64,
    ref_allele: &'a str,
    alt_allele: &'a str,
}

impl<'a> Support<'a> {
    /// The support of a record with these calls, written for a row with this tumor call.
    fn new(calls: &[TumorOnlyVariant], row: &TumorOnlyVariant<'a>) -> Self {
        let reads = |call: &TumorOnlyVariant| call.alt_depth.max(0) as i64;
        Support {
            tumor_reads: Reverse(reads(&calls[0])),
            reads: Reverse(calls.iter().map(reads).sum()),
            tumor_depth: Reverse(calls[0].depth.max(0) as i64),
            start: row.start,
            ref_allele: row.ref_allele,
            alt_allele: row.alt_allele,
        }
    }
}

/// A normalized call's contig, start, END, REF and ALT, which no two written records share.
#[derive(Hash, PartialEq, Eq)]
struct CallKey<'n> {
    contig: &'n str,
    start: u64,
    end: u64,
    ref_allele: &'n str,
    alt_allele: &'n str,
}

impl<'n> CallKey<'n> {
    /// The key of a normalized call on a contig.
    fn new(contig: &'n str, call: &'n Normalized) -> Self {
        CallKey {
            contig,
            start: call.start,
            end: call.end,
            ref_allele: &call.ref_allele,
            alt_allele: &call.alt_allele,
        }
    }
}

/// Writes one VCF record with a sample column per call, the tumor (or only) sample first; FILTER and
/// the site's INFO come from the first call, and a tumor-normal record carries its VARDICT_STATUS.
fn write_record(
    writer: &mut VcfWriter,
    calls: &[TumorOnlyVariant],
    status: Option<&str>,
    filters: &FilterThresholds,
    merged: Option<&Merged>,
) -> Result<(), Box<dyn error::Error>> {
    let var = &calls[0];
    let rid = writer.header().name2rid(var.contig.as_bytes()).unwrap();

    let mut variant = writer.empty_record();
    variant.set_rid(Some(rid));
    variant.set_pos(var.start as i64 - 1);
    variant.set_alleles(&[
        var.ref_allele.as_bytes(),
        var.alt_allele_for_vcf().as_bytes(),
    ])?;

    variant.set_qual(f32::missing());

    if filters.any() && var.has_alt_reads() {
        let labels = filters.labels(var);
        if labels.is_empty() {
            variant.push_filter("PASS".as_bytes())?;
        }
        for label in labels {
            variant.push_filter(label.as_bytes())?;
        }
    }

    if let Some(class) = var.variant_class() {
        variant.push_info_string(b"TYPE", &[class.as_bytes()])?;
    }

    if let Some(shift) = var.indel_3p_shift_value() {
        variant.push_info_integer(b"INDEL_3P_SHIFT", &[shift])?;
    }

    if let (Some(copies), Some(unit_length)) =
        (var.repeat_unit_copies_value(), var.repeat_unit_len_value())
    {
        variant.push_info_float(b"REPEAT_UNIT_COPIES", &[copies])?;
        variant.push_info_integer(b"REPEAT_UNIT_LEN", &[unit_length])?;
    }

    if VALID_SV_TYPES.contains(&var.variant_type) {
        variant.push_info_integer(b"END", &[var.end as i32])?;
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
        variant.push_info_integer(b"SVLEN", &[sv_length])?;
        variant.push_info_string(b"SVTYPE", &[var.variant_type.as_bytes()])?;
    }

    if let (Some(status), [tumor, normal]) = (status, calls) {
        variant.push_info_string(b"VARDICT_STATUS", &[status.as_bytes()])?;
        variant.push_info_float(
            b"TUMOR_NORMAL_FISHER_P",
            &[tumor_normal_fisher_p_value(tumor, normal)],
        )?;
    }

    if let Some(merged) = merged {
        variant.push_info_integer(b"MNV_PARTIALS", &[merged.partials as i32])?;
    }

    let genotypes = Values::<_, 4>::collect(
        GenotypeAllele::UnphasedMissing,
        calls
            .iter()
            .flat_map(|call| call.gt_value(MIN_HOM_ALT_AF).iter().copied()),
    );
    variant.push_genotypes(&genotypes)?;
    variant.push_format_integer(b"AD", &each_allele(calls, TumorOnlyVariant::ad_value))?;
    variant.push_format_integer(b"ADF", &each_allele(calls, TumorOnlyVariant::adf_value))?;
    variant.push_format_integer(b"ADR", &each_allele(calls, TumorOnlyVariant::adr_value))?;
    variant.push_format_float(
        b"STRAND_BIAS_FISHER_P",
        &each(calls, TumorOnlyVariant::strand_bias_fisher_p_value),
    )?;
    variant.push_format_integer(b"DP", &each(calls, |call| call.depth))?;
    variant.push_format_float(b"AF", &each(calls, TumorOnlyVariant::af_value))?;
    variant.push_format_integer(b"HICNT", &each(calls, TumorOnlyVariant::hicnt_value))?;
    variant.push_format_float(
        b"REALIGNED_FRAC_OF_DP",
        &each(calls, TumorOnlyVariant::realigned_frac_of_dp_value),
    )?;
    variant.push_format_float(
        b"MEAN_DIST_TO_READ_END",
        &each(calls, TumorOnlyVariant::mean_dist_to_read_end_value),
    )?;
    variant.push_format_integer(
        b"ALT_READ_POS_VARIES",
        &each(calls, TumorOnlyVariant::alt_read_pos_varies_value),
    )?;
    variant.push_format_float(b"QMEAN", &each(calls, TumorOnlyVariant::qmean_value))?;
    variant.push_format_float(
        b"MEAN_MAPQ",
        &each(calls, TumorOnlyVariant::mean_mapq_value),
    )?;
    variant.push_format_float(
        b"MEAN_MISMATCHES",
        &each(calls, TumorOnlyVariant::mean_mismatches_value),
    )?;

    if let Some(merged) = merged {
        variant.push_format_integer(b"MNV_PARTIAL_READS", &merged.partial_reads)?;
    }

    if let Some(counts) = calls
        .iter()
        .map(TumorOnlyVariant::sv_read_counts)
        .collect::<Option<Vec<_>>>()
    {
        variant.push_format_integer(b"SV_SOFTCLIP_READS", &each(&counts, |count| count.0))?;
        variant.push_format_integer(b"SV_DISCORDANT_READS", &each(&counts, |count| count.1))?;
    }

    writer.write(&variant)?;
    Ok(())
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

/// The calls written so far, each held as a fixed-size key, for skipping a call written before.
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

/// Reads the next row, refusing one whose layout differs from the first row's.
fn read_row<R: Read>(
    reader: &mut Reader<R>,
    record: &mut StringRecord,
    layout: Layout,
) -> Result<bool, Box<dyn error::Error>> {
    if !reader.read_record(record)? {
        return Ok(false);
    }
    if !is_segment(record, layout.segment_column()) {
        return Err(layout_error(record).into());
    }
    check_width(record, layout)?;
    Ok(true)
}

/// Refuses a row with fewer columns than its layout has.
fn check_width(record: &StringRecord, layout: Layout) -> Result<(), String> {
    if record.len() >= layout.min_columns() {
        return Ok(());
    }
    let line = record.position().map_or(0, |position| position.line());
    Err(format!(
        "Expected at least {} columns on line {line}, found {}!",
        layout.min_columns(),
        record.len()
    ))
}

/// Detects the layout of a row from the column holding its segment.
fn detect_layout(record: &StringRecord) -> Result<Layout, String> {
    [Layout::TumorOnly, Layout::TumorNormal]
        .into_iter()
        .find(|layout| is_segment(record, layout.segment_column()))
        .ok_or_else(|| layout_error(record))
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

/// Describes why a row matches neither layout.
fn layout_error(record: &StringRecord) -> String {
    let line = record.position().map_or(0, |position| position.line());
    if NO_FISHER_SEGMENT_COLUMNS
        .iter()
        .any(|&index| is_segment(record, index))
    {
        format!(
            "Expected the strand-bias p-value and odds-ratio columns on line {line}: run VarDictJava with --fisher!"
        )
    } else {
        let contig = record.get(2).unwrap_or_default();
        format!(
            "Expected a {contig}:start-end segment in column 35 (tumor-only) or 53 (tumor-normal) on line {line}!"
        )
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

/// Names the 1-based line and column of a record that could not be deserialized.
fn describe_parse_error(error: csv::Error) -> String {
    match error.kind() {
        csv::ErrorKind::Deserialize {
            pos: Some(pos),
            err,
        } => match err.field() {
            Some(index) => format!(
                "Could not parse column {} on line {}: {}!",
                index + 1,
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
        "tests/calls.tumor-normal.var",
        57,
        "Expected at least 61 columns on line 1, found 57!"
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

    #[rstest]
    #[case("tests/calls.no-fisher.reference.var")]
    #[case("tests/calls.no-fisher.variant.var")]
    fn test_refuses_rows_without_fisher_columns(
        #[case] path: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let input = BufReader::new(File::open(path)?);
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
            "Expected the strand-bias p-value and odds-ratio columns on line 1: run VarDictJava with --fisher!"
        );
        assert_eq!(output.as_file().metadata()?.len(), 0);
        Ok(())
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
            "Expected a chr1:start-end segment in column 35 (tumor-only) or 53 (tumor-normal) on line 1!"
        );
    }

    fn paired_row_without_fisher() -> String {
        let rows = std::fs::read_to_string("tests/calls.tumor-normal.var").unwrap();
        let mut fields: Vec<&str> = rows.lines().next().unwrap().split('\t').collect();
        fields.truncate(59);
        fields.remove(46);
        fields.remove(45);
        fields.remove(26);
        fields.remove(25);
        format!("{}\n", fields.join("\t"))
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

    #[test]
    fn test_detect_layout() {
        let tumor_only = std::fs::read_to_string("tests/calls.var").unwrap();
        let paired = std::fs::read_to_string("tests/calls.tumor-normal.var").unwrap();
        assert_eq!(
            detect_layout(&first_record(&tumor_only)),
            Ok(Layout::TumorOnly)
        );
        assert_eq!(
            detect_layout(&first_record(&paired)),
            Ok(Layout::TumorNormal)
        );
    }

    #[test]
    fn test_detect_layout_refuses_paired_rows_without_fisher_columns() {
        assert_eq!(
            paired_row_without_fisher().trim_end().split('\t').count(),
            55
        );
        assert_eq!(
            detect_layout(&first_record(&paired_row_without_fisher())),
            Err("Expected the strand-bias p-value and odds-ratio columns on line 1: run VarDictJava with --fisher!".to_string())
        );
    }

    #[test]
    fn test_detect_layout_refuses_rows_without_a_segment() {
        let input = "dna00001\tBRINP3\tchr1\t190098265\t190098265\tA\tA\n";
        assert_eq!(
            detect_layout(&first_record(input)),
            Err("Expected a chr1:start-end segment in column 35 (tumor-only) or 53 (tumor-normal) on line 1!".to_string())
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
        use rust_htslib::bcf::{Read, Reader as VcfReader};
        let reader = VcfReader::from_path(output.path()).expect("Error opening output file!");
        let samples: Vec<String> = reader
            .header()
            .samples()
            .iter()
            .map(|sample| String::from_utf8_lossy(sample).to_string())
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
                fix_leaky_mnv,
                ..Default::default()
            },
        )?;
        use rust_htslib::bcf::{Read, Reader as VcfReader};
        let mut reader = VcfReader::from_path(output.path()).expect("Error opening output file!");
        assert_eq!(
            reader.header().info_type(b"MNV_PARTIALS").is_ok(),
            fix_leaky_mnv
        );
        assert_eq!(
            reader.header().format_type(b"MNV_PARTIAL_READS").is_ok(),
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
        use rust_htslib::bcf::{Read, Reader as VcfReader};
        let mut reader = VcfReader::from_path(output.path()).expect("Error opening output file!");

        let mut record_count = 0;
        for record_result in reader.records() {
            let record = record_result?;
            let alleles = record.alleles();

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
