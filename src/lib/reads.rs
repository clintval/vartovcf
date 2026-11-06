//! Follows a partial row's reads to its MNV's other changed bases, to tell the MNV's own molecules
//! from another allele's.
//!
//! VarDict counts a read toward a partial row when it cannot grow the MNV through the read: the read
//! ends inside the MNV, or has a no-call (N) or a base below `-q` plus 5 at a changed base. The
//! read's template usually still holds that base, in the read itself at a lower quality or in its
//! mate, and that base says whose molecule the read is. A template with the MNV's ALT base at every
//! other changed position, at a good quality, is the MNV's own molecule; one with a good-quality base
//! that is not the ALT at any of them is another allele's; one with no good-quality base at some of
//! them, a no-call, a base below `-q` or no coverage, is consistent with leakage and proves nothing.
//!
//! The reads of a sample are piled up over the MNV's span with `streampile`, from an indexed query,
//! and each read's walk through the span replays VarDict's `CigarParser` for substitutions, so a
//! read feeds exactly the rows VarDict counted it toward: the partial its growth stopped at, or the
//! MNV. Each template with a read feeding a partial is then classified by its called bases, the
//! bases of its reads called into one as fgbio's `CallOverlappingConsensusBases` calls them, at the
//! MNV's changed positions outside the partial.
//!
//! Each read is also walked as VarDict clips it, and a partial splits when either walk numbers VarDict's.
//!
//! The same pileup recounts an SNV's high-quality ALT reads: those carrying the ALT base aligned at
//! the site, and those whose soft clip starts there with the ALT base and matches the reference
//! beyond it. VarDict's `adjSNV` adds a soft clip of up to 5 bases to an SNV after comparing only
//! its second base to the reference, so a clip whose later bases are foreign lifts a one-read SNV
//! over its `-r` floor; the recount leaves such clips out.
//!
//! The pileup also looks for the reads that carry an insertion whole, as an insertion operation of
//! its bases: an aligner soft-clips a read at an insertion too long to open as a gap, so VarDict
//! assembles such an insertion from the soft clips at its ends (`VariationRealigner.realignlgins`)
//! and no read holds it. A read whose insertion operation sits where the insertion slides to, in
//! a repeat, carries it too.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use bstr::BString;
use log::debug;
use noodles::bam;
use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
use noodles::sam::header::record::value::map::read_group::tag::SAMPLE;
use noodles::sam::header::record::value::{Map, map};
use noodles::sam::{self};
use streampile::{
    AgreementStrategy, DisagreementStrategy, EntryKind, RecordSource, StreamingPileupBuilder,
};

use crate::mnv::{ReadFate, ReadFates, is_sub_event, is_substitution};
use crate::record::TumorOnlyVariant;
use crate::reference::Reference;

/// The indexed BAMs VarDict read, one per sample with the tumor (or only) sample's first, and how
/// VarDict read them.
#[derive(Clone, Debug, PartialEq)]
pub struct ReadFiles {
    /// The BAMs, each indexed beside it.
    pub bams: Vec<PathBuf>,
    /// How VarDict read them.
    pub options: ReadOptions,
}

/// How VarDict read the BAM, so the reads followed are the reads it counted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReadOptions {
    /// The quality a base needs to be a good call, VarDict's `-q`: an MNV grows through a changed
    /// base only at this plus 5, and a base at another changed position is evidence only at this.
    pub base_quality: f64,
    /// The mapping quality a read needs, VarDict's `-Q`.
    pub min_mapping_quality: u8,
    /// The flags that leave a read out, VarDict's `-F`.
    pub exclude_flags: u16,
    /// The most mismatches a read may have, VarDict's `-m`: its NM minus its inserted and deleted
    /// bases.
    pub max_mismatches: i64,
    /// The most reference-matching bases VarDict bridges between an MNV's mismatches, its `-X`.
    pub mnv_gap: usize,
}

impl Default for ReadOptions {
    /// VarDict's own defaults.
    fn default() -> Self {
        ReadOptions {
            base_quality: 22.5,
            min_mapping_quality: 0,
            exclude_flags: 0x504,
            max_mismatches: 8,
            mnv_gap: 2,
        }
    }
}

impl ReadOptions {
    /// Whether a base quality is a good call.
    fn is_good(self, quality: u8) -> bool {
        f64::from(quality) >= self.base_quality
    }

    /// Whether a base quality lets an MNV grow through the base.
    fn grows(self, quality: u8) -> bool {
        f64::from(quality) >= self.base_quality + 5.0
    }

    /// The least integer quality that is a good call.
    fn floor(self) -> u8 {
        self.base_quality.ceil().clamp(0.0, 255.0) as u8
    }

    /// The bases on each side of an MNV that a read's walk must cover: a mismatch chain that reaches
    /// the MNV has a link within this many bases of it, and VarDict's clip scan stops within them.
    fn flank(self) -> usize {
        (self.mnv_gap + 1).max(CLIP_SCAN_MATCHES)
    }
}

/// How many partial rows the reads decided, and how many of them reconciled with VarDict's counts.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    /// MNV rows whose span was piled up, counted once per sample.
    pub pileups: usize,
    /// Partial rows decided, counted once per sample with reads of the row.
    pub decisions: usize,
    /// Decisions where the reads feeding the row number exactly VarDict's ALT reads.
    pub reconciled: usize,
    /// Reconciled decisions whose reads numbered VarDict's only once its end clipping was replayed.
    pub reconciled_clipped: usize,
    /// Decisions the reads could change: the row has fewer reads than the MNV in the sample and
    /// some read found would move.
    pub decisive: usize,
    /// Decisive decisions whose reads did not number VarDict's, so the sample kept every read.
    pub decisive_unreconciled: usize,
    /// Decisions where no read in the pileup fed the row VarDict counted reads toward.
    pub empty: usize,
}

/// How many SNV rows were recounted from the tumor (or only) sample's reads, and what was found.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RecountSummary {
    /// SNV rows recounted.
    pub rows: usize,
    /// Rows found short of the floor they were recounted against.
    pub short: usize,
    /// Rows no read in the pileup supports at all.
    pub unsupported: usize,
}

/// The high-quality ALT reads of an SNV found in one sample's reads, by how each carries the ALT
/// base: aligned at the site, or first in a soft clip at the site that matches the reference beyond
/// it, which both support the call; or first in a soft clip that does not, which VarDict counted
/// anyway when the clip's second base matched.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Recount {
    /// Reads with the ALT base aligned at the site, once VarDict's end clipping is replayed.
    pub aligned: i32,
    /// Reads whose soft clip starts at the site with the ALT base and matches the reference past it.
    pub clean_clips: i32,
    /// Reads whose soft clip starts at the site with the ALT base and does not match the reference
    /// past it, though its second base does.
    pub foreign_clips: i32,
}

impl Recount {
    /// The reads that support the call: aligned, or in a clip that matches the reference.
    pub fn supporting(&self) -> i32 {
        self.aligned + self.clean_clips
    }
}

/// What the reads say about one call, for the labels that ask them: an SNV's high-quality ALT
/// reads recounted, and the reads carrying an insertion whole; `None` where the label is not
/// applied, the call is not one it tests, or no reads were given.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ReadEvidence {
    /// An SNV's recount, for FOREIGN_CLIP_SUPPORT.
    pub recount: Option<Recount>,
    /// The reads carrying an insertion whole as an insertion operation, for UNSPANNED_INSERTION.
    pub spanning: Option<i32>,
}

/// How many insertion rows were looked for in the tumor (or only) sample's reads, and how many a
/// read carried whole.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SpanningSummary {
    /// Insertion rows looked for.
    pub rows: usize,
    /// Rows some read carried whole.
    pub spanned: usize,
}

/// One sample's indexed BAM.
struct Sample {
    path: PathBuf,
    reader: bam::io::IndexedReader<bgzf::io::Reader<File>>,
    header: sam::Header,
}

/// The indexed BAMs VarDict read, one per sample, the tumor (or only) sample first.
pub struct ReadSources {
    samples: Vec<Sample>,
    options: ReadOptions,
    /// What the reads decided so far.
    pub summary: Summary,
    /// What the recounts found so far.
    pub recounts: RecountSummary,
    /// What the searches for reads carrying an insertion whole found so far.
    pub spans: SpanningSummary,
}

impl ReadSources {
    /// Opens each BAM and its index, which must be beside it as `<bam>.bai` or `<bam>.csi`.
    ///
    /// An indexed BAM is in coordinate order whatever its header says, so a header without that
    /// order is given it.
    pub fn open(paths: &[PathBuf], options: ReadOptions) -> Result<Self, String> {
        let samples = paths
            .iter()
            .map(|path| {
                let describe = |error| {
                    format!(
                        "Could not open the reads {path:?} or their index (.bai or .csi beside them): {error}!"
                    )
                };
                let mut reader = bam::io::indexed_reader::Builder::default()
                    .build_from_path(path)
                    .map_err(describe)?;
                let mut header = reader.read_header().map_err(describe)?;
                header
                    .header_mut()
                    .get_or_insert_with(|| Map::<map::Header>::new(Default::default()))
                    .other_fields_mut()
                    .insert(SORT_ORDER, BString::from("coordinate"));
                Ok(Sample {
                    path: path.clone(),
                    reader,
                    header,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        Ok(ReadSources {
            samples,
            options,
            summary: Summary::default(),
            recounts: RecountSummary::default(),
            spans: SpanningSummary::default(),
        })
    }

    /// The number of BAMs, one per sample.
    pub fn samples(&self) -> usize {
        self.samples.len()
    }

    /// Warnings for BAMs that look swapped, given the sample names in BAM order: a BAM whose read
    /// groups name none of its own sample's names but do name another sample's, while the BAM
    /// standing for that other sample does not name it either. A read group sample that is simply
    /// another name for the sample is logged, not warned about, as pipelines name the two freely,
    /// and two BAMs naming one sample are not told apart.
    pub fn sample_warnings(&self, names: &[&str]) -> Vec<String> {
        let read_group_samples: Vec<Vec<String>> = self
            .samples
            .iter()
            .map(|sample| {
                sample
                    .header
                    .read_groups()
                    .values()
                    .filter_map(|group| group.other_fields().get(&SAMPLE))
                    .map(|sample| sample.to_string())
                    .collect()
            })
            .collect();
        let names_sample = |index: usize, name: &str| {
            read_group_samples[index]
                .iter()
                .any(|sample| sample == name)
        };
        self.samples
            .iter()
            .zip(names)
            .enumerate()
            .filter_map(|(index, (sample, name))| {
                if read_group_samples[index].is_empty() || names_sample(index, name) {
                    return None;
                }
                let swapped = names.iter().enumerate().find(|&(position, other)| {
                    position != index
                        && position < read_group_samples.len()
                        && names_sample(index, other)
                        && !names_sample(position, other)
                });
                match swapped {
                    Some((position, other)) => Some(format!(
                        "The reads {} name the sample '{other}' in their read groups, the sample of the BAM at position {}, not '{name}': give --reads the BAMs in sample order, the tumor (or only) sample's first!",
                        sample.path.display(),
                        position + 1
                    )),
                    None => {
                        debug!(
                            "The reads {} name the sample(s) {} in their read groups; they stand for '{name}'.",
                            sample.path.display(),
                            read_group_samples[index].join(", ")
                        );
                        None
                    }
                }
            })
            .collect()
    }

    /// Decides the fate of every partial row of a region under every MNV row it is a sub-event of,
    /// from the reads of each sample, keyed by (MNV row, partial row). Rows already written are not
    /// partials here, as they are never folded again. The fold applies a fate only when the MNV
    /// claims the row, so a row visited before its MNV, as a better supported one is, keeps every
    /// read whatever its fate says.
    pub fn fates<'a>(
        &mut self,
        rows: &[&[TumorOnlyVariant<'a>]],
        written: &[bool],
        reference: &mut Reference,
    ) -> Result<ReadFates<'a>, String> {
        let options = self.options;
        let mut fates = ReadFates::new();
        let with_reads = |row: usize| {
            is_substitution(&rows[row][0]) && rows[row].iter().any(|call| call.alt_depth > 0)
        };
        for mnv in
            (0..rows.len()).filter(|&row| with_reads(row) && rows[row][0].ref_allele.len() > 1)
        {
            let partials: Vec<usize> = (0..rows.len())
                .filter(|&row| {
                    !written[row] && with_reads(row) && is_sub_event(&rows[mnv][0], &rows[row][0])
                })
                .collect();
            if partials.is_empty() {
                continue;
            }
            let window = Window::around(&rows[mnv][0], options, reference)?;
            let surveys = (0..rows[mnv].len())
                .map(|sample| self.survey_sample(sample, &window))
                .collect::<Result<Vec<Survey>, String>>()?;
            for row in partials {
                let partial_calls: Vec<&TumorOnlyVariant<'a>> = rows[row].iter().collect();
                let classified: Vec<Classified> = surveys
                    .iter()
                    .zip(&partial_calls)
                    .map(|(survey, partial)| {
                        classify(survey, &window, &rows[mnv][0], partial, options)
                    })
                    .collect();
                let (fate, outcomes) = resolve(&classified, rows[mnv], &partial_calls, options);
                for (sample, (classified, outcome)) in classified.iter().zip(&outcomes).enumerate()
                {
                    let partial = partial_calls[sample];
                    if partial.alt_depth == 0 {
                        continue;
                    }
                    self.summary.decisions += 1;
                    let reconciled = classified.decision.reconciled();
                    self.summary.reconciled += usize::from(reconciled);
                    self.summary.reconciled_clipped +=
                        usize::from(reconciled && classified.decision.unmatched != 0);
                    self.summary.empty += usize::from(classified.feeding.is_empty());
                    if partial.alt_depth < rows[mnv][sample].alt_depth && outcome.would_move > 0 {
                        self.summary.decisive += 1;
                        self.summary.decisive_unreconciled += usize::from(!reconciled);
                    }
                    debug!(
                        "{contig}:{} {}>{} under {contig}:{} {}>{} in sample {}: VarDict counted {} reads; {} are the MNV's molecules, {} another allele's, {} prove nothing, {} unmatched; the reads would move {} and keep {}; {:.3} of the templates prove nothing there; {} unmatched once VarDict's end clipping is replayed.",
                        partial.start,
                        partial.ref_allele,
                        partial.alt_allele,
                        rows[mnv][0].start,
                        rows[mnv][0].ref_allele,
                        rows[mnv][0].alt_allele,
                        sample + 1,
                        partial.alt_depth,
                        classified.decision.proven,
                        classified.decision.separate,
                        classified.decision.ambiguous,
                        classified.decision.unmatched,
                        outcome.moved,
                        outcome.kept,
                        classified.decision.ambiguity,
                        classified.decision.unmatched_clipped,
                        contig = partial.contig,
                    );
                }
                fates.insert((mnv, row), fate);
            }
        }
        Ok(fates)
    }

    /// Piles up one sample's reads over a window.
    fn survey_sample(&mut self, sample: usize, window: &Window) -> Result<Survey, String> {
        let options = self.options;
        let source = &mut self.samples[sample];
        let contig = window.contig.as_str();
        let (start, end) = (window.start + 1, window.start + window.reference.len());
        let describe = |error: &dyn std::fmt::Display| {
            format!(
                "Could not read the reads of {} over {contig}:{start}-{end}: {error}!",
                source.path.display()
            )
        };
        let header = &source.header;
        if header
            .reference_sequences()
            .get(contig.as_bytes())
            .is_none()
        {
            return Err(format!(
                "Contig '{contig}' is not in the reads {}: give --reads the BAMs VarDict was run on!",
                source.path.display()
            ));
        }
        let region = Region::new(
            contig,
            Position::try_from(start).map_err(|error| error.to_string())?
                ..=Position::try_from(end).map_err(|error| error.to_string())?,
        );
        let query = source
            .reader
            .query(header, &region)
            .map_err(|error| describe(&error))?;
        let survey =
            survey(Query(query), header, window, options).map_err(|error| describe(&error))?;
        self.summary.pileups += 1;
        Ok(survey)
    }

    /// Recounts an SNV's high-quality ALT reads from the tumor (or only) sample's reads, noting in
    /// the summary whether they fall short of `min_hicnt`.
    pub fn recount(
        &mut self,
        call: &TumorOnlyVariant,
        min_hicnt: i32,
        reference: &mut Reference,
    ) -> Result<Recount, String> {
        let options = self.options;
        let window = Window::around(call, options, reference)?;
        let survey = self.survey_sample(0, &window)?;
        let recount = recount(&survey, &window, call, options, reference)?;
        self.recounts.rows += 1;
        self.recounts.short += usize::from(recount.supporting() < min_hicnt);
        self.recounts.unsupported += usize::from(recount.supporting() == 0);
        debug!(
            "{}:{} {}>{} in sample 1: VarDict counted {} high-quality ALT reads; the reads hold {} aligned, {} first in a soft clip that matches the reference and {} first in one that does not.",
            call.contig,
            call.start,
            call.ref_allele,
            call.alt_allele,
            call.hicnt_value(),
            recount.aligned,
            recount.clean_clips,
            recount.foreign_clips,
        );
        Ok(recount)
    }

    /// Counts the reads of the tumor (or only) sample that carry an insertion whole: those with an
    /// insertion operation of its bases, at or above `-q`, at its position or one it slides to.
    pub fn spanning(
        &mut self,
        call: &TumorOnlyVariant,
        reference: &mut Reference,
    ) -> Result<i32, String> {
        let options = self.options;
        let window = Window::sliding(call, options, reference)?;
        let survey = self.survey_sample(0, &window)?;
        let mut spanning = 0;
        for read in &survey.reads {
            spanning += i32::from(read.carries(call, options, reference)?);
        }
        self.spans.rows += 1;
        self.spans.spanned += usize::from(spanning > 0);
        debug!(
            "{}:{} {}>{} in sample 1: VarDict counted {} high-quality ALT reads; {} read(s) carry the insertion whole.",
            call.contig,
            call.start,
            call.ref_allele,
            call.alt_allele,
            call.hicnt_value(),
            spanning,
        );
        Ok(spanning)
    }
}

/// The records an indexed query returns, read one at a time into a builder's record.
struct Query<'r, R>(bam::io::reader::Query<'r, R>);

impl<R: bgzf::io::BufRead + bgzf::io::Seek> RecordSource for Query<'_, R> {
    type Record = bam::Record;

    fn read_record(&mut self, record: &mut bam::Record) -> io::Result<bool> {
        self.0.read_record(record).map(|read| read > 0)
    }
}

/// The reference bases over the span piled up, starting at a 0-based position.
struct Window {
    contig: String,
    start: usize,
    reference: Vec<u8>,
}

impl Window {
    /// The window around a row: its span plus a flank on each side, within the contig.
    fn around(
        row: &TumorOnlyVariant,
        options: ReadOptions,
        reference: &mut Reference,
    ) -> Result<Self, String> {
        let length = reference.length(row.contig)?;
        let start = (row.start - 1).saturating_sub(options.flank() as u64);
        let end = (row.end + options.flank() as u64).min(length);
        Ok(Window {
            contig: row.contig.to_string(),
            start: start as usize,
            reference: reference.bases(row.contig, start, end)?.to_vec(),
        })
    }

    /// The window around an indel row's position and every position its indel slides to, its 3'
    /// shift, plus a flank on each side, within the contig, so a read placing the indel anywhere
    /// along the slide is piled up.
    fn sliding(
        row: &TumorOnlyVariant,
        options: ReadOptions,
        reference: &mut Reference,
    ) -> Result<Self, String> {
        let length = reference.length(row.contig)?;
        let shift = row.num_bases_3_prime_shift_for_deletions.max(0) as u64;
        let start = (row.start - 1).saturating_sub(options.flank() as u64);
        let end = (row.end + shift + options.flank() as u64).min(length);
        Ok(Window {
            contig: row.contig.to_string(),
            start: start as usize,
            reference: reference.bases(row.contig, start, end)?.to_vec(),
        })
    }

    /// The column of a 1-based position, if it is in the window.
    fn column(&self, position: u64) -> Option<usize> {
        let column = position.checked_sub(self.start as u64 + 1)? as usize;
        (column < self.reference.len()).then_some(column)
    }

    /// The 1-based position of a column, outside the window included.
    fn position(&self, column: i64) -> i64 {
        self.start as i64 + column + 1
    }
}

/// One read's base at one column.
#[derive(Clone, Copy, Debug)]
struct Base {
    base: u8,
    quality: u8,
    query_position: usize,
}

/// What one read holds over the window.
#[derive(Debug)]
struct ReadView {
    template: usize,
    bases: Vec<Option<Base>>,
    /// The columns where the read has an insertion, a deletion or a skip.
    indels: Vec<usize>,
    reverse: bool,
    mapping_quality: u8,
    mismatches: f64,
    aligned_length: usize,
    leading_clip: usize,
    /// The window columns of the alignment's first and last base, outside the window included.
    span: Option<(i64, i64)>,
    /// Whether the CIGAR starts and ends with a match operation, the only ends VarDict clips.
    starts_with_match: bool,
    ends_with_match: bool,
    /// The bases the aligner soft-clipped at each end, with their qualities, outward from the
    /// alignment.
    head_clip: Vec<(u8, u8)>,
    tail_clip: Vec<(u8, u8)>,
    /// The read's insertion operations.
    insertions: Vec<Inserted>,
    /// The alleles the read feeds as aligned, and once VarDict's end clipping is replayed.
    aligned: Vec<Found>,
    clipped: Vec<Found>,
}

/// One insertion operation of a read: the 1-based position of the reference base before the
/// inserted bases, the bases, and their mean quality.
#[derive(Clone, Debug, PartialEq)]
struct Inserted {
    anchor: u64,
    bases: Vec<u8>,
    quality: f64,
}

/// The bases VarDict's `CigarModifier` soft-clips from a read's ends before walking it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Clips {
    head: usize,
    tail: usize,
}

/// The matches in a row that stop VarDict's scan for a mismatch near a read's end.
const CLIP_SCAN_MATCHES: usize = 3;
/// The most bases from a read's end VarDict soft-clips a mismatch at.
const CLIP_MAX: usize = 3;

/// What one template's reads, called into one, hold over the window: `None` where no read of the
/// template holds a base that votes, with the called base and quality otherwise, and whether any
/// read of the template is present at all.
#[derive(Debug, Default)]
struct TemplateView {
    called: Vec<Option<(u8, u8)>>,
    present: Vec<bool>,
}

/// The reads and templates of a sample over a window.
#[derive(Debug, Default)]
struct Survey {
    reads: Vec<ReadView>,
    templates: Vec<TemplateView>,
}

/// Piles up a sample's reads over the window and records what each read and each template holds
/// at each column, filtering reads as VarDict does: by the flags, and every supplementary record
/// whenever a flag filter is set at all, by mapping quality, and by mismatches.
fn survey<S: RecordSource<Record = bam::Record>>(
    source: S,
    header: &sam::Header,
    window: &Window,
    options: ReadOptions,
) -> streampile::Result<Survey> {
    let columns = window.reference.len();
    let max_mismatches = options.max_mismatches;
    let mut exclude_flags = Flags::from_bits_retain(options.exclude_flags);
    if !exclude_flags.is_empty() {
        exclude_flags |= Flags::SUPPLEMENTARY;
    }
    let mut builder = StreamingPileupBuilder::new(source, header)?
        .min_mapping_quality(options.min_mapping_quality)
        .exclude_flags(exclude_flags)
        .min_base_quality(options.floor())
        .read_filter(move |record: &bam::Record| {
            substitution_mismatches(record).is_none_or(|mismatches| mismatches <= max_mismatches)
        });
    let mut survey = Survey::default();
    let mut templates: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut reads: HashMap<(usize, u16, i64), usize> = HashMap::new();
    let mut sweep = builder.columns(&window.contig, window.start, window.start + columns)?;
    for column in 0..columns {
        let Some(pileup) = sweep.next_pileup() else {
            break;
        };
        let pileup = pileup?;
        for entry in pileup.iter() {
            let record = entry.record();
            let name: &[u8] = record.name().map_or(&b"*"[..], |name| name.as_ref());
            let template = match templates.get(name) {
                Some(&template) => template,
                None => {
                    templates.insert(name.to_vec(), survey.templates.len());
                    survey.templates.push(TemplateView {
                        called: vec![None; columns],
                        present: vec![false; columns],
                    });
                    survey.templates.len() - 1
                }
            };
            survey.templates[template].present[column] = true;
            let start = record
                .alignment_start()
                .and_then(Result::ok)
                .map_or(-1, |start| usize::from(start) as i64);
            let key = (template, record.flags().bits(), start);
            let next = survey.reads.len();
            let read = *reads.entry(key).or_insert(next);
            if read == survey.reads.len() {
                survey
                    .reads
                    .push(ReadView::new(record, template, window.start, columns));
            }
            let view = &mut survey.reads[read];
            match entry.kind() {
                EntryKind::Base => {
                    if let (Some(base), Some(quality), Some(query_position)) =
                        (entry.base(), entry.quality(), entry.query_position())
                    {
                        view.bases[column] = Some(Base {
                            base: base.to_ascii_uppercase(),
                            quality,
                            query_position,
                        });
                    }
                }
                EntryKind::Insertion | EntryKind::Deletion | EntryKind::Skip => {
                    view.indels.push(column)
                }
                _ => {}
            }
        }
        for template in pileup.templates(
            AgreementStrategy::Consensus,
            DisagreementStrategy::Consensus,
        ) {
            if let Some(&index) = templates.get(template.name().as_ref() as &[u8]) {
                survey.templates[index].called[column] = template.base().zip(template.quality());
            }
        }
    }
    builder.close()?;
    for read in &mut survey.reads {
        read.walk_both(window, options);
    }
    Ok(survey)
}

/// A read's NM minus its inserted and deleted bases, as VarDict counts its mismatches, or `None`
/// for a read without an NM tag, which VarDict keeps.
fn substitution_mismatches(record: &bam::Record) -> Option<i64> {
    let nm = record
        .data()
        .get(&Tag::EDIT_DISTANCE)
        .and_then(|value| value.ok())
        .and_then(|value| value.as_int())?;
    let indels: i64 = record
        .cigar()
        .iter()
        .filter_map(|op| op.ok())
        .filter(|op| matches!(op.kind(), Kind::Insertion | Kind::Deletion))
        .map(|op| op.len() as i64)
        .sum();
    Some(nm - indels)
}

impl ReadView {
    /// A read with no bases recorded yet, with what VarDict's per-read measures need from it, over
    /// a window starting at a 0-based position.
    fn new(record: &bam::Record, template: usize, window_start: usize, columns: usize) -> Self {
        let sequence = record.sequence();
        let qualities = record.quality_scores();
        let clip = |range: std::ops::Range<usize>| -> Vec<(u8, u8)> {
            sequence
                .iter()
                .zip(qualities.iter())
                .skip(range.start)
                .take(range.len())
                .collect()
        };
        let alignment_start = record.alignment_start().and_then(Result::ok);
        let mut leading_clip = 0;
        let mut trailing_clip = 0;
        let mut aligned_length = 0;
        let mut reference_span = 0;
        let mut query_position = 0;
        let mut before_alignment = true;
        let mut first_op = None;
        let mut last_op = None;
        let mut insertions = Vec::new();
        for op in record.cigar().iter().filter_map(|op| op.ok()) {
            first_op.get_or_insert(op.kind());
            last_op = Some(op.kind());
            match op.kind() {
                Kind::SoftClip if before_alignment => leading_clip += op.len(),
                Kind::SoftClip => trailing_clip += op.len(),
                Kind::HardClip | Kind::Pad => {}
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion => {
                    before_alignment = false;
                    aligned_length += op.len();
                    if op.kind() != Kind::Insertion {
                        reference_span += op.len();
                    } else if let Some(start) = alignment_start {
                        let inserted = clip(query_position..query_position + op.len());
                        let quality: f64 = inserted
                            .iter()
                            .map(|&(_, quality)| f64::from(quality))
                            .sum();
                        insertions.push(Inserted {
                            anchor: (usize::from(start) + reference_span) as u64 - 1,
                            bases: inserted
                                .iter()
                                .map(|&(base, _)| base.to_ascii_uppercase())
                                .collect(),
                            quality: quality / op.len() as f64,
                        });
                    }
                }
                Kind::Deletion | Kind::Skip => {
                    before_alignment = false;
                    reference_span += op.len();
                }
            }
            if op.kind().consumes_read() {
                query_position += op.len();
            }
        }
        let span = alignment_start.filter(|_| reference_span > 0).map(|start| {
            let first = usize::from(start) as i64 - 1 - window_start as i64;
            (first, first + reference_span as i64 - 1)
        });
        let mut head_clip = clip(0..leading_clip);
        head_clip.reverse();
        let tail_clip = clip(sequence.len().saturating_sub(trailing_clip)..sequence.len());
        ReadView {
            template,
            bases: vec![None; columns],
            indels: Vec::new(),
            reverse: record.flags().is_reverse_complemented(),
            mapping_quality: record.mapping_quality().map_or(255, u8::from),
            mismatches: substitution_mismatches(record).unwrap_or(0) as f64,
            aligned_length,
            leading_clip,
            span,
            starts_with_match: first_op == Some(Kind::Match),
            ends_with_match: last_op == Some(Kind::Match),
            head_clip,
            tail_clip,
            insertions,
            aligned: Vec::new(),
            clipped: Vec::new(),
        }
    }

    /// Whether the read carries an insertion row's inserted bases whole: in an insertion operation
    /// whose bases' mean quality reaches `-q`, at the row's position or one the insertion slides to.
    fn carries(
        &self,
        row: &TumorOnlyVariant,
        options: ReadOptions,
        reference: &mut Reference,
    ) -> Result<bool, String> {
        let inserted = row.alt_allele.as_bytes()[row.ref_allele.len()..].to_ascii_uppercase();
        for insertion in &self.insertions {
            if insertion.quality >= options.base_quality
                && same_insertion(
                    row.contig,
                    (row.start, &inserted),
                    (insertion.anchor, &insertion.bases),
                    reference,
                )?
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// The soft clip starting at a column, if the read has one there: VarDict's own, from a
    /// mismatch near a match-operation end, or the aligner's, after VarDict extends the alignment
    /// through its clipped bases that match the reference at a quality above 10.
    fn clip_at(
        &self,
        column: usize,
        window: &Window,
        reference: &mut Reference,
    ) -> Result<Option<Clip>, String> {
        let Some((start, end)) = self.span else {
            return Ok(None);
        };
        let site = column as i64;
        let clips = self.clips(window);
        let held = |columns: Vec<usize>| -> Vec<(u8, u8)> {
            columns
                .into_iter()
                .filter_map(|column| self.bases[column])
                .map(|base| (base.base, base.quality))
                .collect()
        };
        if clips.tail > 0 {
            if end - clips.tail as i64 + 1 == site {
                let bases = held((column..=end as usize).collect());
                return Ok(Some(Clip {
                    bases,
                    direction: 1,
                }));
            }
        } else if !self.tail_clip.is_empty() && end < site {
            let innermost = window.position(end + 1);
            let extended = extension(&self.tail_clip, &window.contig, innermost, 1, reference)?;
            if end + 1 + extended as i64 == site {
                return Ok(Some(Clip {
                    bases: self.tail_clip[extended..].to_vec(),
                    direction: 1,
                }));
            }
        }
        if clips.head > 0 {
            if start + clips.head as i64 - 1 == site {
                let bases = held((start as usize..=column).rev().collect());
                return Ok(Some(Clip {
                    bases,
                    direction: -1,
                }));
            }
        } else if !self.head_clip.is_empty() && start > site {
            let innermost = window.position(start - 1);
            let extended = extension(&self.head_clip, &window.contig, innermost, -1, reference)?;
            if start - 1 - extended as i64 == site {
                return Ok(Some(Clip {
                    bases: self.head_clip[extended..].to_vec(),
                    direction: -1,
                }));
            }
        }
        Ok(None)
    }

    /// The read's first and last columns within the window, if any.
    fn bounds(&self) -> Option<(usize, usize)> {
        let (first, last) = self.span?;
        let columns = self.bases.len() as i64;
        (last >= 0 && first < columns)
            .then(|| (first.max(0) as usize, last.min(columns - 1) as usize))
    }

    /// Walks the read as aligned and as VarDict clips it, keeping the alleles each walk feeds.
    fn walk_both(&mut self, window: &Window, options: ReadOptions) {
        let clips = self.clips(window);
        self.aligned = self.walk(window, options, Clips::default());
        self.clipped = if clips == Clips::default() {
            self.aligned.clone()
        } else {
            self.walk(window, options, clips)
        };
    }

    /// The bases VarDict soft-clips: a mismatch within three of a match-operation end, 3' end first.
    fn clips(&self, window: &Window) -> Clips {
        let (Some((first, last)), Some((start, end))) = (self.bounds(), self.span) else {
            return Clips::default();
        };
        if !self.indels.is_empty() {
            return Clips::default();
        }
        let columns = self.bases.len() as i64;
        let scan = |iter: &mut dyn Iterator<Item = usize>, ran_out: bool| -> usize {
            let mut clip = 0;
            let mut matches = 0;
            let mut scanned = 0;
            for column in iter {
                scanned += 1;
                let Some(base) = self.bases[column] else {
                    return 0;
                };
                if base.base != window.reference[column] {
                    clip = scanned;
                    matches = 0;
                } else {
                    matches += 1;
                    if matches == CLIP_SCAN_MATCHES {
                        return if clip <= CLIP_MAX { clip } else { 0 };
                    }
                }
            }
            if ran_out || clip > CLIP_MAX { 0 } else { clip }
        };
        let tail = if self.ends_with_match && end < columns {
            scan(&mut (first..=last).rev(), start < 0)
        } else {
            0
        };
        let head = if self.starts_with_match && start >= 0 && first + tail <= last {
            scan(&mut (first..=last - tail), end >= columns)
        } else {
            0
        };
        Clips { head, tail }
    }

    /// The alleles VarDict counts this read toward over the window, by replaying its walk: an
    /// allele starts at each mismatch and grows through the following bases as an MNV while they
    /// are contiguous in the read, not no-calls, at least `-q` plus 5, and either mismatches or
    /// matches bridged within the gap to a mismatch. An allele with an insertion, deletion or skip
    /// within the gap of either end feeds no substitution row, as VarDict folds it into a complex
    /// allele with the indel. The bases under the clips are left out, as soft-clipped bases are.
    fn walk(&self, window: &Window, options: ReadOptions, clips: Clips) -> Vec<Found> {
        let mut found = Vec::new();
        let reference = &window.reference;
        let columns = self.bases.len();
        let (first, last) = self.bounds().unwrap_or((0, 0));
        let at = |column: usize| -> Option<Base> {
            let clipped = clips != Clips::default()
                && (first..=last).contains(&column)
                && (column < first + clips.head || column + clips.tail > last);
            if clipped { None } else { self.bases[column] }
        };
        let aligned_length = self
            .aligned_length
            .saturating_sub(clips.head)
            .saturating_sub(clips.tail);
        let mut column = 0;
        while column < columns {
            let Some(first) = at(column) else {
                column += 1;
                continue;
            };
            if first.base == b'N' || reference[column] == b'N' || first.base == reference[column] {
                column += 1;
                continue;
            }
            let mut alt = vec![first.base];
            let mut quality = f64::from(first.quality);
            let mut last = column;
            let mut appended_mismatches = 0;
            let mut dropped = false;
            let contiguous = |from: usize, to: usize| -> Option<Base> {
                let base = at(to)?;
                (base.query_position == at(from)?.query_position + (to - from)).then_some(base)
            };
            if options.is_good(first.quality) {
                loop {
                    let next = last + 1;
                    let Some(base) = (next < columns).then(|| contiguous(last, next)).flatten()
                    else {
                        break;
                    };
                    if !options.grows(base.quality) || base.base == b'N' || reference[next] == b'N'
                    {
                        break;
                    }
                    if base.base != reference[next] {
                        alt.push(base.base);
                        quality += f64::from(base.quality);
                        appended_mismatches += 1;
                        last = next;
                        continue;
                    }
                    let bridge = (1..=options.mnv_gap)
                        .map(|gap| next + gap)
                        .take_while(|&to| to < columns)
                        .map_while(|to| contiguous(last, to).map(|base| (to, base)))
                        .find(|&(to, base)| base.base != reference[to]);
                    let Some((to, mismatch)) = bridge else {
                        break;
                    };
                    if !options.grows(mismatch.quality) {
                        break;
                    }
                    for base in (next..=to)
                        .map(|at_column| at(at_column).expect("a bridged base is contiguous"))
                    {
                        if base.base == b'N' {
                            dropped = true;
                        }
                        alt.push(base.base);
                        quality += f64::from(base.quality);
                    }
                    last = to;
                }
            }
            let near_indel = self.indels.iter().any(|&indel| {
                (indel < column && indel + options.mnv_gap >= column)
                    || (indel >= last && indel <= last + options.mnv_gap)
            });
            if !dropped && !near_indel {
                let in_read = at(last).expect("the last base is held").query_position
                    - self.leading_clip
                    - clips.head;
                let from_end = aligned_length.saturating_sub(in_read);
                found.push(Found {
                    start: window.start as u64 + column as u64 + 1,
                    quality: quality / alt.len() as f64,
                    position_in_read: (in_read + 1).min(from_end) as i32,
                    mismatches: self.mismatches - f64::from(appended_mismatches),
                    alt,
                });
            }
            column = last + 1;
        }
        found
    }
}

/// One allele a read feeds, with VarDict's per-read measures of it.
#[derive(Clone, Debug, PartialEq)]
struct Found {
    /// The allele's 1-based start.
    start: u64,
    /// The allele's bases as VarDict writes ALT, with any matching bases it bridged.
    alt: Vec<u8>,
    /// The mean quality of the allele's bases, VarDict's per-read variant quality.
    quality: f64,
    /// The 1-based distance of the allele's last base from the nearer end of the aligned read.
    position_in_read: i32,
    /// The read's substitution mismatches outside the allele, as VarDict counts them.
    mismatches: f64,
}

impl Found {
    /// Whether this is a row's allele.
    fn is(&self, row: &TumorOnlyVariant) -> bool {
        self.start == row.start && self.alt == row.alt_allele.as_bytes()
    }
}

/// A soft clip starting at a site: its bases with their qualities, outward from the alignment, and
/// the direction it runs in along the contig, 1 at a read's 3' end and -1 at its 5' end.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Clip {
    bases: Vec<(u8, u8)>,
    direction: i64,
}

/// How a soft clip starting at an SNV with the ALT base reads past it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ClipSupport {
    /// The rest matches the reference, so the read supports the call.
    Clean,
    /// The rest does not match the reference, though its second base does, so VarDict counted the
    /// read anyway.
    Foreign,
}

/// The most clipped bases VarDict's `adjSNV` reads into an SNV.
const ADJ_SNV_MAX_CLIP: usize = 5;
/// The base quality at or below which VarDict reads a clipped base as low.
const CLIP_LOW_QUALITY: u8 = 12;
/// The base quality VarDict's `CigarParser` needs to extend an alignment through a clipped base.
const CLIP_EXTENSION_QUALITY: u8 = 10;

/// Recounts an SNV's high-quality ALT reads in a surveyed sample: a read supports the call when its
/// walk, as VarDict clips it, feeds the SNV at or above `-q`, or when its soft clip at the site is
/// clean; a foreign clip is counted apart.
fn recount(
    survey: &Survey,
    window: &Window,
    call: &TumorOnlyVariant,
    options: ReadOptions,
    reference: &mut Reference,
) -> Result<Recount, String> {
    let mut recount = Recount::default();
    let Some(column) = window.column(call.start) else {
        return Ok(recount);
    };
    let alt = call.alt_allele.as_bytes()[0].to_ascii_uppercase();
    for read in &survey.reads {
        if read
            .clipped
            .iter()
            .any(|found| found.is(call) && found.quality >= options.base_quality)
        {
            recount.aligned += 1;
            continue;
        }
        let Some(clip) = read.clip_at(column, window, reference)? else {
            continue;
        };
        match classify_clip(
            &clip,
            alt,
            &window.contig,
            call.start as i64,
            options,
            reference,
        )? {
            Some(ClipSupport::Clean) => recount.clean_clips += 1,
            Some(ClipSupport::Foreign) => recount.foreign_clips += 1,
            None => {}
        }
    }
    Ok(recount)
}

/// How a soft clip at an SNV reads, as VarDict reads it: the clip is cut at its first no-call and
/// before its second base at a low quality, as `CigarParser` cuts it, and is read only when it
/// keeps a base and more bases than low ones; it is a high-quality ALT read only when it starts
/// with the ALT base and its mean quality reaches `-q`. It is clean when the rest matches the
/// reference by VarDict's `ismatchref`, with at most 3 mismatches and under 15%, and foreign when
/// the rest does not match but its first base does and the clip is short enough for `adjSNV`.
fn classify_clip(
    clip: &Clip,
    alt: u8,
    contig: &str,
    site: i64,
    options: ReadOptions,
    reference: &mut Reference,
) -> Result<Option<ClipSupport>, String> {
    let mut kept: Vec<(u8, u8)> = Vec::new();
    let mut low = 0;
    for &(base, quality) in &clip.bases {
        if base.eq_ignore_ascii_case(&b'N') {
            break;
        }
        if quality <= CLIP_LOW_QUALITY {
            low += 1;
            if low > 1 {
                break;
            }
        }
        kept.push((base.to_ascii_uppercase(), quality));
    }
    if kept.is_empty() || kept.len() <= low || kept[0].0 != alt {
        return Ok(None);
    }
    let mean = kept
        .iter()
        .map(|&(_, quality)| f64::from(quality))
        .sum::<f64>()
        / kept.len() as f64;
    if mean < options.base_quality {
        return Ok(None);
    }
    let rest: Vec<u8> = kept[1..].iter().map(|&(base, _)| base).collect();
    let outward = outward(
        reference,
        contig,
        site + clip.direction,
        clip.direction,
        rest.len(),
    )?;
    if matches_reference(&rest, &outward) {
        Ok(Some(ClipSupport::Clean))
    } else if kept.len() <= ADJ_SNV_MAX_CLIP && outward.first() == rest.first() {
        Ok(Some(ClipSupport::Foreign))
    } else {
        Ok(None)
    }
}

/// VarDict's `ismatchref`: whether clipped bases match the reference with at most 3 mismatches and
/// under 15% of them, and never where the reference runs out; nothing to compare matches.
fn matches_reference(bases: &[u8], reference: &[u8]) -> bool {
    if reference.len() < bases.len() {
        return false;
    }
    let mismatches = bases
        .iter()
        .zip(reference)
        .filter(|(base, reference)| base != reference)
        .count();
    mismatches <= 3 && (bases.is_empty() || (mismatches as f64) < 0.15 * bases.len() as f64)
}

/// How many of a clip's innermost bases VarDict's `CigarParser` aligns before reading the clip:
/// those matching the reference, from the 1-based position of the innermost outward, at a quality
/// above 10.
fn extension(
    clip: &[(u8, u8)],
    contig: &str,
    innermost: i64,
    direction: i64,
    reference: &mut Reference,
) -> Result<usize, String> {
    let outward = outward(reference, contig, innermost, direction, clip.len())?;
    Ok(clip
        .iter()
        .zip(&outward)
        .take_while(|((base, quality), reference)| {
            base.to_ascii_uppercase() == **reference && *quality > CLIP_EXTENSION_QUALITY
        })
        .count())
}

/// Whether two insertions, each as the bases inserted after a 1-based anchor position, give the
/// same haplotype: the same bases at one anchor, or, where the reference between the anchors
/// repeats the inserted bases, the bases rotated by the distance between them.
fn same_insertion(
    contig: &str,
    (anchor, bases): (u64, &[u8]),
    (other_anchor, other): (u64, &[u8]),
    reference: &mut Reference,
) -> Result<bool, String> {
    if bases.len() != other.len() || bases.is_empty() {
        return Ok(false);
    }
    let ((left, left_bases), (right, right_bases)) = if anchor <= other_anchor {
        ((anchor, bases), (other_anchor, other))
    } else {
        ((other_anchor, other), (anchor, bases))
    };
    if right > reference.length(contig)? {
        return Ok(false);
    }
    let between = reference.bases(contig, left, right)?;
    let repeats = between
        .iter()
        .enumerate()
        .all(|(offset, &base)| base == left_bases[offset % left_bases.len()]);
    let rotation = between.len() % left_bases.len();
    let rotated = left_bases[rotation..].iter().chain(&left_bases[..rotation]);
    Ok(repeats && rotated.eq(right_bases.iter()))
}

/// Up to `count` uppercase reference bases from a 1-based position outward along `direction`,
/// fewer where the contig ends.
fn outward(
    reference: &mut Reference,
    contig: &str,
    from: i64,
    direction: i64,
    count: usize,
) -> Result<Vec<u8>, String> {
    let length = reference.length(contig)? as i64;
    let (start, end) = if direction > 0 {
        (from - 1, (from - 1 + count as i64).min(length))
    } else {
        ((from - count as i64).max(0), from.min(length))
    };
    if start < 0 || start >= end {
        return Ok(Vec::new());
    }
    let mut bases = reference.bases(contig, start as u64, end as u64)?.to_vec();
    if direction < 0 {
        bases.reverse();
    }
    Ok(bases)
}

/// Whose molecule a template with the partial's ALT bases is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Molecule {
    /// The MNV's ALT at every other changed position, at a good quality.
    Mnv,
    /// A good-quality base that is not the MNV's ALT at some other changed position.
    Separate,
    /// No good-quality base at some other changed position.
    Ambiguous,
}

impl TemplateView {
    /// Whose molecule this template is, by its called bases at the MNV's other changed columns.
    /// A base under the pileup's floor never votes, so the quality check here only catches two
    /// disagreeing mates, whose called quality is the difference of theirs.
    fn molecule(&self, others: &[(usize, u8)], options: ReadOptions) -> Molecule {
        let mut ambiguous = false;
        for &(column, alt) in others {
            match self.called[column] {
                Some((base, quality)) if base != b'N' && options.is_good(quality) => {
                    if base != alt {
                        return Molecule::Separate;
                    }
                }
                _ => ambiguous = true,
            }
        }
        if ambiguous {
            Molecule::Ambiguous
        } else {
            Molecule::Mnv
        }
    }
}

/// The reads of one fate, summed as VarDict sums an allele's reads.
#[derive(Clone, Debug, Default, PartialEq)]
struct Reads {
    count: i32,
    forward: i32,
    high_quality: i32,
    positions: Vec<i32>,
    quality: f64,
    mapping_quality: f64,
    mismatches: f64,
}

impl Reads {
    /// Adds one read's allele.
    fn add(&mut self, read: &ReadView, found: &Found, options: ReadOptions) {
        self.count += 1;
        self.forward += i32::from(!read.reverse);
        self.high_quality += i32::from(found.quality >= options.base_quality);
        self.positions.push(found.position_in_read);
        self.quality += found.quality;
        self.mapping_quality += f64::from(read.mapping_quality);
        self.mismatches += found.mismatches;
    }

    /// The row's call holding exactly these reads, with their counts and means.
    fn call<'a>(&self, row: &TumorOnlyVariant<'a>) -> TumorOnlyVariant<'a> {
        if self.count == 0 {
            return without_reads(row);
        }
        let reads = f64::from(self.count);
        let positions_vary = self.positions.windows(2).any(|pair| pair[0] != pair[1]);
        TumorOnlyVariant {
            alt_depth: self.count,
            alt_forward: self.forward,
            alt_reverse: self.count - self.forward,
            high_quality_variant_reads: self.high_quality,
            mean_position_in_read: (self.positions.iter().map(|&p| f64::from(p)).sum::<f64>()
                / reads) as f32,
            stdev_position_in_read: if positions_vary { 1.0 } else { 0.0 },
            base_quality_mean: (self.quality / reads) as f32,
            mean_mapping_quality: (self.mapping_quality / reads) as f32,
            mean_mismatches_in_reads: (self.mismatches / reads) as f32,
            ..row.clone()
        }
    }
}

/// The row's call with none of its ALT reads.
fn without_reads<'a>(row: &TumorOnlyVariant<'a>) -> TumorOnlyVariant<'a> {
    TumorOnlyVariant {
        alt_depth: 0,
        alt_forward: 0,
        alt_reverse: 0,
        high_quality_variant_reads: 0,
        ..row.clone()
    }
}

/// How one partial's reads in one sample were classified.
#[derive(Clone, Debug, Default, PartialEq)]
struct Decision {
    /// Reads whose template is the MNV's molecule.
    proven: i32,
    /// Reads whose template is another allele's.
    separate: i32,
    /// Reads whose template proves nothing.
    ambiguous: i32,
    /// Reads VarDict counted that no read in the pileup feeds, or, below zero, reads in the
    /// pileup feeding the row beyond VarDict's count.
    unmatched: i32,
    /// The same, for the reads that still feed the row once VarDict's end clipping is replayed.
    unmatched_clipped: i32,
    /// The fraction of the templates at the MNV's other changed positions that prove nothing there.
    ambiguity: f64,
}

impl Decision {
    /// Whether the reads found number VarDict's, as aligned or as VarDict clips them.
    fn reconciled(&self) -> bool {
        self.unmatched == 0 || self.unmatched_clipped == 0
    }
}

/// One read feeding a partial, the allele it feeds and whose molecule its template is.
type Feeding<'s> = (&'s ReadView, Found, Molecule);

/// One sample's reads feeding a partial: as VarDict walks them aligned, and as it walks them after
/// clipping their ends. The decision tallies the walk that numbers VarDict's reads, else the aligned.
struct Classified<'s> {
    feeding: Vec<Feeding<'s>>,
    clipped: Vec<Feeding<'s>>,
    decision: Decision,
}

impl<'s> Classified<'s> {
    /// The walk whose reads number VarDict's: as aligned first, else as VarDict clips them.
    fn reconciled(&self) -> Option<&[Feeding<'s>]> {
        if self.decision.unmatched == 0 {
            Some(&self.feeding)
        } else if self.decision.unmatched_clipped == 0 {
            Some(&self.clipped)
        } else {
            None
        }
    }

    /// The reads the decision tallies.
    fn tallied(&self) -> &[Feeding<'s>] {
        self.reconciled().unwrap_or(&self.feeding)
    }
}

/// Classifies one sample's reads of a partial under an MNV by their templates' bases at the MNV's
/// changed positions outside the partial.
fn classify<'s>(
    survey: &'s Survey,
    window: &Window,
    mnv: &TumorOnlyVariant,
    partial: &TumorOnlyVariant,
    options: ReadOptions,
) -> Classified<'s> {
    let others: Vec<(usize, u8)> = mnv
        .ref_allele
        .bytes()
        .zip(mnv.alt_allele.bytes())
        .enumerate()
        .filter(|&(offset, (reference, alt))| {
            reference != alt && {
                let position = mnv.start + offset as u64;
                position < partial.start || position > partial.end
            }
        })
        .filter_map(|(offset, (_, alt))| {
            window
                .column(mnv.start + offset as u64)
                .map(|column| (column, alt.to_ascii_uppercase()))
        })
        .collect();
    let others = &others;
    let feed = |alleles: fn(&'s ReadView) -> &'s [Found]| -> Vec<Feeding<'s>> {
        survey
            .reads
            .iter()
            .flat_map(|read| {
                alleles(read)
                    .iter()
                    .filter(|allele| allele.is(partial))
                    .map(move |allele| {
                        let molecule = survey.templates[read.template].molecule(others, options);
                        (read, allele.clone(), molecule)
                    })
            })
            .collect()
    };
    let feeding = feed(|read| &read.aligned);
    let clipped = feed(|read| &read.clipped);
    let mut decision = Decision {
        unmatched: partial.alt_depth - feeding.len() as i32,
        unmatched_clipped: partial.alt_depth - clipped.len() as i32,
        ..Decision::default()
    };
    let mut classified = Classified {
        feeding,
        clipped,
        decision: decision.clone(),
    };
    for (_, _, molecule) in classified.tallied() {
        match molecule {
            Molecule::Mnv => decision.proven += 1,
            Molecule::Separate => decision.separate += 1,
            Molecule::Ambiguous => decision.ambiguous += 1,
        }
    }
    let (covering, ambiguous) = survey
        .templates
        .iter()
        .filter(|template| others.iter().any(|&(column, _)| template.present[column]))
        .fold((0, 0), |(covering, ambiguous), template| {
            let proves_nothing = template.molecule(others, options) == Molecule::Ambiguous;
            (covering + 1, ambiguous + i64::from(proves_nothing))
        });
    decision.ambiguity = if covering > 0 {
        ambiguous as f64 / covering as f64
    } else {
        1.0
    };
    classified.decision = decision;
    classified
}

/// What one sample's reads of a partial would do, and did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct Outcome {
    /// Reads that moved to the MNV.
    moved: i32,
    /// Reads that stayed.
    kept: i32,
    /// Reads found that would move by their class, whether or not the row reconciled.
    would_move: i32,
}

/// Decides a partial's reads in every sample at once, the tumor (or only) sample first, from each
/// sample's classified reads: the fate of the row, and each sample's outcome.
///
/// A read whose template is the MNV's molecule moves whatever the MNV's count in the sample. A read
/// whose template is another allele's stays. A template that proves nothing is as likely the MNV's as
/// another allele's as those are common among the molecules with the partial's bases: the MNV's
/// reads, and the row's templates shown to be the MNV's or another allele's, over the samples whose
/// reads numbered VarDict's, since an allele in one sample is in the individual. Of a sample's `u`
/// unknown templates, `u s / (m + p + s)` rounded stay, with `m` the MNV's reads, `p` the proven and
/// `s` the separate templates, and the rest move: a lone separate template beside a deep MNV keeps
/// none, while one beside a shallow MNV keeps them all. The fraction of all templates that prove
/// nothing at those positions bounds how often an allele's `k` templates all prove nothing: that
/// fraction to the power `k`. Unknown templates move only when the MNV has reads in their sample,
/// since nothing leaks from an MNV with none; and, in a sample after the tumor, only when the tumor
/// keeps no read of the row, so a normal never gives up the reads that say it carries an allele the
/// tumor keeps.
///
/// A sample's reads are split between the MNV and the row only when the reads found are exactly the
/// reads VarDict counted, as aligned or else as VarDict clips them, so every count and mean written
/// comes from a read seen. When the pileup feeds the row more reads than VarDict counted, the row
/// moves whole, with VarDict's own values, only if every read found as aligned would move;
/// otherwise, and whenever it feeds fewer, the sample keeps VarDict's row whole, since which of its
/// reads VarDict counted is not known.
fn resolve<'a>(
    classified: &[Classified],
    mnv: &[TumorOnlyVariant<'a>],
    partials: &[&TumorOnlyVariant<'a>],
    options: ReadOptions,
) -> (ReadFate<'a>, Vec<Outcome>) {
    let pooled = classified
        .iter()
        .filter(|sample| sample.reconciled().is_some())
        .map(|sample| templates(sample.tallied()))
        .fold((0, 0), |(proven, separate), (p, s, _)| {
            (proven + p, separate + s)
        });
    let molecules = mnv
        .iter()
        .map(|call| call.alt_depth.max(0) as usize)
        .sum::<usize>()
        + pooled.0
        + pooled.1;
    let mut fate = ReadFate {
        moved: Vec::new(),
        kept: Vec::new(),
    };
    let mut outcomes = Vec::new();
    let mut tumor_keeps = false;
    for (sample, (classified, partial)) in classified.iter().zip(partials).enumerate() {
        let unknown = templates(classified.tallied()).2;
        let keep_unknown = if mnv[sample].alt_depth > 0 && (sample == 0 || !tumor_keeps) {
            Some((2 * unknown * pooled.1 + molecules) / (2 * molecules.max(1)))
        } else {
            None
        };
        let (moved, kept, would_move) = settle(classified, partial, keep_unknown, options);
        if sample == 0 {
            tumor_keeps = kept.alt_depth > 0;
        }
        outcomes.push(Outcome {
            moved: moved.alt_depth,
            kept: kept.alt_depth,
            would_move,
        });
        fate.moved.push(moved);
        fate.kept.push(kept);
    }
    (fate, outcomes)
}

/// The templates among these reads of each class: (proven, separate, unknown), each counted once.
fn templates(feeding: &[Feeding]) -> (usize, usize, usize) {
    let mut seen: [HashSet<usize>; 3] = Default::default();
    for (read, _, molecule) in feeding {
        let class = match molecule {
            Molecule::Mnv => 0,
            Molecule::Separate => 1,
            Molecule::Ambiguous => 2,
        };
        seen[class].insert(read.template);
    }
    (seen[0].len(), seen[1].len(), seen[2].len())
}

/// One sample's calls of a partial's reads that move and that stay, as [`resolve`] decides them,
/// and the reads found that would move by their class; `keep_unknown` is how many unknown templates
/// stay, the first in the pileup, or every one when `None`.
fn settle<'a>(
    classified: &Classified,
    partial: &TumorOnlyVariant<'a>,
    keep_unknown: Option<usize>,
    options: ReadOptions,
) -> (TumorOnlyVariant<'a>, TumorOnlyVariant<'a>, i32) {
    let mut staying: HashSet<usize> = HashSet::new();
    for (read, _, molecule) in classified.tallied() {
        if *molecule == Molecule::Ambiguous && keep_unknown.is_none_or(|keep| staying.len() < keep)
        {
            staying.insert(read.template);
        }
    }
    let goes = |read: &ReadView, molecule: Molecule| match molecule {
        Molecule::Mnv => true,
        Molecule::Ambiguous => !staying.contains(&read.template),
        Molecule::Separate => false,
    };
    let movable = |feeding: &[Feeding]| -> i32 {
        feeding
            .iter()
            .filter(|(read, _, molecule)| goes(read, *molecule))
            .count() as i32
    };
    let would_move = movable(classified.tallied());
    if partial.alt_depth > 0 {
        if let Some(feeding) = classified.reconciled() {
            let mut moved = Reads::default();
            let mut kept = Reads::default();
            for (read, allele, molecule) in feeding {
                if goes(read, *molecule) {
                    moved.add(read, allele, options);
                } else {
                    kept.add(read, allele, options);
                }
            }
            return (moved.call(partial), kept.call(partial), would_move);
        }
        let found = classified.feeding.len() as i32;
        if found > partial.alt_depth && would_move == found {
            return (partial.clone(), without_reads(partial), would_move);
        }
    }
    (without_reads(partial), partial.clone(), would_move)
}

/// Parses VarDict's `-F` flag value: hexadecimal with a `0x` prefix, or decimal.
pub fn parse_flags(text: &str) -> Result<u16, String> {
    let parsed = match text.strip_prefix("0x").or_else(|| text.strip_prefix("0X")) {
        Some(hex) => u16::from_str_radix(hex, 16),
        None => text.parse(),
    };
    parsed.map_err(|_| {
        format!("Flags must be a hexadecimal value like 0x504 or a decimal one, not {text}!")
    })
}

/// Splits `tumor.bam|normal.bam` into its paths, as VarDict's `-b` takes them.
pub fn split_bams(text: &str) -> Vec<PathBuf> {
    text.split('|')
        .map(Path::new)
        .map(Path::to_path_buf)
        .collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use streampile::testing::{Frag, Pair, SamBuilder, Strand};

    use super::*;

    /// A window of a reference over 1-based positions 101 to 120 of `chr1`.
    const REFERENCE: &[u8] = b"ACGTACGTACGTACGTACGT";

    fn window() -> Window {
        Window {
            contig: "chr1".to_string(),
            start: 100,
            reference: REFERENCE.to_vec(),
        }
    }

    /// Read bases over the window with the given substitutions at 1-based positions.
    fn edited(edits: &[(u64, u8)]) -> String {
        let mut bases = REFERENCE.to_vec();
        for &(position, base) in edits {
            bases[(position - 101) as usize] = base;
        }
        String::from_utf8(bases).unwrap()
    }

    fn quals(edits: &[(u64, u8)]) -> Vec<u8> {
        let mut quals = vec![40u8; REFERENCE.len()];
        for &(position, quality) in edits {
            quals[(position - 101) as usize] = quality;
        }
        quals
    }

    /// A VarDict row in the window with these counts in one sample.
    fn row(
        start: u64,
        ref_allele: &'static str,
        alt_allele: &'static str,
        alt_depth: i32,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            contig: "chr1",
            start,
            end: start + ref_allele.len() as u64 - 1,
            ref_allele,
            alt_allele,
            depth: 100,
            alt_depth,
            alt_forward: alt_depth,
            alt_reverse: 0,
            high_quality_variant_reads: alt_depth,
            mean_position_in_read: 7.0,
            base_quality_mean: 33.0,
            mean_mapping_quality: 50.0,
            mean_mismatches_in_reads: 1.5,
            ..Default::default()
        }
    }

    /// Surveys fragments with these bases and qualities, one per read, over the window.
    fn survey_frags(reads: &[(&str, Vec<u8>)]) -> Survey {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        for (bases, quals) in reads {
            builder.add_frag(Frag::at(101).bases(*bases).quals(quals.clone()));
        }
        survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap()
    }

    fn alleles(survey: &Survey) -> Vec<Vec<(u64, String)>> {
        survey
            .reads
            .iter()
            .map(|read| {
                read.aligned
                    .clone()
                    .into_iter()
                    .map(|found| (found.start, String::from_utf8(found.alt).unwrap()))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn test_a_read_with_contiguous_mismatches_feeds_one_mnv() {
        let survey = survey_frags(&[(&edited(&[(105, b'T'), (106, b'A')]), quals(&[]))]);
        assert_eq!(alleles(&survey), vec![vec![(105, "TA".to_string())]]);
    }

    #[test]
    fn test_mismatches_within_the_gap_are_bridged_into_one_mnv() {
        let survey = survey_frags(&[(&edited(&[(105, b'T'), (108, b'A')]), quals(&[]))]);
        assert_eq!(alleles(&survey), vec![vec![(105, "TCGA".to_string())]]);
    }

    #[test]
    fn test_mismatches_beyond_the_gap_are_two_alleles() {
        let survey = survey_frags(&[(&edited(&[(105, b'T'), (109, b'G')]), quals(&[]))]);
        assert_eq!(
            alleles(&survey),
            vec![vec![(105, "T".to_string()), (109, "G".to_string())]]
        );
    }

    #[rstest]
    #[case::below_plus_five(27, vec![(105, "T"), (106, "A")])]
    #[case::at_plus_five(28, vec![(105, "TA")])]
    fn test_an_mnv_grows_only_through_a_base_at_q_plus_five(
        #[case] quality: u8,
        #[case] expected: Vec<(u64, &str)>,
    ) {
        let survey = survey_frags(&[(
            &edited(&[(105, b'T'), (106, b'A')]),
            quals(&[(106, quality)]),
        )]);
        let expected: Vec<(u64, String)> = expected
            .into_iter()
            .map(|(s, a)| (s, a.to_string()))
            .collect();
        assert_eq!(alleles(&survey), vec![expected]);
    }

    #[test]
    fn test_a_first_base_below_q_does_not_grow() {
        let survey = survey_frags(&[(&edited(&[(105, b'T'), (106, b'A')]), quals(&[(105, 20)]))]);
        assert_eq!(
            alleles(&survey),
            vec![vec![(105, "T".to_string()), (106, "A".to_string())]]
        );
    }

    #[test]
    fn test_a_no_call_stops_the_mnv_and_feeds_no_allele() {
        let survey = survey_frags(&[(
            &edited(&[(105, b'T'), (106, b'N'), (107, b'A')]),
            quals(&[]),
        )]);
        assert_eq!(
            alleles(&survey),
            vec![vec![(105, "T".to_string()), (107, "A".to_string())]]
        );
    }

    #[test]
    fn test_a_no_call_bridged_into_an_mnv_drops_the_allele() {
        let survey = survey_frags(&[(&edited(&[(105, b'T'), (107, b'N')]), quals(&[]))]);
        assert_eq!(alleles(&survey), vec![Vec::<(u64, String)>::new()]);
    }

    #[test]
    fn test_a_read_ending_inside_the_mnv_feeds_the_prefix() {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        let bases = edited(&[(105, b'T'), (106, b'A')]);
        builder.add_frag(
            Frag::at(101)
                .bases(&bases[..5])
                .quals(vec![40; 5])
                .cigar("5M"),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(alleles(&survey), vec![vec![(105, "T".to_string())]]);
    }

    #[test]
    fn test_a_read_with_an_indel_past_the_gap_feeds_its_allele() {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        let bases = edited(&[(105, b'T'), (106, b'A')]);
        builder.add_frag(
            Frag::at(101)
                .bases(&bases[..19])
                .quals(vec![40; 19])
                .cigar("10M1D9M"),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(survey.reads[0].indels, vec![10]);
        assert_eq!(alleles(&survey), vec![vec![(105, "TA".to_string())]]);
    }

    #[test]
    fn test_the_per_read_measures_are_vardicts() {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        let bases = edited(&[(105, b'T'), (106, b'A'), (115, b'C')]);
        builder.add_frag(
            Frag::at(103)
                .bases(&bases[..18])
                .quals(vec![40; 18])
                .cigar("2S16M")
                .attr(Tag::EDIT_DISTANCE, 3),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        let found = survey.reads[0].aligned.clone();
        assert_eq!(found[0].alt, b"TA");
        assert_eq!(found[0].quality, 40.0);
        // The MNV's last base is the 4th aligned base of 16, so 4 from the start and 12 from the end.
        assert_eq!(found[0].position_in_read, 4);
        // NM 3 minus the one mismatch VarDict appends after the first.
        assert_eq!(found[0].mismatches, 2.0);
    }

    #[test]
    fn test_a_read_over_the_mismatch_limit_is_left_out() {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        let bases = edited(&[(105, b'T')]);
        builder.add_frag(
            Frag::at(101)
                .bases(&bases)
                .quals(vec![40; 20])
                .attr(Tag::EDIT_DISTANCE, 9),
        );
        builder.add_frag(
            Frag::at(101)
                .bases(&bases)
                .quals(vec![40; 20])
                .attr(Tag::EDIT_DISTANCE, 8),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(survey.reads.len(), 1);
    }

    #[test]
    fn test_a_read_under_the_mapping_quality_floor_is_left_out() {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        builder.add_frag(
            Frag::at(101)
                .bases(edited(&[]))
                .quals(vec![40; 20])
                .mapq(19),
        );
        let options = ReadOptions {
            min_mapping_quality: 20,
            ..ReadOptions::default()
        };
        let survey = survey(builder.to_source(), builder.header(), &window(), options).unwrap();
        assert!(survey.reads.is_empty());
    }

    /// A pair of mates, each given as (bases, quals).
    type Mates<'a> = ((&'a str, Vec<u8>), (&'a str, Vec<u8>));

    /// Builds paired reads over the window, mate 2 reverse, named `q0`, `q1`, ...
    fn pairs_of(pairs: &[Mates]) -> SamBuilder {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        for (index, (read1, read2)) in pairs.iter().enumerate() {
            builder.add_pair(
                Pair::at(101, 101)
                    .name(format!("q{index}"))
                    .bases1(read1.0)
                    .quals1(read1.1.clone())
                    .bases2(read2.0)
                    .quals2(read2.1.clone())
                    .strand2(Strand::Minus),
            );
        }
        builder
    }

    /// Decides one partial under one MNV in one sample from a survey.
    fn decide_survey(
        survey: &Survey,
        mnv: &TumorOnlyVariant<'static>,
        partial: &TumorOnlyVariant<'static>,
    ) -> (
        TumorOnlyVariant<'static>,
        TumorOnlyVariant<'static>,
        Decision,
    ) {
        let options = ReadOptions::default();
        let classified = classify(survey, &window(), mnv, partial, options);
        let (fate, _) = resolve(
            std::slice::from_ref(&classified),
            std::slice::from_ref(mnv),
            &[partial],
            options,
        );
        (
            fate.moved[0].clone(),
            fate.kept[0].clone(),
            classified.decision,
        )
    }

    /// Decides one partial under one MNV from paired reads, with whether the reads reconciled.
    fn decide_pairs(
        pairs: &[Mates],
        mnv: &TumorOnlyVariant<'static>,
        partial: &TumorOnlyVariant<'static>,
    ) -> (TumorOnlyVariant<'static>, TumorOnlyVariant<'static>, bool) {
        let (moved, kept, decision) = decide_pairs_fully(pairs, mnv, partial);
        (moved, kept, decision.unmatched == 0)
    }

    fn decide_pairs_fully(
        pairs: &[Mates],
        mnv: &TumorOnlyVariant<'static>,
        partial: &TumorOnlyVariant<'static>,
    ) -> (
        TumorOnlyVariant<'static>,
        TumorOnlyVariant<'static>,
        Decision,
    ) {
        let builder = pairs_of(pairs);
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        decide_survey(&survey, mnv, partial)
    }

    #[test]
    fn test_a_mate_with_the_alt_at_the_other_changed_base_proves_the_mnv_molecule() {
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 1);
        // Read 1 has the second changed base at a quality too low to grow; read 2 holds it well.
        let read1 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let read2 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let (moved, kept, reconciled) = decide_pairs(
            &[((&read1.0, read1.1), (&read2.0, read2.1))],
            &mnv,
            &partial,
        );
        assert_eq!((moved.alt_depth, kept.alt_depth, reconciled), (1, 0, true));
        assert_eq!((moved.alt_forward, moved.alt_reverse), (1, 0));
    }

    #[test]
    fn test_a_mate_with_the_reference_at_the_other_changed_base_keeps_the_read() {
        // Both unclipped mates stop at the first base, so VarDict counted two reads of the partial.
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 2);
        let read1 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let read2 = (edited(&[(105, b'T')]), quals(&[]));
        let (moved, kept, reconciled) = decide_pairs(
            &[((&read1.0, read1.1), (&read2.0, read2.1))],
            &mnv,
            &partial,
        );
        assert_eq!((moved.alt_depth, kept.alt_depth, reconciled), (0, 2, true));
        assert_eq!(kept.base_quality_mean, 40.0);
        assert_eq!((kept.alt_forward, kept.alt_reverse), (1, 1));
    }

    #[test]
    fn test_a_no_call_at_the_other_changed_base_is_leak_only_under_an_mnv_with_reads() {
        let partial = row(105, "A", "T", 1);
        let read1 = (edited(&[(105, b'T'), (106, b'N')]), quals(&[]));
        let read2 = (edited(&[(105, b'T'), (106, b'N')]), quals(&[]));
        let pairs = [(
            (&read1.0[..], read1.1.clone()),
            (&read2.0[..], read2.1.clone()),
        )];
        let (moved, kept, _) = decide_pairs(&pairs, &row(105, "AC", "TA", 10), &partial);
        assert_eq!((moved.alt_depth, kept.alt_depth), (1, 0));
        let (moved, kept, _) = decide_pairs(&pairs, &row(105, "AC", "TA", 0), &partial);
        assert_eq!((moved.alt_depth, kept.alt_depth), (0, 1));
    }

    #[test]
    fn test_a_read_that_grew_through_the_mnv_feeds_the_mnv_not_the_partial() {
        let mnv = row(105, "AC", "TA", 1);
        let partial = row(105, "A", "T", 1);
        let full = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let (moved, kept, reconciled) = decide_pairs(
            &[((&full.0, full.1.clone()), (&full.0, full.1.clone()))],
            &mnv,
            &partial,
        );
        assert_eq!((moved.alt_depth, kept.alt_depth, reconciled), (0, 1, false));
    }

    #[test]
    fn test_the_alt_at_one_changed_base_and_the_reference_at_another_keeps_the_read() {
        // The template holds the MNV's second base in its mate but the reference at its last.
        let mnv = row(105, "ACGT", "TAGA", 10);
        let partial = row(105, "A", "T", 1);
        let read1 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let read2 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let (moved, kept, reconciled) = decide_pairs(
            &[((&read1.0, read1.1), (&read2.0, read2.1))],
            &mnv,
            &partial,
        );
        assert_eq!((moved.alt_depth, kept.alt_depth, reconciled), (0, 1, true));
    }

    /// Pairs of a partial's reads beside `others` templates that hold the reference, `masked` of
    /// them with a no-call at the MNV's other changed base.
    fn beside<'a>(partial_pairs: &[Mates<'a>], others: usize, masked: usize) -> Vec<Mates<'a>> {
        let reference: &'static str = Box::leak(edited(&[]).into_boxed_str());
        let with_no_call: &'static str = Box::leak(edited(&[(106, b'N')]).into_boxed_str());
        let mut pairs: Vec<Mates<'a>> = partial_pairs.to_vec();
        for index in 0..others {
            let bases = if index < masked {
                with_no_call
            } else {
                reference
            };
            pairs.push(((bases, quals(&[])), (bases, quals(&[]))));
        }
        pairs
    }

    #[test]
    fn test_ambiguous_reads_move_beside_a_lone_separate_template_under_a_deep_mnv() {
        let mnv = row(105, "AC", "TA", 10);
        let ambiguous = (edited(&[(105, b'T'), (106, b'N')]), quals(&[]));
        let separate = (edited(&[(105, b'T')]), quals(&[]));
        let own = [
            (
                (&ambiguous.0[..], ambiguous.1.clone()),
                (&ambiguous.0[..], ambiguous.1.clone()),
            ),
            (
                (&ambiguous.0[..], ambiguous.1.clone()),
                (&ambiguous.0[..], ambiguous.1.clone()),
            ),
            (
                (&separate.0[..], separate.1.clone()),
                (&separate.0[..], separate.1.clone()),
            ),
        ];
        let pairs = beside(&own, 20, 4);
        let (moved, kept, decision) = decide_pairs_fully(&pairs, &mnv, &row(105, "A", "T", 6));
        assert_eq!((decision.separate, decision.ambiguous), (2, 4));
        // The partial's two ambiguous templates and the four masked ones, of 23 at the other base.
        assert!((decision.ambiguity - 6.0 / 23.0).abs() < 1e-9);
        // Of 2 unknown templates, 2 * 1 / (10 + 1) rounds to none staying.
        assert_eq!((moved.alt_depth, kept.alt_depth), (4, 2));
    }

    #[test]
    fn test_a_proven_mnv_molecule_moves_beside_a_read_shown_to_be_another_alleles() {
        let mnv = row(105, "AC", "TA", 10);
        let proven = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let proven_mate = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let separate = (edited(&[(105, b'T')]), quals(&[]));
        let pairs = [
            (
                (&proven.0[..], proven.1.clone()),
                (&proven_mate.0[..], proven_mate.1.clone()),
            ),
            (
                (&separate.0[..], separate.1.clone()),
                (&separate.0[..], separate.1.clone()),
            ),
        ];
        let (moved, kept, decision) = decide_pairs_fully(&pairs, &mnv, &row(105, "A", "T", 3));
        assert_eq!(
            (decision.proven, decision.separate, decision.unmatched),
            (1, 2, 0)
        );
        assert_eq!((moved.alt_depth, kept.alt_depth), (1, 2));
    }

    #[test]
    fn test_a_row_with_reads_no_template_feeds_keeps_vardicts_row_whole() {
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 3);
        let read1 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let read2 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let (moved, kept, reconciled) = decide_pairs(
            &[((&read1.0, read1.1), (&read2.0, read2.1))],
            &mnv,
            &partial,
        );
        assert_eq!((moved.alt_depth, kept.alt_depth, reconciled), (0, 3, false));
        assert_eq!(kept, partial);
    }

    #[test]
    fn test_a_row_fed_more_reads_than_vardict_counted_keeps_them_all_when_one_is_another_alleles() {
        // Two templates prove the MNV molecule and one carries the reference, but VarDict counted
        // only two reads, so which two is not known and nothing moves.
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 2);
        let proven = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let mate = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let separate = (edited(&[(105, b'T')]), quals(&[]));
        let pairs = [
            (
                (&proven.0[..], proven.1.clone()),
                (&mate.0[..], mate.1.clone()),
            ),
            (
                (&proven.0[..], proven.1.clone()),
                (&mate.0[..], mate.1.clone()),
            ),
            (
                (&separate.0[..], separate.1.clone()),
                (&separate.0[..], separate.1.clone()),
            ),
        ];
        let (moved, kept, decision) = decide_pairs_fully(&pairs, &mnv, &partial);
        assert_eq!(
            (decision.proven, decision.separate, decision.unmatched),
            (2, 2, -2)
        );
        assert_eq!((moved.alt_depth, kept.alt_depth), (0, 2));
        assert_eq!(kept, partial);
    }

    #[test]
    fn test_a_row_fed_more_reads_than_vardict_counted_moves_whole_when_every_read_would() {
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 2);
        let ambiguous = (edited(&[(105, b'T'), (106, b'N')]), quals(&[]));
        let pairs = [
            (
                (&ambiguous.0[..], ambiguous.1.clone()),
                (&ambiguous.0[..], ambiguous.1.clone()),
            ),
            (
                (&ambiguous.0[..], ambiguous.1.clone()),
                (&ambiguous.0[..], ambiguous.1.clone()),
            ),
        ];
        let (moved, kept, decision) = decide_pairs_fully(&pairs, &mnv, &partial);
        assert_eq!((decision.ambiguous, decision.unmatched), (4, -2));
        // VarDict's own row moves, not a count made from the reads found.
        assert_eq!(moved, partial);
        assert_eq!(kept.alt_depth, 0);
    }

    #[test]
    fn test_a_kept_read_below_q_is_not_a_high_quality_read() {
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 2);
        let low = (edited(&[(105, b'T')]), quals(&[(105, 20)]));
        let pairs = [((&low.0[..], low.1.clone()), (&low.0[..], low.1.clone()))];
        let (moved, kept, _) = decide_pairs_fully(&pairs, &mnv, &partial);
        assert_eq!((moved.alt_depth, kept.alt_depth), (0, 2));
        assert_eq!(kept.high_quality_variant_reads, 0);
        assert_eq!(kept.base_quality_mean, 20.0);
    }

    #[test]
    fn test_disagreeing_mates_are_called_to_the_higher_quality_base_at_the_difference() {
        // At the other changed base, ALT at 90 against the reference at 40 calls ALT at 50, a good
        // call, and the reference at 90 against ALT at 40 calls the reference at 50.
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 4);
        let alt_90 = (
            edited(&[(105, b'T'), (106, b'A')]),
            quals(&[(105, 20), (106, 90)]),
        );
        let alt_40 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(105, 20)]));
        let ref_90 = (edited(&[(105, b'T')]), quals(&[(106, 90)]));
        let ref_40 = (edited(&[(105, b'T')]), quals(&[]));
        let pairs = [
            (
                (&ref_40.0[..], ref_40.1.clone()),
                (&alt_90.0[..], alt_90.1.clone()),
            ),
            (
                (&ref_90.0[..], ref_90.1.clone()),
                (&alt_40.0[..], alt_40.1.clone()),
            ),
        ];
        let (moved, kept, decision) = decide_pairs_fully(&pairs, &mnv, &partial);
        assert_eq!(
            (decision.proven, decision.separate, decision.ambiguous),
            (2, 2, 0)
        );
        assert_eq!((moved.alt_depth, kept.alt_depth), (2, 2));
    }

    #[test]
    fn test_the_ambiguity_rate_counts_only_templates_at_the_other_changed_base() {
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 1);
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        builder.add_frag(
            Frag::at(101)
                .bases(edited(&[(105, b'T'), (106, b'N')]))
                .quals(quals(&[])),
        );
        for _ in 0..3 {
            builder.add_frag(Frag::at(101).bases(edited(&[])).quals(quals(&[])));
        }
        // Reads ending before the other changed base are not at it, so they do not count.
        for _ in 0..4 {
            builder.add_frag(
                Frag::at(101)
                    .bases(&edited(&[])[..5])
                    .quals(vec![40; 5])
                    .cigar("5M"),
            );
        }
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        let (_, _, decision) = decide_survey(&survey, &mnv, &partial);
        assert!((decision.ambiguity - 0.25).abs() < 1e-9);
    }

    #[test]
    fn test_an_unchanged_base_inside_the_mnv_is_not_read() {
        // The MNV ACG>TCA keeps its middle base, so the reference there proves nothing against it.
        let mnv = row(105, "ACG", "TCA", 10);
        let partial = row(105, "A", "T", 1);
        let read1 = (edited(&[(105, b'T'), (107, b'A')]), quals(&[(107, 20)]));
        let read2 = (edited(&[(105, b'T'), (107, b'A')]), quals(&[]));
        let (moved, kept, decision) = decide_pairs_fully(
            &[(
                (&read1.0[..], read1.1.clone()),
                (&read2.0[..], read2.1.clone()),
            )],
            &mnv,
            &partial,
        );
        assert_eq!(decision.proven, 1);
        assert_eq!((moved.alt_depth, kept.alt_depth), (1, 0));
    }

    #[test]
    fn test_a_shorter_mnv_under_an_mnv_is_decided_by_the_bases_outside_it() {
        let mnv = row(105, "ACGT", "TAGA", 10);
        let sub_mnv = row(105, "AC", "TA", 1);
        // Read 1 grows through the first two bases and stops at the last, which its mate holds.
        let read1 = (
            edited(&[(105, b'T'), (106, b'A'), (108, b'A')]),
            quals(&[(108, 20)]),
        );
        let read2 = (edited(&[(105, b'T'), (106, b'A'), (108, b'A')]), quals(&[]));
        let (moved, kept, decision) = decide_pairs_fully(
            &[(
                (&read1.0[..], read1.1.clone()),
                (&read2.0[..], read2.1.clone()),
            )],
            &mnv,
            &sub_mnv,
        );
        assert_eq!((decision.proven, decision.unmatched), (1, 0));
        assert_eq!((moved.alt_depth, kept.alt_depth), (1, 0));
    }

    /// Surveys fragments at 101 with these flag bits added, under these options.
    fn survey_flagged(bits: &[u16], options: ReadOptions) -> Survey {
        let mut built = SamBuilder::new().read_length(REFERENCE.len());
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        for (index, &bits) in bits.iter().enumerate() {
            let record = built
                .add_frag(
                    Frag::at(101)
                        .name(format!("r{index}"))
                        .bases(edited(&[(105, b'T')]))
                        .quals(quals(&[])),
                )
                .remove(0);
            builder.extend([SamBuilder::with_flags(record, bits)]);
        }
        survey(builder.to_source(), builder.header(), &window(), options).unwrap()
    }

    #[rstest]
    #[case::duplicate(0x400, 0x504, 0)]
    #[case::qc_fail_counted_by_default(0x200, 0x504, 1)]
    #[case::qc_fail_excluded(0x200, 0xF00, 0)]
    #[case::supplementary_skipped_with_any_filter(0x800, 0x504, 0)]
    #[case::supplementary_counted_without_a_filter(0x800, 0, 1)]
    fn test_reads_are_left_out_as_vardict_leaves_them_out(
        #[case] bits: u16,
        #[case] exclude_flags: u16,
        #[case] expected: usize,
    ) {
        let options = ReadOptions {
            exclude_flags,
            ..ReadOptions::default()
        };
        assert_eq!(survey_flagged(&[bits], options).reads.len(), expected);
    }

    #[test]
    fn test_two_records_of_one_name_are_two_reads() {
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        builder.add_frag(
            Frag::at(101)
                .name("r")
                .bases(edited(&[(105, b'T')]))
                .quals(quals(&[])),
        );
        builder.add_frag(
            Frag::at(102)
                .name("r")
                .bases(&edited(&[(105, b'T')])[1..])
                .quals(vec![40; 19])
                .cigar("19M"),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(survey.reads.len(), 2);
        assert_eq!(survey.templates.len(), 1);
    }

    #[rstest]
    #[case::adjacent(5, 0)]
    #[case::within_the_gap(6, 0)]
    #[case::past_the_gap(7, 1)]
    fn test_an_indel_within_the_gap_of_an_allele_makes_it_complex(
        #[case] deleted: usize,
        #[case] expected: usize,
    ) {
        // The read carries T at 105 (column 4) and lacks the reference base at `deleted`.
        let mut bases = edited(&[(105, b'T')]).into_bytes();
        bases.remove(deleted);
        bases.push(b'A');
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        builder.add_frag(
            Frag::at(101)
                .bases(String::from_utf8(bases).unwrap())
                .quals(quals(&[]))
                .cigar(format!("{deleted}M1D{}M", 20 - deleted)),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(survey.reads[0].indels, vec![deleted]);
        let found = survey.reads[0].aligned.clone();
        assert_eq!(found.len(), expected, "{found:?}");
    }

    /// Resolves one partial across samples from each sample's classes, given as
    /// (proven, separate, ambiguous) reads, under MNVs with these reads.
    fn resolve_classes(samples: &[((i32, i32, i32), i32)]) -> Vec<(i32, i32)> {
        let options = ReadOptions::default();
        let mut surveys = Vec::new();
        let mut partials = Vec::new();
        let mut mnvs = Vec::new();
        for &((proven, separate, ambiguous), mnv_reads) in samples {
            let proven_read = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
            let proven_mate = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
            let separate_read = (edited(&[(105, b'T')]), quals(&[]));
            let ambiguous_read = (edited(&[(105, b'T'), (106, b'N')]), quals(&[]));
            let mut pairs: Vec<Mates> = Vec::new();
            for _ in 0..proven {
                pairs.push((
                    (&proven_read.0, proven_read.1.clone()),
                    (&proven_mate.0, proven_mate.1.clone()),
                ));
            }
            for _ in 0..separate {
                pairs.push((
                    (&separate_read.0, separate_read.1.clone()),
                    (&separate_read.0, separate_read.1.clone()),
                ));
            }
            for _ in 0..ambiguous {
                pairs.push((
                    (&ambiguous_read.0, ambiguous_read.1.clone()),
                    (&ambiguous_read.0, ambiguous_read.1.clone()),
                ));
            }
            let builder = pairs_of(&pairs);
            surveys
                .push(survey(builder.to_source(), builder.header(), &window(), options).unwrap());
            partials.push(row(105, "A", "T", proven + 2 * separate + 2 * ambiguous));
            mnvs.push(row(105, "AC", "TA", mnv_reads));
        }
        let partial_refs: Vec<&TumorOnlyVariant> = partials.iter().collect();
        let classified: Vec<Classified> = surveys
            .iter()
            .zip(&partial_refs)
            .map(|(survey, partial)| classify(survey, &window(), &mnvs[0], partial, options))
            .collect();
        resolve(&classified, &mnvs, &partial_refs, options)
            .1
            .into_iter()
            .map(|outcome| (outcome.moved, outcome.kept))
            .collect()
    }

    #[test]
    fn test_a_normal_keeps_its_ambiguous_reads_while_the_tumor_keeps_reads_of_the_row() {
        // The tumor's unknown template moves (1 * 1 / 21 rounds to 0) but its separate one stays.
        assert_eq!(
            resolve_classes(&[((0, 1, 1), 10), ((0, 0, 1), 10)]),
            vec![(2, 2), (0, 2)]
        );
    }

    #[rstest]
    #[case::deep_mnv_moves_it(((0, 1, 1), 10), vec![(2, 2)])]
    #[case::shallow_mnv_keeps_it(((0, 1, 1), 1), vec![(0, 4)])]
    #[case::half_rounds_up_and_keeps_it(((0, 1, 1), 1), vec![(0, 4)])]
    #[case::under_half_moves_it(((0, 1, 1), 2), vec![(2, 2)])]
    #[case::two_separate_templates(((0, 2, 1), 10), vec![(2, 4)])]
    #[case::three_unknown_templates_all_move(((0, 1, 3), 10), vec![(6, 2)])]
    #[case::six_unknown_templates_keep_one(((0, 1, 6), 10), vec![(10, 4)])]
    #[case::thirty_unknown_templates_keep_three(((0, 1, 30), 10), vec![(54, 8)])]
    #[case::proven_templates_count_for_the_mnv(((5, 1, 1), 1), vec![(7, 2)])]
    fn test_unknown_templates_stay_in_the_share_expected_to_be_another_alleles(
        #[case] sample: ((i32, i32, i32), i32),
        #[case] expected: Vec<(i32, i32)>,
    ) {
        // Each class unit is one template; `u s / (m + p + s)` rounded stay: 1 * 1 / 2 keeps one.
        assert_eq!(resolve_classes(&[sample]), expected);
    }

    #[test]
    fn test_a_normal_moves_its_ambiguous_reads_once_the_tumor_keeps_none() {
        assert_eq!(
            resolve_classes(&[((1, 0, 0), 10), ((0, 0, 1), 10)]),
            vec![(1, 0), (2, 0)]
        );
        // A tumor with no reads of the row keeps none either.
        assert_eq!(
            resolve_classes(&[((0, 0, 0), 10), ((0, 0, 1), 10)]),
            vec![(0, 0), (2, 0)]
        );
    }

    #[test]
    fn test_a_separate_read_in_the_normal_counts_against_the_tumors_ambiguous_reads() {
        // Beside 20 MNV reads, 2 * 1 / 21 rounds to no template staying; beside 2, 2 * 1 / 3 to one.
        assert_eq!(
            resolve_classes(&[((0, 0, 2), 10), ((0, 1, 0), 10)]),
            vec![(4, 0), (0, 2)]
        );
        assert_eq!(
            resolve_classes(&[((0, 0, 2), 1), ((0, 1, 0), 1)]),
            vec![(2, 2), (0, 2)]
        );
    }

    #[test]
    fn test_a_proven_molecule_in_the_normal_moves_whatever_the_tumor_keeps() {
        assert_eq!(
            resolve_classes(&[((0, 1, 0), 10), ((1, 0, 0), 10)]),
            vec![(0, 2), (1, 0)]
        );
    }

    /// Surveys one fragment aligned over 1-based positions `start` to `end` of the window with this
    /// CIGAR, whose bases are the window's reference with these edits and qualities.
    fn survey_span(
        start: u64,
        end: u64,
        cigar: &str,
        edits: &[(u64, u8)],
        quality_edits: &[(u64, u8)],
    ) -> Survey {
        let (from, to) = ((start - 101) as usize, (end - 100) as usize);
        let ops: Vec<(usize, char)> = cigar
            .split_inclusive(|c: char| c.is_ascii_alphabetic())
            .map(|op| {
                (
                    op[..op.len() - 1].parse().unwrap(),
                    op.chars().last().unwrap(),
                )
            })
            .collect();
        let soft = |op: Option<&(usize, char)>| {
            op.filter(|&&(_, kind)| kind == 'S')
                .map_or(0, |&(len, _)| len)
        };
        let (leading, trailing) = (soft(ops.first()), soft(ops.last()));
        let bases = format!(
            "{}{}{}",
            "A".repeat(leading),
            &edited(edits)[from..to],
            "A".repeat(trailing)
        );
        let mut quals = vec![40u8; leading];
        quals.extend_from_slice(&self::quals(quality_edits)[from..to]);
        quals.extend(std::iter::repeat_n(40u8, trailing));
        let mut builder = SamBuilder::new().read_length(to - from);
        builder.add_frag(
            Frag::at(start as usize)
                .bases(bases)
                .quals(quals)
                .cigar(cigar),
        );
        survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap()
    }

    #[rstest]
    #[case::at_the_end(105, "5M", 1)]
    #[case::one_before_the_end(106, "6M", 2)]
    #[case::two_before_the_end(107, "7M", 3)]
    #[case::three_before_the_end(108, "8M", 0)]
    #[case::hard_clipped_end(106, "6M14H", 0)]
    #[case::soft_clipped_end(106, "6M2S", 0)]
    fn test_vardict_clips_a_mismatch_within_three_bases_of_a_match_ended_read(
        #[case] end: u64,
        #[case] cigar: &str,
        #[case] tail: usize,
    ) {
        let survey = survey_span(101, end, cigar, &[(105, b'T')], &[]);
        let read = &survey.reads[0];
        assert_eq!(read.clips(&window()), Clips { head: 0, tail });
        assert_eq!(read.aligned.clone().len(), 1);
        assert_eq!(read.clipped.clone().len(), usize::from(tail == 0));
    }

    #[rstest]
    #[case::at_the_start(105, "16M", 1)]
    #[case::one_after_the_start(104, "17M", 2)]
    #[case::two_after_the_start(103, "18M", 3)]
    #[case::three_after_the_start(102, "19M", 0)]
    #[case::hard_clipped_start(104, "3H17M", 0)]
    #[case::soft_clipped_start(104, "2S17M", 0)]
    fn test_vardict_clips_a_mismatch_within_three_bases_of_a_match_started_read(
        #[case] start: u64,
        #[case] cigar: &str,
        #[case] head: usize,
    ) {
        let survey = survey_span(start, 120, cigar, &[(105, b'T')], &[]);
        let read = &survey.reads[0];
        assert_eq!(read.clips(&window()), Clips { head, tail: 0 });
        assert_eq!(read.clipped.clone().len(), usize::from(head == 0));
    }

    #[test]
    fn test_the_clip_scan_stops_at_three_matches_in_a_row() {
        // The mismatch at 105 is met before three matches follow the one at 108: a 4-base clip, too long.
        let survey = survey_span(101, 108, "8M", &[(105, b'T'), (108, b'A')], &[]);
        assert_eq!(survey.reads[0].clips(&window()), Clips::default());
        let survey = survey_span(101, 108, "8M", &[(104, b'T'), (108, b'A')], &[]);
        assert_eq!(survey.reads[0].clips(&window()), Clips { head: 0, tail: 1 });
    }

    #[test]
    fn test_a_no_call_counts_as_a_mismatch_in_the_clip_scan() {
        let survey = survey_span(101, 107, "7M", &[(105, b'T'), (107, b'N')], &[]);
        assert_eq!(survey.reads[0].clips(&window()), Clips { head: 0, tail: 3 });
    }

    #[test]
    fn test_the_clip_scan_compares_bases_whatever_their_quality() {
        let survey = survey_span(101, 106, "6M", &[(105, b'T')], &[(104, 10), (105, 10)]);
        assert_eq!(survey.reads[0].clips(&window()), Clips { head: 0, tail: 2 });
    }

    #[test]
    fn test_a_read_with_an_indel_in_the_window_is_walked_as_aligned() {
        let reference = edited(&[(105, b'T')]);
        let mut builder = SamBuilder::new().read_length(6);
        builder.add_frag(
            Frag::at(101)
                .bases(format!("{}{}", &reference[..3], &reference[4..7]))
                .quals(vec![40u8; 6])
                .cigar("3M1D3M"),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(survey.reads[0].clips(&window()), Clips::default());
    }

    #[test]
    fn test_a_clipped_end_shortens_the_read_for_the_position_in_read() {
        // With the last 2 bases clipped from the mismatch at 119, the allele at 110 is 9 from the end.
        let survey = survey_span(101, 120, "20M", &[(110, b'T'), (119, b'A')], &[]);
        let read = &survey.reads[0];
        assert_eq!(read.aligned.clone()[0].position_in_read, 10);
        let clipped = read.clipped.clone();
        assert_eq!(clipped.len(), 1);
        assert_eq!(clipped[0].position_in_read, 9);
    }

    #[test]
    fn test_a_partial_reconciles_once_vardicts_clipping_is_replayed() {
        // One read ends at 105, which VarDict clips away; it counted the two whose mates prove the MNV.
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 2);
        let proven_read = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let proven_mate = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let mut builder = pairs_of(&[
            (
                (&proven_read.0, proven_read.1.clone()),
                (&proven_mate.0, proven_mate.1.clone()),
            ),
            (
                (&proven_read.0, proven_read.1.clone()),
                (&proven_mate.0, proven_mate.1.clone()),
            ),
        ]);
        builder.add_frag(
            Frag::at(101)
                .bases(edited(&[(105, b'T')])[..5].to_string())
                .quals(quals(&[])[..5].to_vec())
                .cigar("5M"),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        let (moved, kept, decision) = decide_survey(&survey, &mnv, &partial);
        assert_eq!((decision.unmatched, decision.unmatched_clipped), (-1, 0));
        assert!(decision.reconciled());
        assert_eq!((moved.alt_depth, kept.alt_depth), (2, 0));
    }

    #[test]
    fn test_a_short_read_clipped_at_both_ends_is_walked_without_panic() {
        let survey = survey_span(105, 107, "3M", &[(105, b'T'), (107, b'A')], &[]);
        let read = &survey.reads[0];
        assert_eq!(read.clips(&window()), Clips { head: 0, tail: 3 });
        assert_eq!(read.aligned.len(), 1);
        assert!(read.clipped.is_empty());
    }

    #[rstest]
    #[case(0, 3)]
    #[case(2, 3)]
    #[case(4, 5)]
    fn test_the_flank_covers_the_clip_scan(#[case] mnv_gap: usize, #[case] flank: usize) {
        let options = ReadOptions {
            mnv_gap,
            ..ReadOptions::default()
        };
        assert_eq!(options.flank(), flank);
    }

    #[test]
    fn test_a_separate_read_only_the_clipping_removes_keeps_the_row_whole() {
        // VarDict counted 1 read; as aligned 2 feed the row, one of them another allele's.
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 1);
        let proven_read = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let proven_mate = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let mut builder = pairs_of(&[(
            (&proven_read.0, proven_read.1.clone()),
            (&proven_mate.0, proven_mate.1.clone()),
        )]);
        builder.add_frag(
            Frag::at(101)
                .bases(edited(&[(105, b'T')])[..5].to_string())
                .quals(quals(&[])[..5].to_vec())
                .cigar("5M"),
        );
        builder.add_frag(
            Frag::at(101)
                .bases(edited(&[(105, b'T')]))
                .quals(quals(&[])),
        );
        let survey = survey(
            builder.to_source(),
            builder.header(),
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        let (moved, kept, decision) = decide_survey(&survey, &mnv, &partial);
        assert_eq!((decision.unmatched, decision.unmatched_clipped), (-2, -1));
        assert_eq!((moved.alt_depth, kept.alt_depth), (0, 1));
    }

    /// Writes one indexed BAM per read group sample name and opens them in that order.
    /// Whether a read carries the row's insertion whole, for each of these fragments, in a
    /// reference of `ACGT` repeated, whose bases at 101-120 are the window's.
    fn carries_frags(frags: Vec<Frag>, row: &TumorOnlyVariant) -> Vec<bool> {
        let (_dir, path) = crate::reference::tests::fasta(&[("chr1", &"ACGT".repeat(50))]);
        let mut reference = Reference::from_path(&path).unwrap();
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        for frag in frags {
            builder.add_frag(frag);
        }
        let options = ReadOptions::default();
        let survey = survey(builder.to_source(), builder.header(), &window(), options).unwrap();
        survey
            .reads
            .iter()
            .map(|read| read.carries(row, options, &mut reference).unwrap())
            .collect()
    }

    /// A fragment at 101 with these bases inserted after this position, at this quality, and the
    /// reference elsewhere over 101-120.
    fn inserting(anchor: usize, inserted: &str, quality: u8) -> Frag {
        let before = anchor - 100;
        let bases = format!(
            "{}{inserted}{}",
            &edited(&[])[..before],
            &edited(&[])[before..]
        );
        let mut quals = vec![40u8; bases.len()];
        quals[before..before + inserted.len()].fill(quality);
        Frag::at(101)
            .bases(bases)
            .quals(quals)
            .cigar(format!("{before}M{}I{}M", inserted.len(), 20 - before))
            .attr(Tag::EDIT_DISTANCE, inserted.len() as i32)
    }

    #[rstest]
    #[case::at_the_rows_position(110, "GTAC", 40, true)]
    #[case::slid_along_the_repeat(112, "ACGT", 40, true)]
    #[case::slid_to_a_place_that_is_not_the_same_haplotype(111, "ACGT", 40, false)]
    #[case::other_bases(110, "GTAA", 40, false)]
    #[case::shorter(110, "GTA", 40, false)]
    #[case::longer(110, "GTACG", 40, false)]
    #[case::below_q(110, "GTAC", 20, false)]
    fn test_a_read_carries_an_insertion_in_an_insertion_operation(
        #[case] anchor: usize,
        #[case] inserted: &str,
        #[case] quality: u8,
        #[case] expected: bool,
    ) {
        let frags = vec![
            inserting(anchor, inserted, quality),
            Frag::at(101).bases(edited(&[])).quals(quals(&[])),
        ];
        assert_eq!(
            carries_frags(frags, &row(110, "C", "CGTAC", 1)),
            vec![expected, false]
        );
    }

    #[rstest]
    #[case::same_place("GTAC", 110, "GTAC", true)]
    #[case::same_place_other_bases("GTAC", 110, "GTAA", false)]
    #[case::one_base_right_along_the_repeat("GTAC", 111, "TACG", true)]
    #[case::one_repeat_right("GTAC", 114, "GTAC", true)]
    #[case::one_repeat_and_two_right("GTAC", 116, "ACGT", true)]
    #[case::one_base_left("GTAC", 109, "CGTA", true)]
    #[case::two_right_but_not_rotated("GTAC", 112, "GTAC", false)]
    #[case::two_right_along_a_shorter_repeat("GTAA", 112, "AAGT", true)]
    #[case::not_a_repeat_of_the_reference("TTTT", 112, "TTTT", false)]
    #[case::past_the_contig("GTAC", 400, "GTAC", false)]
    fn test_insertions_slid_along_a_repeat_are_the_same(
        #[case] bases: &str,
        #[case] other_anchor: u64,
        #[case] other: &str,
        #[case] expected: bool,
    ) {
        let (_dir, path) = crate::reference::tests::fasta(&[("chr1", &"ACGT".repeat(50))]);
        let mut reference = Reference::from_path(&path).unwrap();
        let same = same_insertion(
            "chr1",
            (110, bases.as_bytes()),
            (other_anchor, other.as_bytes()),
            &mut reference,
        )
        .unwrap();
        assert_eq!(same, expected);
    }

    fn open_bams(dir: &std::path::Path, read_group_samples: &[&str]) -> ReadSources {
        use noodles::sam::header::record::value::map::read_group::tag::SAMPLE;
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        builder.add_frag(Frag::at(101).bases(edited(&[])).quals(quals(&[])));
        let paths: Vec<PathBuf> = read_group_samples
            .iter()
            .enumerate()
            .map(|(index, name)| {
                let mut header = builder.header().clone();
                for group in header.read_groups_mut().values_mut() {
                    group
                        .other_fields_mut()
                        .insert(SAMPLE, BString::from(*name));
                }
                let path = dir.join(format!("{index}.bam"));
                let mut writer = bam::io::Writer::new(std::fs::File::create(&path).unwrap());
                writer.write_header(&header).unwrap();
                for record in builder.to_bam_records() {
                    writer.write_record(&header, &record).unwrap();
                }
                writer.try_finish().unwrap();
                let bai = bam::fs::index(&path).unwrap();
                bam::bai::fs::write(dir.join(format!("{index}.bam.bai")), &bai).unwrap();
                path
            })
            .collect();
        ReadSources::open(&paths, ReadOptions::default()).unwrap()
    }

    #[rstest]
    #[case::own_names(&["X", "Y"], &["X", "Y"], 0)]
    #[case::other_names(&["X", "Y"], &["T", "N"], 0)]
    #[case::swapped(&["X", "Y"], &["Y", "X"], 2)]
    #[case::one_name_for_both(&["S", "S"], &["S", "N"], 0)]
    fn test_only_bams_that_look_swapped_are_warned_about(
        #[case] read_group_samples: &[&str],
        #[case] names: &[&str],
        #[case] warnings: usize,
    ) {
        let dir = tempfile::tempdir().unwrap();
        let sources = open_bams(dir.path(), read_group_samples);
        let found = sources.sample_warnings(names);
        assert_eq!(found.len(), warnings, "{found:?}");
        if warnings > 0 {
            assert!(
                found[0].contains("the sample of the BAM at position 2"),
                "{found:?}"
            );
        }
    }

    #[rstest]
    #[case("0x504", 0x504)]
    #[case("0xF00", 0xF00)]
    #[case("3840", 3840)]
    fn test_flags_parse_as_vardict_reads_them(#[case] text: &str, #[case] expected: u16) {
        assert_eq!(parse_flags(text).unwrap(), expected);
    }

    #[test]
    fn test_flags_that_are_not_a_number_are_refused() {
        assert!(parse_flags("dup").is_err());
    }

    #[test]
    fn test_bams_split_on_the_pipe() {
        assert_eq!(
            split_bams("t.bam|n.bam"),
            vec![PathBuf::from("t.bam"), PathBuf::from("n.bam")]
        );
    }

    /// A fragment aligned over 1-based `start..=end` of the window, carrying the reference with
    /// these edits, with these soft-clipped bases and qualities at its head and tail, in the order
    /// the BAM holds them.
    fn clipped(
        start: u64,
        end: u64,
        edits: &[(u64, u8)],
        head: &[(u8, u8)],
        tail: &[(u8, u8)],
    ) -> Frag {
        let (from, to) = ((start - 101) as usize, (end - 100) as usize);
        let mut bases: Vec<u8> = head.iter().map(|&(base, _)| base).collect();
        bases.extend_from_slice(&edited(edits).as_bytes()[from..to]);
        bases.extend(tail.iter().map(|&(base, _)| base));
        let mut quals: Vec<u8> = head.iter().map(|&(_, quality)| quality).collect();
        quals.extend(std::iter::repeat_n(40u8, to - from));
        quals.extend(tail.iter().map(|&(_, quality)| quality));
        let clip = |count: usize| {
            if count > 0 {
                format!("{count}S")
            } else {
                String::new()
            }
        };
        let cigar = format!("{}{}M{}", clip(head.len()), to - from, clip(tail.len()));
        Frag::at(start as usize)
            .bases(String::from_utf8(bases).unwrap())
            .quals(quals)
            .cigar(cigar)
            .attr(Tag::EDIT_DISTANCE, edits.len() as i32)
    }

    /// Bases at quality 40 as a clip, with `n` for a no-call at quality 2 and any other lowercase base
    /// at quality 12.
    fn clip(text: &str) -> Vec<(u8, u8)> {
        text.bytes()
            .map(|base| match base {
                b'n' => (b'N', 2),
                base if base.is_ascii_lowercase() => (base.to_ascii_uppercase(), 12),
                base => (base, 40),
            })
            .collect()
    }

    /// Recounts the SNV C>T at 110 from these fragments at the window, under these options, in a
    /// reference of `ACGT` repeated, whose bases at 101-120 are the window's.
    fn recount_frags(frags: Vec<Frag>, options: ReadOptions) -> Recount {
        let (_dir, path) = crate::reference::tests::fasta(&[("chr1", &"ACGT".repeat(50))]);
        let mut reference = Reference::from_path(&path).unwrap();
        let mut builder = SamBuilder::new().read_length(REFERENCE.len());
        for frag in frags {
            builder.add_frag(frag);
        }
        let survey = survey(builder.to_source(), builder.header(), &window(), options).unwrap();
        recount(
            &survey,
            &window(),
            &row(110, "C", "T", 2),
            options,
            &mut reference,
        )
        .unwrap()
    }

    fn counts(recount: Recount) -> (i32, i32, i32) {
        (recount.aligned, recount.clean_clips, recount.foreign_clips)
    }

    #[rstest]
    #[case::no_remainder("T", (0, 1, 0))]
    #[case::one_matching_base("TG", (0, 1, 0))]
    #[case::two_matching_bases("TGT", (0, 1, 0))]
    #[case::three_matching_bases("TGTA", (0, 1, 0))]
    #[case::four_matching_bases("TGTAC", (0, 1, 0))]
    #[case::a_mismatch_past_a_matching_second_base("TGA", (0, 0, 1))]
    #[case::mismatches_past_a_matching_second_base("TGACA", (0, 0, 1))]
    #[case::a_mismatch_at_the_second_base("TAT", (0, 0, 0))]
    #[case::another_base_first("AGT", (0, 0, 0))]
    #[case::six_matching_bases("TGTACG", (0, 1, 0))]
    #[case::six_bases_with_a_mismatch("TGTACT", (0, 0, 0))]
    #[case::six_bases_with_mismatches_past_a_matching_second_base("TGACAC", (0, 0, 0))]
    #[case::cut_at_a_no_call("TGnTA", (0, 1, 0))]
    #[case::cut_at_a_no_call_after_the_alt("TnACA", (0, 1, 0))]
    #[case::a_no_call_first("nTGT", (0, 0, 0))]
    #[case::a_low_quality_base_kept("TgT", (0, 1, 0))]
    #[case::cut_before_a_second_low_quality_base("TGtaAC", (0, 1, 0))]
    #[case::low_quality_bases_outnumbering_the_rest("Tgt", (0, 0, 0))]
    fn test_a_tail_clip_at_the_site_is_read_as_vardict_reads_it(
        #[case] bases: &str,
        #[case] expected: (i32, i32, i32),
    ) {
        let recount = recount_frags(
            vec![clipped(101, 109, &[], &[], &clip(bases))],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), expected);
    }

    #[rstest]
    #[case::no_remainder("T", (0, 1, 0))]
    #[case::matching_bases("GTAT", (0, 1, 0))]
    #[case::a_mismatch_past_a_matching_second_base("GAT", (0, 0, 1))]
    #[case::a_mismatch_at_the_second_base("TGT", (0, 0, 0))]
    #[case::cut_at_a_no_call("CnAT", (0, 1, 0))]
    fn test_a_head_clip_at_the_site_is_read_outward(
        #[case] bases: &str,
        #[case] expected: (i32, i32, i32),
    ) {
        let recount = recount_frags(
            vec![clipped(111, 120, &[], &clip(bases), &[])],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), expected);
    }

    #[test]
    fn test_an_aligned_alt_read_and_a_clean_clip_both_support_the_call() {
        let recount = recount_frags(
            vec![
                clipped(101, 120, &[(110, b'T')], &[], &[]),
                clipped(101, 109, &[], &[], &clip("TGT")),
            ],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), (1, 1, 0));
        assert_eq!(recount.supporting(), 2);
    }

    #[test]
    fn test_clips_piled_at_one_site_are_each_read_on_their_own() {
        let recount = recount_frags(
            vec![
                clipped(101, 109, &[], &[], &clip("TGT")),
                clipped(101, 109, &[], &[], &clip("TGA")),
                clipped(101, 109, &[], &[], &clip("TAC")),
                clipped(111, 120, &[], &clip("GAT"), &[]),
            ],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), (0, 1, 2));
    }

    #[test]
    fn test_reverse_strand_reads_are_read_the_same() {
        let recount = recount_frags(
            vec![
                clipped(101, 120, &[(110, b'T')], &[], &[]).strand(Strand::Minus),
                clipped(101, 109, &[], &[], &clip("TGA")).strand(Strand::Minus),
                clipped(111, 120, &[], &clip("TAT"), &[]).strand(Strand::Minus),
            ],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), (1, 1, 1));
    }

    #[rstest]
    #[case::at_the_end(110, &[(110, b'T')], (0, 1, 0))]
    #[case::one_matching_base_after(111, &[(110, b'T')], (0, 1, 0))]
    #[case::two_matching_bases_after(112, &[(110, b'T')], (0, 1, 0))]
    #[case::a_mismatch_at_the_end(111, &[(110, b'T'), (111, b'A')], (0, 0, 0))]
    #[case::a_mismatch_two_after(112, &[(110, b'T'), (112, b'A')], (0, 0, 1))]
    #[case::three_matching_bases_after(113, &[(110, b'T')], (1, 0, 0))]
    fn test_a_clip_vardict_makes_at_a_read_end_is_read_like_the_aligners(
        #[case] end: u64,
        #[case] edits: &[(u64, u8)],
        #[case] expected: (i32, i32, i32),
    ) {
        let recount = recount_frags(
            vec![clipped(101, end, edits, &[], &[])],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), expected);
    }

    #[rstest]
    #[case::matching_bases_before_the_alt_are_aligned(&[(b'T', 40), (b'A', 40)], (0, 1, 0))]
    #[case::not_at_a_low_quality(&[(b'T', 10), (b'A', 40)], (0, 0, 0))]
    #[case::not_through_a_mismatch(&[(b'C', 40), (b'A', 40)], (0, 0, 0))]
    fn test_vardict_extends_an_alignment_through_clipped_bases_that_match(
        #[case] before: &[(u8, u8)],
        #[case] expected: (i32, i32, i32),
    ) {
        // The clip starts at 108, where the reference reads TA before the site.
        let mut tail = before.to_vec();
        tail.extend(clip("TG"));
        let recount = recount_frags(
            vec![clipped(101, 107, &[], &[], &tail)],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), expected);
    }

    #[rstest]
    #[case::mean_at_q(&[(b'T', 20), (b'G', 40)], (0, 1, 0))]
    #[case::mean_below_q(&[(b'T', 20), (b'G', 39)], (0, 0, 0))]
    #[case::aligned_below_q(&[], (0, 0, 0))]
    fn test_a_clip_is_a_high_quality_read_by_its_mean_quality(
        #[case] tail: &[(u8, u8)],
        #[case] expected: (i32, i32, i32),
    ) {
        let options = ReadOptions {
            base_quality: 30.0,
            ..ReadOptions::default()
        };
        let frag = if tail.is_empty() {
            Frag::at(101)
                .bases(edited(&[(110, b'T')]))
                .quals(quals(&[(110, 29)]))
                .attr(Tag::EDIT_DISTANCE, 1)
        } else {
            clipped(101, 109, &[], &[], tail)
        };
        assert_eq!(counts(recount_frags(vec![frag], options)), expected);
    }

    #[test]
    fn test_a_read_feeding_an_mnv_through_the_site_does_not_support_the_snv() {
        let recount = recount_frags(
            vec![clipped(101, 120, &[(110, b'T'), (111, b'A')], &[], &[])],
            ReadOptions::default(),
        );
        assert_eq!(counts(recount), (0, 0, 0));
    }

    #[rstest]
    #[case::nothing_to_compare(b"", b"", true)]
    #[case::all_matching(b"TGTA", b"TGTA", true)]
    #[case::one_of_four(b"TGTX", b"TGTA", false)]
    #[case::one_of_seven(b"TGTACGX", b"TGTACGT", true)]
    #[case::three_of_twenty(b"XXXTTTTTTTTTTTTTTTTT", b"TTTTTTTTTTTTTTTTTTTT", false)]
    #[case::three_of_twenty_one(b"XXXTTTTTTTTTTTTTTTTTT", b"TTTTTTTTTTTTTTTTTTTTT", true)]
    #[case::four_of_forty(
        b"XXXXTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
        false
    )]
    #[case::reference_runs_out(b"TG", b"T", false)]
    fn test_matches_reference_is_vardicts_ismatchref(
        #[case] bases: &[u8],
        #[case] reference: &[u8],
        #[case] expected: bool,
    ) {
        assert_eq!(matches_reference(bases, reference), expected);
    }

    #[test]
    fn test_outward_bases_stop_at_the_contig_ends() {
        let (_dir, path) = crate::reference::tests::fasta(&[("c1", "ACGTAC")]);
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(outward(&mut reference, "c1", 3, 1, 2).unwrap(), b"GT");
        assert_eq!(outward(&mut reference, "c1", 3, -1, 2).unwrap(), b"GC");
        assert_eq!(outward(&mut reference, "c1", 5, 1, 4).unwrap(), b"AC");
        assert_eq!(outward(&mut reference, "c1", 2, -1, 4).unwrap(), b"CA");
        assert_eq!(outward(&mut reference, "c1", 7, 1, 2).unwrap(), b"");
        assert_eq!(outward(&mut reference, "c1", 0, -1, 2).unwrap(), b"");
    }
}
