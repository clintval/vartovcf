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

use std::collections::HashMap;
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
use noodles::sam::header::record::value::{Map, map};
use noodles::sam::{self};
use streampile::{
    AgreementStrategy, DisagreementStrategy, EntryKind, RecordSource, StreamingPileupBuilder,
};

use crate::mnv::{ReadFate, ReadFates, is_sub_event, is_substitution};
use crate::record::TumorOnlyVariant;
use crate::reference::Reference;

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
    /// the MNV has a link within this many bases of it.
    fn flank(self) -> usize {
        self.mnv_gap + 1
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
                let describe = |error| format!("Could not open the reads {path:?}: {error}!");
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
        })
    }

    /// The number of BAMs, one per sample.
    pub fn samples(&self) -> usize {
        self.samples.len()
    }

    /// Decides the fate of every partial row of a region under every MNV row it is a sub-event of,
    /// from the reads of each sample, keyed by (MNV row, partial row). Rows already written are not
    /// partials here, as they are never folded again.
    pub fn fates<'a>(
        &mut self,
        rows: &[&[TumorOnlyVariant<'a>]],
        written: &[bool],
        reference: &mut Reference,
    ) -> Result<ReadFates<'a>, String> {
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
            let mut per_partial: Vec<ReadFate<'a>> = partials
                .iter()
                .map(|_| ReadFate {
                    moved: Vec::new(),
                    kept: Vec::new(),
                })
                .collect();
            for (sample, calls) in rows[mnv].iter().enumerate() {
                let partial_calls: Vec<&TumorOnlyVariant<'a>> =
                    partials.iter().map(|&row| &rows[row][sample]).collect();
                let decided = self.decide(sample, calls, &partial_calls, reference)?;
                for (fate, (moved, kept)) in per_partial.iter_mut().zip(decided) {
                    fate.moved.push(moved);
                    fate.kept.push(kept);
                }
            }
            for (row, fate) in partials.into_iter().zip(per_partial) {
                fates.insert((mnv, row), fate);
            }
        }
        Ok(fates)
    }

    /// Decides one sample's reads of each partial of an MNV: the call of the reads that go to the
    /// MNV and the call of the reads that stay, per partial.
    fn decide<'a>(
        &mut self,
        sample: usize,
        mnv: &TumorOnlyVariant<'a>,
        partials: &[&TumorOnlyVariant<'a>],
        reference: &mut Reference,
    ) -> Result<Vec<(TumorOnlyVariant<'a>, TumorOnlyVariant<'a>)>, String> {
        let options = self.options;
        let contig = mnv.contig;
        let length = reference.length(contig)?;
        let start = (mnv.start - 1).saturating_sub(options.flank() as u64);
        let end = (mnv.end + options.flank() as u64).min(length);
        let window = Window {
            start: start as usize,
            reference: reference.bases(contig, start, end)?.to_vec(),
        };
        let source = &mut self.samples[sample];
        let region = Region::new(
            contig,
            Position::try_from(start as usize + 1).map_err(|error| error.to_string())?
                ..=Position::try_from(end as usize).map_err(|error| error.to_string())?,
        );
        let describe = |error: &dyn std::fmt::Display| {
            format!(
                "Could not read the reads of {} over {contig}:{}-{end}: {error}!",
                source.path.display(),
                start + 1
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
        let query = source
            .reader
            .query(header, &region)
            .map_err(|error| describe(&error))?;
        let survey = survey(Query(query), header, contig, &window, options)
            .map_err(|error| describe(&error))?;
        self.summary.pileups += 1;
        let decided = partials
            .iter()
            .map(|partial| {
                let (moved, kept, decision) = fate(&survey, &window, mnv, partial, options);
                if partial.alt_depth > 0 {
                    self.summary.decisions += 1;
                    self.summary.reconciled += usize::from(decision.unmatched == 0);
                    debug!(
                        "{contig}:{} {}>{} under {contig}:{} {}>{} in sample {}: VarDict counted {} reads; {} are the MNV's molecules, {} another allele's, {} prove nothing, {} unmatched; {} move and {} stay; {:.3} of the templates prove nothing there.",
                        partial.start,
                        partial.ref_allele,
                        partial.alt_allele,
                        mnv.start,
                        mnv.ref_allele,
                        mnv.alt_allele,
                        sample + 1,
                        partial.alt_depth,
                        decision.proven,
                        decision.separate,
                        decision.ambiguous,
                        decision.unmatched,
                        moved.alt_depth,
                        kept.alt_depth,
                        decision.ambiguity
                    );
                }
                (moved, kept)
            })
            .collect();
        Ok(decided)
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
    start: usize,
    reference: Vec<u8>,
}

impl Window {
    /// The column of a 1-based position, if it is in the window.
    fn column(&self, position: u64) -> Option<usize> {
        let column = position.checked_sub(self.start as u64 + 1)? as usize;
        (column < self.reference.len()).then_some(column)
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
    indel: bool,
    reverse: bool,
    mapping_quality: u8,
    mismatches: f64,
    aligned_length: usize,
    leading_clip: usize,
}

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
/// at each column, filtering reads as VarDict does.
fn survey<S: RecordSource<Record = bam::Record>>(
    source: S,
    header: &sam::Header,
    contig: &str,
    window: &Window,
    options: ReadOptions,
) -> streampile::Result<Survey> {
    let columns = window.reference.len();
    let max_mismatches = options.max_mismatches;
    let mut builder = StreamingPileupBuilder::new(source, header)?
        .min_mapping_quality(options.min_mapping_quality)
        .exclude_flags(Flags::from_bits_retain(options.exclude_flags))
        .min_base_quality(options.floor())
        .read_filter(move |record: &bam::Record| {
            substitution_mismatches(record).is_none_or(|mismatches| mismatches <= max_mismatches)
        });
    let mut survey = Survey::default();
    let mut templates: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut reads: HashMap<(usize, bool), usize> = HashMap::new();
    let mut sweep = builder.columns(contig, window.start, window.start + columns)?;
    for column in 0..columns {
        let Some(pileup) = sweep.next_pileup() else {
            break;
        };
        let pileup = pileup?;
        for entry in pileup.iter() {
            let record = entry.record();
            let name = record
                .name()
                .map_or(&b"*"[..], |name| name.as_ref())
                .to_vec();
            let next = templates.len();
            let template = *templates.entry(name).or_insert(next);
            if template == survey.templates.len() {
                survey.templates.push(TemplateView {
                    called: vec![None; columns],
                    present: vec![false; columns],
                });
            }
            survey.templates[template].present[column] = true;
            let key = (template, entry.flags().is_last_segment());
            let next = survey.reads.len();
            let read = *reads.entry(key).or_insert(next);
            if read == survey.reads.len() {
                survey.reads.push(ReadView::new(record, template, columns));
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
                EntryKind::Insertion | EntryKind::Deletion | EntryKind::Skip => view.indel = true,
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
    /// A read with no bases recorded yet, with what VarDict's per-read measures need from it.
    fn new(record: &bam::Record, template: usize, columns: usize) -> Self {
        let mut leading_clip = 0;
        let mut aligned_length = 0;
        let mut before_alignment = true;
        for op in record.cigar().iter().filter_map(|op| op.ok()) {
            match op.kind() {
                Kind::SoftClip if before_alignment => leading_clip += op.len(),
                Kind::HardClip | Kind::Pad | Kind::SoftClip => {}
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Insertion => {
                    before_alignment = false;
                    aligned_length += op.len();
                }
                Kind::Deletion | Kind::Skip => before_alignment = false,
            }
        }
        ReadView {
            template,
            bases: vec![None; columns],
            indel: false,
            reverse: record.flags().is_reverse_complemented(),
            mapping_quality: record.mapping_quality().map_or(255, u8::from),
            mismatches: substitution_mismatches(record).unwrap_or(0) as f64,
            aligned_length,
            leading_clip,
        }
    }

    /// The alleles VarDict counts this read toward over the window, by replaying its walk: an
    /// allele starts at each mismatch and grows through the following bases as an MNV while they
    /// are contiguous in the read, not no-calls, at least `-q` plus 5, and either mismatches or
    /// matches bridged within the gap to a mismatch. A read with an indel in the window feeds no
    /// substitution row, as VarDict folds a nearby mismatch into a complex allele.
    fn alleles(&self, window: &Window, options: ReadOptions) -> Vec<Found> {
        let mut found = Vec::new();
        if self.indel {
            return found;
        }
        let reference = &window.reference;
        let bases = &self.bases;
        let mut column = 0;
        while column < bases.len() {
            let Some(first) = bases[column] else {
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
                let base = bases[to]?;
                (base.query_position == bases[from]?.query_position + (to - from)).then_some(base)
            };
            if options.is_good(first.quality) {
                loop {
                    let next = last + 1;
                    let Some(base) = (next < bases.len())
                        .then(|| contiguous(last, next))
                        .flatten()
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
                        .take_while(|&to| to < bases.len())
                        .map_while(|to| contiguous(last, to).map(|base| (to, base)))
                        .find(|&(to, base)| base.base != reference[to]);
                    let Some((to, mismatch)) = bridge else {
                        break;
                    };
                    if !options.grows(mismatch.quality) {
                        break;
                    }
                    for base in bases[next..=to]
                        .iter()
                        .map(|base| base.expect("a bridged base is contiguous"))
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
            if !dropped {
                let in_read =
                    bases[last].expect("the last base is held").query_position - self.leading_clip;
                let from_end = self.aligned_length.saturating_sub(in_read);
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

    /// The row's call holding `count` of these reads, with their means, or the row's own means
    /// when none of these reads was found.
    fn call<'a>(&self, row: &TumorOnlyVariant<'a>, count: i32) -> TumorOnlyVariant<'a> {
        if self.count == 0 {
            let forward = row.alt_forward.clamp(0, count);
            return TumorOnlyVariant {
                alt_depth: count,
                alt_forward: forward,
                alt_reverse: count - forward,
                high_quality_variant_reads: if row.high_quality_variant_reads < 0 {
                    row.high_quality_variant_reads
                } else {
                    row.high_quality_variant_reads.min(count)
                },
                ..row.clone()
            };
        }
        let reads = f64::from(self.count);
        let forward = self.forward.min(count);
        let positions_vary = self.positions.windows(2).any(|pair| pair[0] != pair[1]);
        TumorOnlyVariant {
            alt_depth: count,
            alt_forward: forward,
            alt_reverse: count - forward,
            high_quality_variant_reads: self.high_quality.min(count),
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

/// How one partial's reads in one sample were decided.
#[derive(Clone, Debug, Default, PartialEq)]
struct Decision {
    /// Reads whose template is the MNV's molecule.
    proven: i32,
    /// Reads whose template is another allele's.
    separate: i32,
    /// Reads whose template proves nothing.
    ambiguous: i32,
    /// Reads VarDict counted that no read in the pileup feeds, or fed by more reads than it counted.
    unmatched: i32,
    /// The fraction of the templates at the MNV's other changed positions that prove nothing there.
    ambiguity: f64,
}

/// Decides one partial's reads in one sample: the call of the reads that go to the MNV, the call of
/// the reads that stay, and how they were decided.
///
/// A read whose template is the MNV's molecule moves whatever the MNV's count in the sample. A read
/// whose template is another allele's stays. A read whose template proves nothing follows the row's
/// other reads: it moves when no template of the row in the sample is shown to be another allele's
/// and the MNV has reads there, since a row with no read shown to be its own is leakage far more
/// often than an allele that lost the base in every read, and nothing can leak from an MNV with no
/// reads; it stays once a template of the row is shown to be another allele's, as the allele is then
/// real and a read of it that lost the other base is not told from a leak. The fraction of all
/// templates that prove nothing at those positions is reported, as it bounds how often an allele's
/// `k` reads all prove nothing: that fraction to the power `k`. A row moves at most as many reads as
/// VarDict counted on it and keeps the rest, so reads VarDict counted that no read in the pileup
/// feeds stay.
fn fate<'a>(
    survey: &Survey,
    window: &Window,
    mnv: &TumorOnlyVariant<'a>,
    partial: &TumorOnlyVariant<'a>,
    options: ReadOptions,
) -> (TumorOnlyVariant<'a>, TumorOnlyVariant<'a>, Decision) {
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
    let feeding: Vec<(&ReadView, Found, Molecule)> = survey
        .reads
        .iter()
        .flat_map(|read| {
            read.alleles(window, options)
                .into_iter()
                .filter(|allele| allele.is(partial))
                .map(move |allele| {
                    let molecule = survey.templates[read.template].molecule(others, options);
                    (read, allele, molecule)
                })
        })
        .collect();
    let mut decision = Decision::default();
    for (_, _, molecule) in &feeding {
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
    let ambiguous_leak = mnv.alt_depth > 0 && decision.separate == 0;
    let mut moved = Reads::default();
    let mut kept = Reads::default();
    for (read, allele, molecule) in &feeding {
        let goes = match molecule {
            Molecule::Mnv => true,
            Molecule::Ambiguous => ambiguous_leak,
            Molecule::Separate => false,
        };
        if goes {
            moved.add(read, allele, options);
        } else {
            kept.add(read, allele, options);
        }
    }
    decision.unmatched = partial.alt_depth - feeding.len() as i32;
    let moved_count = moved.count.min(partial.alt_depth);
    let kept_count = partial.alt_depth - moved_count;
    (
        moved.call(partial, moved_count),
        kept.call(partial, kept_count),
        decision,
    )
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
            "chr1",
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
                read.alleles(&window(), ReadOptions::default())
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
            "chr1",
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert_eq!(alleles(&survey), vec![vec![(105, "T".to_string())]]);
    }

    #[test]
    fn test_a_read_with_an_indel_in_the_window_feeds_nothing() {
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
            "chr1",
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        assert!(survey.reads[0].indel);
        assert_eq!(alleles(&survey), vec![Vec::<(u64, String)>::new()]);
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
            "chr1",
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        let found = survey.reads[0].alleles(&window(), ReadOptions::default());
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
            "chr1",
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
        let survey = survey(
            builder.to_source(),
            builder.header(),
            "chr1",
            &window(),
            options,
        )
        .unwrap();
        assert!(survey.reads.is_empty());
    }

    /// A pair of mates, each given as (bases, quals).
    type Mates<'a> = ((&'a str, Vec<u8>), (&'a str, Vec<u8>));

    /// Decides one partial under one MNV from paired reads.
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
        let survey = survey(
            builder.to_source(),
            builder.header(),
            "chr1",
            &window(),
            ReadOptions::default(),
        )
        .unwrap();
        fate(&survey, &window(), mnv, partial, ReadOptions::default())
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
    fn test_ambiguous_reads_stay_beside_a_read_shown_to_be_another_alleles() {
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
        assert_eq!((moved.alt_depth, kept.alt_depth), (0, 6));
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
    fn test_reads_vardict_counted_that_no_template_feeds_stay_with_the_rows_own_means() {
        let mnv = row(105, "AC", "TA", 10);
        let partial = row(105, "A", "T", 3);
        let read1 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[(106, 20)]));
        let read2 = (edited(&[(105, b'T'), (106, b'A')]), quals(&[]));
        let (moved, kept, reconciled) = decide_pairs(
            &[((&read1.0, read1.1), (&read2.0, read2.1))],
            &mnv,
            &partial,
        );
        assert_eq!((moved.alt_depth, kept.alt_depth, reconciled), (1, 2, false));
        assert_eq!(kept.base_quality_mean, 33.0);
        assert_eq!((kept.alt_forward, kept.alt_reverse), (2, 0));
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
}
