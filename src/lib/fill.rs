//! Fills the positions of VarDict's regions that have no record with reference-only records, so a
//! gVCF of VarDict's pileup (`-p`) output has a record at every position it piled up.
//!
//! VarDict counts each row's reads at one position: the first deleted base of a deletion it writes
//! with the base before, the anchor of an insertion it writes without it, and the first base of
//! every other row. A position of a region has no record when VarDict wrote rows for it that were
//! written elsewhere or folded, when no read has a base VarDict counts there inside a deletion or
//! MNV, or when no read covers it.

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};

use crate::record::{TumorOnlyVariant, variant_class};

/// How a filled position's depth was found, the value of INFO `GVCF_FILL`, ordered from the fill
/// that knows the most about the position to the one that knows the least.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Fill {
    /// From a row VarDict counted at the position, none of which was written there.
    Rows,
    /// Inside a deletion or MNV, with no row VarDict counted at the position.
    Spanned,
    /// No row VarDict counted at the position, and no deletion or MNV spans it.
    NoReads,
}

impl Fill {
    /// The value of INFO `GVCF_FILL` for this fill.
    pub fn as_str(&self) -> &'static str {
        match self {
            Fill::Rows => "ROWS",
            Fill::Spanned => "SPANNED",
            Fill::NoReads => "NO_READS",
        }
    }
}

/// A reference-only record to write at a position without one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Filled {
    /// The 1-based position.
    pub position: u64,
    /// How its depth was found.
    pub fill: Fill,
    /// The forward reads carrying the reference base.
    pub ref_forward: i32,
    /// The reverse reads carrying the reference base.
    pub ref_reverse: i32,
    /// The read depth.
    pub depth: i32,
}

/// The contig and the 1-based start and end of a VarDict region from its `chr:start-end` column.
pub fn parse_segment(segment: &str) -> Option<(&str, u64, u64)> {
    let (contig, range) = segment.rsplit_once(':')?;
    let (start, end) = range.split_once('-')?;
    Some((contig, start.parse().ok()?, end.parse().ok()?))
}

/// A row's insertion or deletion as VarDict's genotype column spells it after the slash: `-` and a
/// deletion's length, `+` and the length of an insertion's description, or `None` for any other
/// allele or a column not in VarDict's notation.
fn vardict_indel<'a>(call: &TumorOnlyVariant<'a>) -> Option<&'a str> {
    let (_, allele) = call.gt.split_once('/')?;
    (allele.starts_with('-') || allele.starts_with('+')).then_some(allele)
}

/// The position VarDict counted a row's reads at. VarDict writes a plain deletion, spelled `-` and a
/// length, from the base before its first deleted base, where it counted it, and an insertion
/// followed by matched bases, spelled `+` and its description's length, which is ALT's plus one,
/// from the base after its anchor; every other row starts where it was counted. A row whose genotype
/// column is not in VarDict's notation is read by its alleles' shape.
pub fn vardict_position(call: &TumorOnlyVariant) -> u64 {
    let digits = |length: &str| length.parse::<usize>().ok();
    match vardict_indel(call) {
        Some(deletion) if deletion.starts_with('-') && digits(&deletion[1..]).is_some() => {
            call.start + 1
        }
        Some(insertion)
            if insertion.starts_with('+')
                && !call.alt_allele.starts_with('<')
                && digits(&insertion[1..]) == Some(call.alt_allele.len() + 1) =>
        {
            call.start - 1
        }
        Some(_) => call.start,
        None => {
            let (ref_allele, alt_allele) = (call.ref_allele.as_bytes(), call.alt_allele.as_bytes());
            let anchored = ref_allele.len() > alt_allele.len()
                && alt_allele.len() == 1
                && ref_allele[0].eq_ignore_ascii_case(&alt_allele[0]);
            if anchored { call.start + 1 } else { call.start }
        }
    }
}

/// Whether VarDict took a row's reads from the REF reads of the base it counted the row at, as it
/// does for an insertion's.
fn inserts(call: &TumorOnlyVariant) -> bool {
    match vardict_indel(call) {
        Some(allele) => allele.starts_with('+') && !call.alt_allele.starts_with('<'),
        None => variant_class(call.ref_allele, call.alt_allele) == Some("Insertion"),
    }
}

/// The REF reads and depth of the row chosen at a position VarDict counted reads at.
#[derive(Clone, Copy)]
struct Counted {
    complex: bool,
    ref_forward: i32,
    ref_reverse: i32,
    depth: i32,
}

/// The fills of the positions `start..=end` of a region missing from `written`, from the region's
/// rows, made one at a time. A position with rows VarDict counted there takes its REF reads and DP
/// from the one that is not Complex, whose REF reads VarDict overcounts, with the most depth, and
/// the reads of each insertion there join its REF reads, as VarDict takes them from the anchor
/// base's; a position inside a deletion, a symbolic `<DEL>` included, or an MNV takes no REF reads
/// and the depth of the one spanning it with the most; any other position has no reads. The rows
/// and spans are held, but no position is until it is made.
pub fn region_fills<'w>(
    calls: &[&TumorOnlyVariant],
    written: &'w AHashSet<u64>,
    start: u64,
    end: u64,
) -> impl Iterator<Item = Filled> + 'w {
    let mut counted: AHashMap<u64, Counted> = AHashMap::new();
    let mut inserted: AHashMap<u64, (i32, i32)> = AHashMap::new();
    let mut spans: Vec<(u64, u64, i32)> = Vec::new();
    for &call in calls {
        let position = vardict_position(call);
        let row = Counted {
            complex: call.variant_type == "Complex",
            ref_forward: call.ref_forward,
            ref_reverse: call.ref_reverse,
            depth: call.depth,
        };
        let better = |row: &Counted| (!row.complex, row.depth);
        counted
            .entry(position)
            .and_modify(|best| {
                if better(&row) > better(best) {
                    *best = row;
                }
            })
            .or_insert(row);
        if inserts(call) {
            let reads = inserted.entry(position).or_default();
            reads.0 += call.alt_forward;
            reads.1 += call.alt_reverse;
        }
        let last = match call.alt_allele {
            "<DEL>" => Some(call.end),
            alt if !alt.starts_with('<') => Some(call.start + call.ref_allele.len() as u64 - 1),
            _ => None,
        };
        if let Some(last) = last.filter(|&last| last > position) {
            spans.push((position + 1, last, call.depth));
        }
    }
    spans.sort_unstable_by_key(|&(from, _, _)| Reverse(from));
    let mut open: BinaryHeap<(i32, u64)> = BinaryHeap::new();
    (start..=end)
        .filter(|position| !written.contains(position))
        .map(move |position| {
            while spans.last().is_some_and(|&(from, _, _)| from <= position) {
                let (_, to, depth) = spans.pop().expect("A span was just seen!");
                open.push((depth, to));
            }
            while open.peek().is_some_and(|&(_, to)| to < position) {
                open.pop();
            }
            match (counted.get(&position), open.peek()) {
                (Some(row), _) => {
                    let (forward, reverse) = inserted.get(&position).copied().unwrap_or_default();
                    Filled {
                        position,
                        fill: Fill::Rows,
                        ref_forward: row.ref_forward + forward,
                        ref_reverse: row.ref_reverse + reverse,
                        depth: row.depth,
                    }
                }
                (None, Some(&(depth, _))) => Filled {
                    position,
                    fill: Fill::Spanned,
                    ref_forward: 0,
                    ref_reverse: 0,
                    depth,
                },
                (None, None) => Filled {
                    position,
                    fill: Fill::NoReads,
                    ref_forward: 0,
                    ref_reverse: 0,
                    depth: 0,
                },
            }
        })
}

/// How VarDict turns a BED line into the region it piles up: widened by `extend` bases on each
/// side, as its `-x` does, then, unless it reads the BED's starts as 1-based, moved one base right.
/// VarDict reads them as 1-based unless it runs with `-z 1`, or reads a 4-column BED without `-c`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BedLayout {
    /// The bases VarDict's `-x` widens each side by.
    pub extend: u64,
    /// Whether VarDict reads the BED's starts as 1-based.
    pub one_based: bool,
}

impl BedLayout {
    /// The 1-based, inclusive region VarDict piles up for a BED line's start and end, as VarDict's
    /// `RegionBuilder` makes it: a 0-based start moves right only while it stays below the end. The
    /// start is 0 or below for a line near its contig's start.
    pub fn region(&self, start: u64, end: u64) -> (i64, i64) {
        let extend = self.extend as i64;
        let (start, end) = (start as i64 - extend, end as i64 + extend);
        if !self.one_based && start < end {
            (start + 1, end)
        } else {
            (start, end)
        }
    }

    /// The layout under which a BED line's start and end give a VarDict region, if one does.
    pub fn building(line: (u64, u64), region: (i64, i64)) -> Option<BedLayout> {
        let extend = u64::try_from(region.1 - line.1 as i64).ok()?;
        [false, true]
            .into_iter()
            .map(|one_based| BedLayout { extend, one_based })
            .find(|layout| layout.region(line.0, line.1) == region)
    }
}

/// One interval of a BED file, as its line gives it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BedLine {
    /// The 1-based line number.
    pub number: usize,
    /// The contig.
    pub contig: String,
    /// The start, read by a [`BedLayout`].
    pub start: u64,
    /// The end.
    pub end: u64,
}

impl std::fmt::Display for BedLayout {
    /// The options that give this layout.
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "--fill-gvcf-bed-extend {}", self.extend)?;
        if self.one_based {
            write!(f, " --fill-gvcf-bed-one-based")?;
        }
        Ok(())
    }
}

/// A BED file's lines and the regions VarDict piles up for them under a layout.
pub struct Bed {
    path: PathBuf,
    layout: BedLayout,
    lines: Vec<BedLine>,
    regions: AHashMap<String, AHashSet<(i64, i64)>>,
}

impl Bed {
    /// The lines of the BED file at `path`, with the region VarDict piles up for each under `layout`.
    pub fn new(path: &Path, layout: BedLayout, lines: Vec<BedLine>) -> Self {
        let mut regions: AHashMap<String, AHashSet<(i64, i64)>> = AHashMap::new();
        for line in &lines {
            regions
                .entry(line.contig.clone())
                .or_default()
                .insert(layout.region(line.start, line.end));
        }
        Bed {
            path: path.to_path_buf(),
            layout,
            lines,
            regions,
        }
    }

    /// Each line's contig and the region VarDict piles up for it.
    pub fn regions(&self) -> impl Iterator<Item = (&str, i64, i64)> {
        self.lines.iter().map(|line| {
            let (start, end) = self.layout.region(line.start, line.end);
            (line.contig.as_str(), start, end)
        })
    }

    /// Refuses a VarDict region that no line gives under the layout, naming the options one line
    /// would give it with, as the BED is then read otherwise than VarDict read it.
    pub fn check(&self, contig: &str, start: u64, end: u64) -> Result<(), String> {
        let region = (start as i64, end as i64);
        if self
            .regions
            .get(contig)
            .is_some_and(|regions| regions.contains(&region))
        {
            return Ok(());
        }
        let hint = self
            .lines
            .iter()
            .filter(|line| line.contig == contig)
            .find_map(|line| {
                BedLayout::building((line.start, line.end), region)
                    .map(|layout| format!("line {} gives it with {layout}", line.number))
            })
            .unwrap_or_else(|| "give it the BED file VarDict was given".to_string());
        Err(format!(
            "VarDict's region {contig}:{start}-{end} is the region of no line of the BED file {:?} with {}: {hint}!",
            self.path, self.layout
        ))
    }
}

/// Sorts 1-based, inclusive ranges and joins those that overlap or touch.
pub fn merge(mut ranges: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    ranges.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ranges.len());
    for (start, end) in ranges {
        match merged.last_mut() {
            Some(last) if start <= last.1 + 1 => last.1 = last.1.max(end),
            _ => merged.push((start, end)),
        }
    }
    merged
}

/// The parts of the 1-based, inclusive range `start..=end` outside every range of `covered`, which
/// [`merge`] has sorted and joined.
pub fn uncovered(start: u64, end: u64, covered: &[(u64, u64)]) -> Vec<(u64, u64)> {
    let mut parts = Vec::new();
    let mut next = start;
    for &(from, to) in covered {
        if to < next {
            continue;
        }
        if from > end {
            break;
        }
        if from > next {
            parts.push((next, from - 1));
        }
        next = to + 1;
    }
    if next <= end {
        parts.push((next, end));
    }
    parts
}

/// Reads a BED file's lines from the first three columns of each, split on whitespace, refusing a
/// line without a start and an end no lower than it, a contig `contig_length` does not know and an
/// end past the contig's; blank, `#`, `track` and `browser` lines are skipped.
pub fn read_bed(
    path: &Path,
    contig_length: impl Fn(&str) -> Option<u64>,
) -> Result<Vec<BedLine>, String> {
    let describe =
        |error: std::io::Error| format!("Could not read the BED file {path:?}: {error}!");
    let reader = BufReader::new(File::open(path).map_err(describe)?);
    let mut lines = Vec::new();
    for (index, line) in reader.lines().enumerate() {
        let line = line.map_err(describe)?;
        let number = index + 1;
        let mut fields = line.split_whitespace();
        let Some(contig) = fields.next() else {
            continue;
        };
        if contig.starts_with('#') || contig == "track" || contig == "browser" {
            continue;
        }
        let bounds = fields
            .next()
            .zip(fields.next())
            .and_then(|(start, end)| Some((start.parse::<u64>().ok()?, end.parse::<u64>().ok()?)));
        let Some((start, end)) = bounds.filter(|(start, end)| start <= end) else {
            return Err(format!(
                "Expected a contig, a start and an end no lower than it on line {number} of the BED file {path:?}, found '{line}'!"
            ));
        };
        let length = contig_length(contig).ok_or_else(|| {
            format!("Contig '{contig}' on line {number} of the BED file {path:?} is not in the reference!")
        })?;
        if end > length {
            return Err(format!(
                "The interval on line {number} of the BED file {path:?} ends at {end}, past the {length} bases of {contig} in the reference!"
            ));
        }
        lines.push(BedLine {
            number,
            contig: contig.to_string(),
            start,
            end,
        });
    }
    Ok(lines)
}

#[cfg(test)]
mod tests {
    use ahash::AHashSet;
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use std::io::Write;

    use super::*;
    use crate::record::TumorOnlyVariant;

    fn call(
        start: u64,
        ref_allele: &'static str,
        alt_allele: &'static str,
        depth: i32,
        ref_reads: (i32, i32),
        variant_type: &'static str,
    ) -> TumorOnlyVariant<'static> {
        TumorOnlyVariant {
            contig: "chr1",
            start,
            end: start + ref_allele.len() as u64 - 1,
            ref_allele,
            alt_allele,
            depth,
            ref_forward: ref_reads.0,
            ref_reverse: ref_reads.1,
            variant_type,
            ..TumorOnlyVariant::default()
        }
    }

    fn filled(position: u64, fill: Fill, ref_reads: (i32, i32), depth: i32) -> Filled {
        Filled {
            position,
            fill,
            ref_forward: ref_reads.0,
            ref_reverse: ref_reads.1,
            depth,
        }
    }

    #[test]
    fn test_parse_segment_splits_the_contig_from_the_range_after_the_last_colon() {
        assert_eq!(parse_segment("chr1:1857-3296"), Some(("chr1", 1857, 3296)));
        assert_eq!(
            parse_segment("HLA-A*01:01:01:01:5-9"),
            Some(("HLA-A*01:01:01:01", 5, 9))
        );
        assert_eq!(parse_segment("chr1"), None);
        assert_eq!(parse_segment("chr1:9-x"), None);
        assert_eq!(parse_segment(""), None);
    }

    #[test]
    fn test_vardict_position_is_the_first_deleted_base_of_an_anchored_deletion() {
        assert_eq!(
            vardict_position(&call(10, "CA", "C", 9, (0, 0), "Deletion")),
            11
        );
        assert_eq!(
            vardict_position(&call(10, "C", "CA", 9, (0, 0), "Insertion")),
            10
        );
        assert_eq!(vardict_position(&call(10, "C", "T", 9, (0, 0), "SNV")), 10);
        assert_eq!(
            vardict_position(&call(10, "CA", "TG", 9, (0, 0), "Complex")),
            10
        );
        assert_eq!(
            vardict_position(&call(10, "CAG", "T", 9, (0, 0), "Complex")),
            10
        );
        assert_eq!(vardict_position(&call(10, "C", "C", 9, (5, 4), "")), 10);
        assert_eq!(
            vardict_position(&call(10, "C", "<DEL>", 9, (0, 0), "DEL")),
            10
        );
    }

    #[rstest]
    #[case("CA", "C", "A/-1", 11)]
    #[case("C", "<DEL>", "C/-1000", 11)]
    #[case("TCG", "T", "T/-2T", 10)]
    #[case("TCG", "T", "T/-2mT", 10)]
    #[case("TCG", "T", "-2T/-2T", 10)]
    #[case("C", "CA", "C/+1", 10)]
    #[case("G", "AT", "C/+3", 9)]
    #[case("CC", "AATG", "C/+5", 9)]
    #[case("C", "<DUP>", "C/+160", 10)]
    fn test_vardict_position_reads_a_deletion_from_vardicts_allele_notation(
        #[case] ref_allele: &'static str,
        #[case] alt_allele: &'static str,
        #[case] gt: &'static str,
        #[case] expected: u64,
    ) {
        let row = TumorOnlyVariant {
            gt,
            ..call(10, ref_allele, alt_allele, 9, (0, 0), "Deletion")
        };
        assert_eq!(vardict_position(&row), expected);
    }

    #[test]
    fn test_region_fills_each_position_without_a_record_by_what_explains_it() {
        let rows = [
            call(2, "G", "G", 50, (20, 25), ""),
            call(3, "G", "A", 60, (30, 28), "SNV"),
            call(4, "TC", "T", 70, (33, 30), "Deletion"),
            call(4, "TCG", "T", 70, (33, 30), "Deletion"),
            call(6, "CGGA", "C", 80, (0, 0), "Deletion"),
            call(6, "CG", "AT", 90, (40, 40), "Complex"),
            call(11, "A", "A", 40, (20, 20), ""),
        ];
        let calls: Vec<&TumorOnlyVariant> = rows.iter().collect();
        let written: AHashSet<u64> = [2, 4, 6, 11].into_iter().collect();
        assert_eq!(
            region_fills(&calls, &written, 1, 12).collect::<Vec<_>>(),
            vec![
                filled(1, Fill::NoReads, (0, 0), 0),
                filled(3, Fill::Rows, (30, 28), 60),
                filled(5, Fill::Rows, (33, 30), 70),
                filled(7, Fill::Rows, (0, 0), 80),
                filled(8, Fill::Spanned, (0, 0), 80),
                filled(9, Fill::Spanned, (0, 0), 80),
                filled(10, Fill::NoReads, (0, 0), 0),
                filled(12, Fill::NoReads, (0, 0), 0),
            ]
        );
    }

    #[test]
    fn test_region_fills_prefer_a_row_that_is_not_complex_then_the_most_depth() {
        let rows = [
            call(5, "CG", "AT", 90, (70, 10), "Complex"),
            call(5, "C", "A", 60, (30, 28), "SNV"),
            call(5, "C", "T", 61, (31, 28), "SNV"),
        ];
        let calls: Vec<&TumorOnlyVariant> = rows.iter().collect();
        assert_eq!(
            region_fills(&calls, &AHashSet::new(), 5, 5).collect::<Vec<_>>(),
            vec![filled(5, Fill::Rows, (31, 28), 61)]
        );
        let complex = [call(5, "CG", "AT", 90, (70, 10), "Complex")];
        let calls: Vec<&TumorOnlyVariant> = complex.iter().collect();
        assert_eq!(
            region_fills(&calls, &AHashSet::new(), 5, 6).collect::<Vec<_>>(),
            vec![
                filled(5, Fill::Rows, (70, 10), 90),
                filled(6, Fill::Spanned, (0, 0), 90),
            ]
        );
    }

    #[test]
    fn test_region_fills_count_the_reads_of_an_insertion_as_ref_at_its_anchor() {
        let rows = [
            TumorOnlyVariant {
                alt_forward: 12,
                alt_reverse: 8,
                ..call(5, "G", "GG", 100, (40, 40), "Insertion")
            },
            TumorOnlyVariant {
                alt_forward: 3,
                alt_reverse: 2,
                ..call(5, "G", "GGG", 100, (40, 40), "Insertion")
            },
            call(5, "GA", "G", 100, (40, 40), "Deletion"),
        ];
        let calls: Vec<&TumorOnlyVariant> = rows.iter().collect();
        assert_eq!(
            region_fills(&calls, &AHashSet::new(), 5, 5).collect::<Vec<_>>(),
            vec![filled(5, Fill::Rows, (55, 50), 100)]
        );
    }

    #[test]
    fn test_region_fills_span_a_symbolic_deletion_to_its_end() {
        let rows = [TumorOnlyVariant {
            gt: "C/-5",
            end: 15,
            ..call(10, "C", "<DEL>", 90, (0, 0), "DEL")
        }];
        let calls: Vec<&TumorOnlyVariant> = rows.iter().collect();
        assert_eq!(
            region_fills(&calls, &AHashSet::new(), 11, 16).collect::<Vec<_>>(),
            vec![
                filled(11, Fill::Rows, (0, 0), 90),
                filled(12, Fill::Spanned, (0, 0), 90),
                filled(13, Fill::Spanned, (0, 0), 90),
                filled(14, Fill::Spanned, (0, 0), 90),
                filled(15, Fill::Spanned, (0, 0), 90),
                filled(16, Fill::NoReads, (0, 0), 0),
            ]
        );
    }

    #[test]
    fn test_region_fills_take_the_most_depth_of_overlapping_spans() {
        let rows = [
            TumorOnlyVariant {
                gt: "C/-4",
                ..call(10, "CGGAT", "C", 50, (0, 0), "Deletion")
            },
            TumorOnlyVariant {
                gt: "G/-1",
                ..call(11, "GG", "G", 80, (0, 0), "Deletion")
            },
            TumorOnlyVariant {
                gt: "A/-1",
                ..call(13, "AT", "A", 60, (0, 0), "Deletion")
            },
        ];
        let calls: Vec<&TumorOnlyVariant> = rows.iter().collect();
        let written: AHashSet<u64> = [11, 12, 14].into_iter().collect();
        assert_eq!(
            region_fills(&calls, &written, 11, 15).collect::<Vec<_>>(),
            vec![
                filled(13, Fill::Spanned, (0, 0), 50),
                filled(15, Fill::NoReads, (0, 0), 0),
            ]
        );
    }

    #[test]
    fn test_region_fills_count_an_insertion_before_matched_bases_at_its_anchor() {
        let rows = [TumorOnlyVariant {
            gt: "C/+3",
            alt_forward: 5,
            alt_reverse: 5,
            ..call(6, "G", "AT", 50, (20, 20), "Complex")
        }];
        let calls: Vec<&TumorOnlyVariant> = rows.iter().collect();
        assert_eq!(
            region_fills(&calls, &AHashSet::new(), 5, 7).collect::<Vec<_>>(),
            vec![
                filled(5, Fill::Rows, (25, 25), 50),
                filled(6, Fill::Spanned, (0, 0), 50),
                filled(7, Fill::NoReads, (0, 0), 0),
            ]
        );
    }

    #[rstest]
    #[case(0, false, (1094, 1112), (1095, 1112))]
    #[case(0, true, (1094, 1112), (1094, 1112))]
    #[case(10, false, (1094, 1112), (1085, 1122))]
    #[case(10, true, (5, 20), (-5, 30))]
    #[case(0, false, (7, 7), (7, 7))]
    fn test_bed_layout_builds_the_region_vardict_piles_up_for_a_line(
        #[case] extend: u64,
        #[case] one_based: bool,
        #[case] line: (u64, u64),
        #[case] expected: (i64, i64),
    ) {
        let layout = BedLayout { extend, one_based };
        assert_eq!(layout.region(line.0, line.1), expected);
    }

    #[rstest]
    #[case((1094, 1112), (1095, 1112), Some(BedLayout { extend: 0, one_based: false }))]
    #[case((1095, 1112), (1095, 1112), Some(BedLayout { extend: 0, one_based: true }))]
    #[case((1099, 1107), (1095, 1112), Some(BedLayout { extend: 5, one_based: false }))]
    #[case((1092, 1112), (1095, 1112), None)]
    #[case((1094, 1113), (1095, 1112), None)]
    fn test_bed_layout_that_builds_a_region_from_a_line(
        #[case] line: (u64, u64),
        #[case] region: (i64, i64),
        #[case] expected: Option<BedLayout>,
    ) {
        assert_eq!(BedLayout::building(line, region), expected);
    }

    fn lines(intervals: &[(&str, u64, u64)]) -> Vec<BedLine> {
        intervals
            .iter()
            .enumerate()
            .map(|(index, &(contig, start, end))| BedLine {
                number: index + 1,
                contig: contig.to_string(),
                start,
                end,
            })
            .collect()
    }

    #[test]
    fn test_bed_check_accepts_a_region_a_line_gives_under_the_layout() {
        let layout = BedLayout {
            extend: 5,
            one_based: false,
        };
        let bed = Bed::new(Path::new("x.bed"), layout, lines(&[("chr1", 1099, 1107)]));
        assert_eq!(bed.check("chr1", 1095, 1112), Ok(()));
        assert_eq!(
            bed.regions().collect::<Vec<_>>(),
            vec![("chr1", 1095, 1112)]
        );
    }

    #[rstest]
    #[case(&[("chr1", 1095, 1112)], "line 1 gives it with --fill-gvcf-bed-extend 0 --fill-gvcf-bed-one-based!")]
    #[case(&[("chr12", 5, 9), ("chr1", 1099, 1107)], "line 2 gives it with --fill-gvcf-bed-extend 5!")]
    #[case(&[("chr1", 1092, 1112)], "give it the BED file VarDict was given!")]
    #[case(&[("chr12", 1094, 1112)], "give it the BED file VarDict was given!")]
    fn test_bed_check_refuses_a_region_no_line_gives_and_hints_the_options(
        #[case] intervals: &[(&str, u64, u64)],
        #[case] hint: &str,
    ) {
        let bed = Bed::new(Path::new("x.bed"), BedLayout::default(), lines(intervals));
        assert_eq!(
            bed.check("chr1", 1095, 1112).unwrap_err(),
            format!(
                "VarDict's region chr1:1095-1112 is the region of no line of the BED file \"x.bed\" with --fill-gvcf-bed-extend 0: {hint}"
            )
        );
    }

    #[test]
    fn test_uncovered_leaves_out_the_covered_ranges() {
        let covered = [(3, 5), (8, 9), (12, 20)];
        assert_eq!(uncovered(1, 14, &covered), vec![(1, 2), (6, 7), (10, 11)]);
        assert_eq!(uncovered(4, 5, &covered), vec![]);
        assert_eq!(uncovered(21, 22, &covered), vec![(21, 22)]);
        assert_eq!(uncovered(1, 2, &[]), vec![(1, 2)]);
    }

    #[test]
    fn test_merge_sorts_and_joins_overlapping_and_adjacent_ranges() {
        assert_eq!(
            merge(vec![(8, 9), (1, 3), (2, 5), (6, 6), (11, 12)]),
            vec![(1, 6), (8, 9), (11, 12)]
        );
    }

    fn contig_length(contig: &str) -> Option<u64> {
        match contig {
            "chr1" => Some(100),
            "chr12" => Some(20),
            _ => None,
        }
    }

    fn bed(text: &str) -> Result<Vec<(String, u64, u64)>, String> {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "{text}").unwrap();
        read_bed(file.path(), contig_length).map(|lines| {
            lines
                .into_iter()
                .map(|line| (line.contig, line.start, line.end))
                .collect()
        })
    }

    #[test]
    fn test_read_bed_numbers_each_line() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "# header\nchr1\t0\t3\n\nchr12\t9\t20\n").unwrap();
        let lines = read_bed(file.path(), contig_length).unwrap();
        assert_eq!(
            lines.iter().map(|line| line.number).collect::<Vec<_>>(),
            vec![2, 4]
        );
    }

    #[test]
    fn test_read_bed_reads_each_lines_contig_start_and_end_and_skips_headers() {
        assert_eq!(
            bed("track name=x\nbrowser position chr1\n# comment\n\nchr1\t0\t3\tA\nchr12 9 20\nchr1\t5\t5\n")
                .unwrap(),
            vec![
                ("chr1".to_string(), 0, 3),
                ("chr12".to_string(), 9, 20),
                ("chr1".to_string(), 5, 5)
            ]
        );
    }

    #[rstest]
    #[case("chr1\t0\t3\nchr1\tten\t20\n", "line 2")]
    #[case("chr1\t9\t3\n", "line 1")]
    #[case("chr1\t9\n", "line 1")]
    #[case("chr1\t0\t3\nchrX\t0\t3\n", "Contig 'chrX' on line 2")]
    #[case("chr12\t10\t21\n", "ends at 21, past the 20 bases of chr12")]
    fn test_read_bed_refuses_a_line_it_cannot_fill(#[case] text: &str, #[case] expected: &str) {
        let error = bed(text).unwrap_err();
        assert!(error.contains(expected), "{error}");
    }
}
