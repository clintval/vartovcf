//! Fills the positions of VarDict's regions that have no record with reference-only records, so a
//! gVCF of VarDict's pileup (`-p`) output has a record at every position it piled up.
//!
//! VarDict counts each row's reads at one position: the first deleted base of a deletion it writes
//! with the base before, and the first base of every other row. A position of a region has no
//! record when VarDict wrote rows for it that were written elsewhere or folded, when no read has a
//! base VarDict counts there inside a deletion or MNV, or when no read covers it.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use ahash::{AHashMap, AHashSet};

use crate::record::TumorOnlyVariant;

/// How a filled position's depth was found, the value of INFO `GVCF_FILL`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fill {
    /// From a row VarDict counted at the position, none of which was written there.
    Rows,
    /// Inside a deletion or MNV, with no read having a base VarDict counted at the position.
    Spanned,
    /// No read covers the position.
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

/// The 1-based start and end of a VarDict region from its `chr:start-end` column.
pub fn segment_bounds(segment: &str) -> Option<(u64, u64)> {
    let (_, range) = segment.rsplit_once(':')?;
    let (start, end) = range.split_once('-')?;
    Some((start.parse().ok()?, end.parse().ok()?))
}

/// The position VarDict counted a row's reads at: the first deleted base of a deletion it writes
/// with the base before as ALT, and the row's start otherwise.
pub fn vardict_position(call: &TumorOnlyVariant) -> u64 {
    let (ref_allele, alt_allele) = (call.ref_allele.as_bytes(), call.alt_allele.as_bytes());
    if ref_allele.len() > alt_allele.len()
        && alt_allele.len() == 1
        && ref_allele[0].eq_ignore_ascii_case(&alt_allele[0])
    {
        call.start + 1
    } else {
        call.start
    }
}

/// The fills of the positions `start..=end` of a region missing from `written`, from the region's
/// rows: a position with rows VarDict counted there takes its REF reads and DP from the one that
/// is not Complex, whose REF reads VarDict overcounts, with the most depth; a position inside a
/// deletion or MNV takes no REF reads and the depth of the one spanning it with the most; any other
/// position has no reads.
pub fn region_fills(
    calls: &[&TumorOnlyVariant],
    written: &AHashSet<u64>,
    start: u64,
    end: u64,
) -> Vec<Filled> {
    let mut counted: AHashMap<u64, &TumorOnlyVariant> = AHashMap::new();
    let mut spanned: AHashMap<u64, i32> = AHashMap::new();
    let better = |call: &TumorOnlyVariant| (call.variant_type != "Complex", call.depth);
    for &call in calls {
        let position = vardict_position(call);
        counted
            .entry(position)
            .and_modify(|best| {
                if better(call) > better(best) {
                    *best = call;
                }
            })
            .or_insert(call);
        if call.ref_allele.len() > 1 && !call.alt_allele.starts_with('<') {
            let last = call.start + call.ref_allele.len() as u64 - 1;
            for inside in position + 1..=last {
                let depth = spanned.entry(inside).or_insert(call.depth);
                *depth = (*depth).max(call.depth);
            }
        }
    }
    (start..=end)
        .filter(|position| !written.contains(position))
        .map(
            |position| match (counted.get(&position), spanned.get(&position)) {
                (Some(row), _) => Filled {
                    position,
                    fill: Fill::Rows,
                    ref_forward: row.ref_forward,
                    ref_reverse: row.ref_reverse,
                    depth: row.depth,
                },
                (None, Some(&depth)) => Filled {
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
            },
        )
        .collect()
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

/// Reads a BED file's intervals as contig, 1-based start and inclusive end, from the first three
/// columns of each line split on whitespace; blank, `#`, `track` and `browser` lines are skipped,
/// and so are empty intervals.
pub fn read_bed(path: &Path) -> Result<Vec<(String, u64, u64)>, String> {
    let describe =
        |error: std::io::Error| format!("Could not read the BED file {path:?}: {error}!");
    let reader = BufReader::new(File::open(path).map_err(describe)?);
    let mut intervals = Vec::new();
    for (index, line) in reader.lines().enumerate() {
        let line = line.map_err(describe)?;
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
        match bounds {
            Some((start, end)) if start < end => {
                intervals.push((contig.to_string(), start + 1, end))
            }
            Some((start, end)) if start == end => {}
            _ => {
                return Err(format!(
                    "Expected a contig, a 0-based start and an end no lower than it on line {} of the BED file {path:?}, found '{line}'!",
                    index + 1
                ));
            }
        }
    }
    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use ahash::AHashSet;
    use pretty_assertions::assert_eq;
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
    fn test_segment_bounds_reads_the_start_and_end_after_the_last_colon() {
        assert_eq!(segment_bounds("chr1:1857-3296"), Some((1857, 3296)));
        assert_eq!(segment_bounds("HLA-A*01:01:01:01:5-9"), Some((5, 9)));
        assert_eq!(segment_bounds("chr1"), None);
        assert_eq!(segment_bounds("chr1:9-x"), None);
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
            region_fills(&calls, &written, 1, 12),
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
            region_fills(&calls, &AHashSet::new(), 5, 5),
            vec![filled(5, Fill::Rows, (31, 28), 61)]
        );
        let complex = [call(5, "CG", "AT", 90, (70, 10), "Complex")];
        let calls: Vec<&TumorOnlyVariant> = complex.iter().collect();
        assert_eq!(
            region_fills(&calls, &AHashSet::new(), 5, 6),
            vec![
                filled(5, Fill::Rows, (70, 10), 90),
                filled(6, Fill::Spanned, (0, 0), 90),
            ]
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

    #[test]
    fn test_read_bed_reads_1_based_inclusive_intervals_and_skips_headers() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(
            file,
            "track name=x\nbrowser position chr1\n# comment\n\nchr1\t0\t3\tA\nchr12 9 12\nchr1\t5\t5\n"
        )
        .unwrap();
        assert_eq!(
            read_bed(file.path()).unwrap(),
            vec![("chr1".to_string(), 1, 3), ("chr12".to_string(), 10, 12)]
        );
    }

    #[test]
    fn test_read_bed_refuses_a_line_without_a_start_and_end() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        write!(file, "chr1\t0\t3\nchr1\tten\t20\n").unwrap();
        let error = read_bed(file.path()).unwrap_err();
        assert!(error.contains("line 2"), "{error}");
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "chr1\t9\t3").unwrap();
        assert!(read_bed(file.path()).is_err());
    }
}
