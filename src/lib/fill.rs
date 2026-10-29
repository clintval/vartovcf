//! How the positions of VarDict's regions that have no reference-only record are filled with one,
//! so a gVCF of VarDict's pileup (`-p`) output has exactly one at every position VarDict piled up,
//! and the regions of a BED file VarDict wrote no row for.
//!
//! A position of a region has none when VarDict wrote only ALT rows for it, as it does at every
//! position with an allele, when the rows it wrote for it were written elsewhere or folded, when no
//! read has a base VarDict counts there inside a deletion or MNV, or when no read covers it; see
//! [`crate::site`] for what VarDict counted at each position.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use ahash::{AHashMap, AHashSet};

/// How a filled position's depth was found, the value of INFO `GVCF_FILL`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Fill {
    /// From the rows VarDict counted at the position, none a reference-only record written there.
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

/// The contig and the 1-based start and end of a VarDict region from its `chr:start-end` column.
pub fn parse_segment(segment: &str) -> Option<(&str, u64, u64)> {
    let (contig, range) = segment.rsplit_once(':')?;
    let (start, end) = range.split_once('-')?;
    Some((contig, start.parse().ok()?, end.parse().ok()?))
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
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use std::io::Write;

    use super::*;

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
