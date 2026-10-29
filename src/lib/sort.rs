//! Sorts output records by coordinate in bounded memory with `spillsort`, spilling to disk, with
//! the VarDict regions and what VarDict counted in each row, so the writer reads them all in one
//! pass in position order.
use std::path::{Path, PathBuf};

use csv::{Position, StringRecord};
use spillsort::{Config, Merged, Sorter, key};

use crate::record::TumorOnlyVariant;

/// The most memory the records being sorted are held in before they spill to the temporary
/// directory; `spillsort` peaks at about a quarter more.
pub const MAX_BYTES: usize = 64 << 20;

/// The most memory the rows sorted to fold are held in before they spill, a quarter of the records'
/// as their last rows are still held while the records fill.
pub const MAX_ROW_BYTES: usize = MAX_BYTES / 4;

/// The kinds of entries, in the order they come at one position.
const REGION: u8 = 0;
const COUNTS: u8 = 1;
const RECORD: u8 = 2;

/// Records sorted by contig, in the header's order, then POS, REF, ALT and END, the call they make,
/// among the VarDict regions starting at each position and what VarDict counted there. Of the
/// records of one call, only the first by [`Rank`] is read back.
pub struct RecordSorter {
    sorter: Sorter,
    key: Vec<u8>,
    value: Vec<u8>,
    counted: u64,
}

/// The sort key of a record: its contig's index in the header, its 0-based POS, REF, ALT, empty
/// for a reference-only record, and END.
#[derive(Clone, Copy, Debug)]
pub struct SortKey<'a> {
    /// The index of the record's contig in the header.
    pub rid: u32,
    /// The record's 0-based position.
    pub pos: u64,
    /// The REF allele.
    pub ref_allele: &'a [u8],
    /// The ALT allele, empty for a reference-only record.
    pub alt_allele: &'a [u8],
    /// The 1-based position of the record's last REF base, or its END.
    pub end: u64,
}

/// Which record of a call is written, ordered best first: the first VarDict region's, then the most
/// ALT reads in the tumor (or only) sample, then in all samples, then the most depth in the tumor,
/// then VarDict's lowest position, REF and ALT for the row, so the choice never depends on the order
/// of the rows within a region.
#[derive(Clone, Copy, Debug)]
pub struct Rank<'a> {
    region: u64,
    tumor_reads: u64,
    reads: u64,
    tumor_depth: u64,
    start: u64,
    ref_allele: &'a [u8],
    alt_allele: &'a [u8],
}

impl<'a> Rank<'a> {
    /// The rank of a record with these calls from the region with this ordinal, written for a row
    /// with this tumor call.
    pub fn new(region: u64, calls: &[TumorOnlyVariant], row: &TumorOnlyVariant<'a>) -> Self {
        let reads = |call: &TumorOnlyVariant| call.alt_depth.max(0) as u64;
        Rank {
            region,
            tumor_reads: reads(&calls[0]),
            reads: calls.iter().map(reads).sum(),
            tumor_depth: calls[0].depth.max(0) as u64,
            start: row.start,
            ref_allele: row.ref_allele.as_bytes(),
            alt_allele: row.alt_allele.as_bytes(),
        }
    }
}

/// One entry read back from the sort, on the contig with index `rid` at a 1-based position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Entry<'a> {
    /// A region starting at the position and ending at `end`, inclusive: a VarDict region, or one
    /// only a BED file names.
    Region {
        /// The contig's index in the header.
        rid: u32,
        /// The region's first position.
        start: u64,
        /// The region's last position.
        end: u64,
    },
    /// What VarDict counted at the position in one row of the region with this ordinal.
    Counts {
        /// The contig's index in the header.
        rid: u32,
        /// The position VarDict counted the row at.
        position: u64,
        /// The ordinal of the row's VarDict region.
        region: u64,
        /// The counts, as the row's writer encoded them.
        counts: &'a [u8],
    },
    /// The VCF line of the first record of a call at the position, its POS.
    Record {
        /// The contig's index in the header.
        rid: u32,
        /// The record's POS.
        position: u64,
        /// Whether the record has no ALT allele.
        reference_only: bool,
        /// The VCF line.
        line: &'a [u8],
    },
}

impl Entry<'_> {
    /// The contig's index in the header and the position of the entry.
    pub fn locus(&self) -> (u32, u64) {
        match *self {
            Entry::Region { rid, start, .. } => (rid, start),
            Entry::Counts { rid, position, .. } | Entry::Record { rid, position, .. } => {
                (rid, position)
            }
        }
    }
}

impl RecordSorter {
    /// Makes a sorter that spills to a private directory inside `tmp_dir`.
    pub fn new(tmp_dir: &Path) -> Result<Self, String> {
        let config = Config::default().max_bytes(MAX_BYTES).tmp_dir(tmp_dir);
        let sorter = Sorter::new(config).map_err(describe)?;
        Ok(RecordSorter {
            sorter,
            key: Vec::new(),
            value: Vec::new(),
            counted: 0,
        })
    }

    /// Starts an entry's key and value at a position of the contig with index `rid`.
    fn begin(&mut self, rid: u32, position: u64, kind: u8) {
        self.key.clear();
        key::push_uint(&mut self.key, rid);
        key::push_uint(&mut self.key, position);
        key::push_uint(&mut self.key, kind);
        self.value.clear();
        self.value.push(kind);
        self.value.extend_from_slice(&rid.to_le_bytes());
        self.value.extend_from_slice(&position.to_le_bytes());
    }

    /// Adds a record's VCF line under the call it makes, ranked against the call's other records.
    pub fn push(&mut self, call: SortKey, rank: Rank, line: &[u8]) -> Result<(), String> {
        self.begin(call.rid, call.pos + 1, RECORD);
        key::push_bytes(&mut self.key, call.ref_allele);
        key::push_bytes(&mut self.key, call.alt_allele);
        key::push_uint(&mut self.key, call.end);
        let call_length = u32::try_from(self.key.len()).map_err(|_| describe_length())?;
        key::push_uint(&mut self.key, rank.region);
        key::push_uint(&mut self.key, u64::MAX - rank.tumor_reads);
        key::push_uint(&mut self.key, u64::MAX - rank.reads);
        key::push_uint(&mut self.key, u64::MAX - rank.tumor_depth);
        key::push_uint(&mut self.key, rank.start);
        key::push_bytes(&mut self.key, rank.ref_allele);
        key::push_bytes(&mut self.key, rank.alt_allele);
        self.value.push(u8::from(call.alt_allele.is_empty()));
        self.value.extend_from_slice(&call_length.to_le_bytes());
        self.value.extend_from_slice(line);
        self.sorter.push(&self.key, &self.value).map_err(describe)
    }

    /// Adds a region of the contig with index `rid` from `start` to `end`, 1-based and inclusive.
    pub fn push_region(&mut self, rid: u32, start: u64, end: u64) -> Result<(), String> {
        self.begin(rid, start, REGION);
        key::push_uint(&mut self.key, end);
        self.value.extend_from_slice(&end.to_le_bytes());
        self.sorter.push(&self.key, &self.value).map_err(describe)
    }

    /// Adds what VarDict counted in one row at a 1-based position of the contig with index `rid`,
    /// in the region with this ordinal; a position's counts come back by region, then as added.
    pub fn push_counts(
        &mut self,
        rid: u32,
        position: u64,
        region: u64,
        counts: &[u8],
    ) -> Result<(), String> {
        self.begin(rid, position, COUNTS);
        key::push_uint(&mut self.key, region);
        key::push_uint(&mut self.key, self.counted);
        self.counted += 1;
        self.value.extend_from_slice(&region.to_le_bytes());
        self.value.extend_from_slice(counts);
        self.sorter.push(&self.key, &self.value).map_err(describe)
    }

    /// Finishes the sort, returning its entries in order.
    pub fn entries(self) -> Result<Entries, String> {
        Ok(Entries {
            merged: self.sorter.finish().map_err(describe)?,
            previous: Vec::new(),
            value: Vec::new(),
        })
    }
}

/// The sorted entries, read back one at a time, with only the first record of each call.
pub struct Entries {
    merged: Merged,
    previous: Vec<u8>,
    value: Vec<u8>,
}

impl Entries {
    /// The next entry, or `None` after the last.
    pub fn next_entry(&mut self) -> Result<Option<Entry<'_>>, String> {
        loop {
            let Some((key, value)) = self.merged.next_record().map_err(describe)? else {
                return Ok(None);
            };
            if value[0] == RECORD {
                let length = u32::from_le_bytes(value[14..18].try_into().unwrap()) as usize;
                let call = &key[..length];
                if self.previous == call {
                    continue;
                }
                self.previous.clear();
                self.previous.extend_from_slice(call);
            }
            self.value.clear();
            self.value.extend_from_slice(value);
            break;
        }
        let value = &self.value;
        let rid = u32::from_le_bytes(value[1..5].try_into().unwrap());
        let position = u64::from_le_bytes(value[5..13].try_into().unwrap());
        let entry = match value[0] {
            REGION => Entry::Region {
                rid,
                start: position,
                end: u64::from_le_bytes(value[13..21].try_into().unwrap()),
            },
            COUNTS => Entry::Counts {
                rid,
                position,
                region: u64::from_le_bytes(value[13..21].try_into().unwrap()),
                counts: &value[21..],
            },
            _ => Entry::Record {
                rid,
                position,
                reference_only: value[13] == 1,
                line: &value[18..],
            },
        };
        Ok(Some(entry))
    }
}

/// VarDict's rows sorted region by region, in the order the regions were read, then by the start
/// VarDict wrote each row at, then in the order they were read, spilling to disk, so a region's rows
/// can be folded a cluster of overlapping rows at a time rather than held whole.
pub struct RowSorter {
    sorter: Sorter,
    key: Vec<u8>,
    value: Vec<u8>,
    rows: u64,
}

impl RowSorter {
    /// Makes a sorter that spills to a private directory inside `tmp_dir`.
    pub fn new(tmp_dir: &Path) -> Result<Self, String> {
        let config = Config::default().max_bytes(MAX_ROW_BYTES).tmp_dir(tmp_dir);
        let sorter = Sorter::new(config).map_err(describe)?;
        Ok(RowSorter {
            sorter,
            key: Vec::new(),
            value: Vec::new(),
            rows: 0,
        })
    }

    /// Adds a row of the region with this ordinal that VarDict wrote at `start`, with its fields
    /// and its place in the input.
    pub fn push(&mut self, region: u64, start: u64, row: &StringRecord) -> Result<(), String> {
        self.key.clear();
        key::push_uint(&mut self.key, region);
        key::push_uint(&mut self.key, start);
        key::push_uint(&mut self.key, self.rows);
        self.rows += 1;
        let position = row.position().cloned().unwrap_or_else(Position::new);
        self.value.clear();
        self.value.extend_from_slice(&region.to_le_bytes());
        for place in [position.byte(), position.line(), position.record()] {
            self.value.extend_from_slice(&place.to_le_bytes());
        }
        for (index, field) in row.iter().enumerate() {
            if index > 0 {
                self.value.push(b'\t');
            }
            self.value.extend_from_slice(field.as_bytes());
        }
        self.sorter.push(&self.key, &self.value).map_err(describe)
    }

    /// Finishes the sort, returning its rows in order.
    pub fn rows(self) -> Result<SortedRows, String> {
        Ok(SortedRows {
            merged: self.sorter.finish().map_err(describe)?,
            row: StringRecord::new(),
        })
    }
}

/// The sorted rows, read back one at a time.
pub struct SortedRows {
    merged: Merged,
    row: StringRecord,
}

impl SortedRows {
    /// The next row and the ordinal of its region, or `None` after the last.
    pub fn next_row(&mut self) -> Result<Option<(u64, &StringRecord)>, String> {
        let Some((_, value)) = self.merged.next_record().map_err(describe)? else {
            return Ok(None);
        };
        let word = |at: usize| u64::from_le_bytes(value[at..at + 8].try_into().unwrap());
        let region = word(0);
        let mut position = Position::new();
        position
            .set_byte(word(8))
            .set_line(word(16))
            .set_record(word(24));
        let fields = std::str::from_utf8(&value[32..])
            .map_err(|error| format!("Could not read back a sorted row: {error}!"))?;
        self.row.clear();
        for field in fields.split('\t') {
            self.row.push_field(field);
        }
        self.row.set_position(Some(position));
        Ok(Some((region, &self.row)))
    }
}

/// Describes an error from sorting the output records.
fn describe(error: spillsort::Error) -> String {
    format!("Could not sort the output records: {error}!")
}

/// Describes a record whose call is too long to sort.
fn describe_length() -> String {
    "Could not sort the output records: a call's alleles are 4 GiB or longer!".to_string()
}

/// Resolves the directory to spill sorted records to: `tmp_dir` when given, else a non-empty
/// `TMPDIR`, else the system's temporary directory, made absolute and checked to be a directory a
/// file can be made in.
pub fn resolve_tmp_dir(tmp_dir: Option<&Path>) -> Result<PathBuf, String> {
    let dir = match tmp_dir {
        Some(dir) => dir.to_path_buf(),
        None => match std::env::var_os("TMPDIR").filter(|dir| !dir.is_empty()) {
            Some(dir) => PathBuf::from(dir),
            None => Some(std::env::temp_dir())
                .filter(|dir| !dir.as_os_str().is_empty())
                .unwrap_or_else(|| PathBuf::from("/tmp")),
        },
    };
    let describe = |error: &dyn std::fmt::Display| {
        format!("Could not use the temporary directory {dir:?}: {error}!")
    };
    let resolved = dir.canonicalize().map_err(|error| describe(&error))?;
    if !resolved.is_dir() {
        return Err(describe(&"it is not a directory"));
    }
    tempfile::NamedTempFile::new_in(&resolved).map_err(|error| describe(&error))?;
    Ok(resolved)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use tempfile::TempDir;

    use super::*;

    fn key<'a>(rid: u32, pos: u64, ref_allele: &'a str, alt_allele: &'a str) -> SortKey<'a> {
        SortKey {
            rid,
            pos,
            ref_allele: ref_allele.as_bytes(),
            alt_allele: alt_allele.as_bytes(),
            end: pos + ref_allele.len() as u64,
        }
    }

    fn rank(region: u64, tumor_reads: u64, tumor_depth: u64, start: u64) -> Rank<'static> {
        Rank {
            region,
            tumor_reads,
            reads: tumor_reads,
            tumor_depth,
            start,
            ref_allele: b"A",
            alt_allele: b"C",
        }
    }

    fn sorted(records: &[(SortKey, &str)]) -> String {
        ranked(
            &records
                .iter()
                .map(|(key, line)| (*key, rank(0, 0, 0, 0), *line))
                .collect::<Vec<_>>(),
        )
    }

    /// The lines of the records read back from a sorter, in order.
    fn lines(sorter: RecordSorter) -> String {
        let mut entries = sorter.entries().unwrap();
        let mut out = Vec::new();
        while let Some(entry) = entries.next_entry().unwrap() {
            if let Entry::Record { line, .. } = entry {
                out.extend_from_slice(line);
            }
        }
        String::from_utf8(out).unwrap()
    }

    fn ranked(records: &[(SortKey, Rank, &str)]) -> String {
        let dir = TempDir::new().unwrap();
        let mut sorter = RecordSorter::new(dir.path()).unwrap();
        for (key, rank, line) in records {
            sorter.push(*key, *rank, line.as_bytes()).unwrap();
        }
        lines(sorter)
    }

    #[test]
    fn test_sorts_by_contig_then_position_then_alleles() {
        let lines = sorted(&[
            (key(1, 5, "A", "C"), "c1:5\n"),
            (key(0, 900, "A", "C"), "c0:900\n"),
            (key(0, 100, "AC", "A"), "c0:100:AC>A\n"),
            (key(0, 100, "A", "C"), "c0:100:A>C\n"),
            (key(0, 100, "A", ""), "c0:100:A\n"),
            (key(0, 100, "A", "<DEL>"), "c0:100:A><DEL>\n"),
            (key(0, 99, "G", "T"), "c0:99\n"),
        ]);
        assert_eq!(
            lines,
            "c0:99\nc0:100:A\nc0:100:A><DEL>\nc0:100:A>C\nc0:100:AC>A\nc0:900\nc1:5\n"
        );
    }

    #[test]
    fn test_writes_only_the_first_record_of_a_key() {
        let lines = sorted(&[
            (key(0, 100, "A", "C"), "first\n"),
            (key(0, 50, "A", "C"), "other\n"),
            (key(0, 100, "A", "C"), "second\n"),
        ]);
        assert_eq!(lines, "other\nfirst\n");
    }

    #[test]
    fn test_writes_the_first_regions_best_supported_record_of_a_call() {
        let call = key(0, 100, "A", "C");
        let lines = ranked(&[
            (call, rank(2, 90, 100, 1), "later region\n"),
            (call, rank(1, 5, 100, 1), "fewer tumor reads\n"),
            (call, rank(1, 50, 100, 3), "less depth\n"),
            (call, rank(1, 50, 400, 7), "later position\n"),
            (call, rank(1, 50, 400, 3), "best\n"),
        ]);
        assert_eq!(lines, "best\n");
    }

    #[test]
    fn test_reads_regions_then_counts_then_records_at_each_position() {
        let dir = TempDir::new().unwrap();
        let mut sorter = RecordSorter::new(dir.path()).unwrap();
        sorter
            .push(key(0, 99, "A", "C"), rank(0, 0, 0, 0), b"record\n")
            .unwrap();
        sorter
            .push(key(0, 99, "A", ""), rank(0, 0, 0, 0), b"reference\n")
            .unwrap();
        sorter.push_counts(0, 100, 2, b"second region").unwrap();
        sorter.push_counts(0, 100, 1, b"first region").unwrap();
        sorter
            .push_counts(0, 100, 1, b"first region again")
            .unwrap();
        sorter.push_region(0, 100, 120).unwrap();
        sorter.push_region(0, 101, 103).unwrap();
        sorter.push_counts(0, 99, 1, b"before").unwrap();
        let mut entries = sorter.entries().unwrap();
        let mut read = Vec::new();
        while let Some(entry) = entries.next_entry().unwrap() {
            read.push(format!("{entry:?}"));
        }
        let counts = |position: u64, region: u64, counts: &[u8]| {
            format!(
                "{:?}",
                Entry::Counts {
                    rid: 0,
                    position,
                    region,
                    counts
                }
            )
        };
        let record = |reference_only: bool, line: &[u8]| {
            format!(
                "{:?}",
                Entry::Record {
                    rid: 0,
                    position: 100,
                    reference_only,
                    line
                }
            )
        };
        let region = |start: u64, end: u64| format!("{:?}", Entry::Region { rid: 0, start, end });
        assert_eq!(
            read,
            vec![
                counts(99, 1, b"before"),
                region(100, 120),
                counts(100, 1, b"first region"),
                counts(100, 1, b"first region again"),
                counts(100, 2, b"second region"),
                record(true, b"reference\n"),
                record(false, b"record\n"),
                region(101, 103),
            ]
        );
    }

    #[test]
    fn test_reads_rows_back_by_region_then_start_then_as_added() {
        let dir = TempDir::new().unwrap();
        let mut sorter = RowSorter::new(dir.path()).unwrap();
        let row = |fields: &[&str], line: u64| {
            let mut row = csv::StringRecord::from(fields.to_vec());
            let mut position = csv::Position::new();
            position
                .set_line(line)
                .set_byte(line * 10)
                .set_record(line - 1);
            row.set_position(Some(position));
            row
        };
        sorter.push(2, 5, &row(&["region 2"], 1)).unwrap();
        sorter.push(1, 9, &row(&["region 1", "at 9"], 2)).unwrap();
        sorter
            .push(1, 3, &row(&["region 1", "at 3", ""], 3))
            .unwrap();
        sorter
            .push(1, 9, &row(&["region 1", "at 9 again"], 4))
            .unwrap();
        let mut rows = sorter.rows().unwrap();
        let mut read = Vec::new();
        while let Some((region, row)) = rows.next_row().unwrap() {
            let position = row.position().unwrap();
            read.push((
                region,
                row.iter().collect::<Vec<_>>().join("|"),
                (position.line(), position.byte(), position.record()),
            ));
        }
        assert_eq!(
            read,
            vec![
                (1, "region 1|at 3|".to_string(), (3, 30, 2)),
                (1, "region 1|at 9".to_string(), (2, 20, 1)),
                (1, "region 1|at 9 again".to_string(), (4, 40, 3)),
                (2, "region 2".to_string(), (1, 10, 0)),
            ]
        );
    }

    #[test]
    fn test_spills_rows_and_reads_them_back_in_order_with_their_places() {
        let dir = TempDir::new().unwrap();
        let mut sorter = RowSorter::new(dir.path()).unwrap();
        let field = "x".repeat(1 << 16);
        for line in 1..=512u64 {
            let mut row = csv::StringRecord::from(vec![line.to_string(), field.clone()]);
            let mut position = csv::Position::new();
            position
                .set_line(line)
                .set_byte(line * 7)
                .set_record(line - 1);
            row.set_position(Some(position));
            sorter.push(line % 2, 1000 - line, &row).unwrap();
        }
        assert!(sorter.sorter.runs() > 0);
        let mut rows = sorter.rows().unwrap();
        let mut read = Vec::new();
        while let Some((region, row)) = rows.next_row().unwrap() {
            let line: u64 = row[0].parse().unwrap();
            let position = row.position().unwrap();
            assert_eq!(&row[1], field);
            assert_eq!(
                (position.line(), position.byte(), position.record()),
                (line, line * 7, line - 1)
            );
            read.push((region, line));
        }
        let mut expected: Vec<(u64, u64)> = (1..=512).map(|line| (line % 2, line)).collect();
        expected.sort_by_key(|&(region, line)| (region, 1000 - line));
        assert_eq!(read, expected);
        drop(rows);
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn test_symbolic_alleles_with_different_ends_are_both_written() {
        let short = SortKey {
            end: 200,
            ..key(0, 100, "A", "<DEL>")
        };
        let long = SortKey {
            end: 300,
            ..key(0, 100, "A", "<DEL>")
        };
        let lines = sorted(&[(long, "long\n"), (short, "short\n"), (long, "again\n")]);
        assert_eq!(lines, "short\nlong\n");
    }

    #[test]
    fn test_spills_and_merges_more_records_than_fit_in_memory() {
        let dir = TempDir::new().unwrap();
        let mut sorter = RecordSorter::new(dir.path()).unwrap();
        let line = vec![b'x'; 1 << 16];
        for pos in (0..2048u64).rev() {
            sorter
                .push(key(0, pos, "A", "C"), rank(0, 0, 0, 0), &line)
                .unwrap();
        }
        assert!(sorter.sorter.runs() > 0);
        assert_eq!(lines(sorter).len(), 2048 * line.len());
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn test_resolves_a_given_directory_to_an_absolute_path() {
        let dir = TempDir::new().unwrap();
        let resolved = resolve_tmp_dir(Some(dir.path())).unwrap();
        assert!(resolved.is_absolute());
        assert_eq!(resolved, dir.path().canonicalize().unwrap());
    }

    #[cfg(unix)]
    #[test]
    fn test_refuses_a_directory_it_cannot_write_to() {
        use std::os::unix::fs::PermissionsExt;
        let dir = TempDir::new().unwrap();
        let locked = dir.path().join("locked");
        std::fs::create_dir(&locked).unwrap();
        std::fs::set_permissions(&locked, std::fs::Permissions::from_mode(0o500)).unwrap();
        let result = resolve_tmp_dir(Some(&locked));
        std::fs::set_permissions(&locked, std::fs::Permissions::from_mode(0o700)).unwrap();
        assert!(
            result
                .unwrap_err()
                .starts_with(&format!("Could not use the temporary directory {locked:?}"))
        );
    }

    #[test]
    fn test_refuses_a_missing_directory_or_a_file() {
        let dir = TempDir::new().unwrap();
        let missing = dir.path().join("missing");
        assert!(
            resolve_tmp_dir(Some(&missing))
                .unwrap_err()
                .starts_with("Could not use the temporary directory")
        );
        let file = dir.path().join("file");
        std::fs::write(&file, "").unwrap();
        assert!(
            resolve_tmp_dir(Some(&file))
                .unwrap_err()
                .ends_with("it is not a directory!")
        );
    }
}
