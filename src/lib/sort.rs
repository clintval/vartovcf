//! Sorts output records by coordinate in bounded memory with `spillsort`, spilling to disk.
use std::io::Write;
use std::path::{Path, PathBuf};

use spillsort::{Config, Merged, Sorter, key};

use crate::record::TumorOnlyVariant;

/// The most memory the records being sorted are held in before they spill to the temporary
/// directory; `spillsort` peaks at about a quarter more.
pub const MAX_BYTES: usize = 64 << 20;

/// Records sorted by contig, in the header's order, then POS, REF, ALT and END, the call they make.
/// Of the records of one call, only the first by [`Rank`] is written.
pub struct RecordSorter {
    sorter: Sorter,
    key: Vec<u8>,
    value: Vec<u8>,
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

/// Which record of a call is written, ordered best first: one VarDict wrote before one filled in
/// for a position without a record, then the first VarDict region's, then the most ALT reads in the
/// tumor (or only) sample, then in all samples, then the most depth in the tumor, then VarDict's
/// lowest position, REF and ALT for the row, so the choice never depends on the order of the rows
/// within a region.
#[derive(Clone, Copy, Debug)]
pub struct Rank<'a> {
    filled: bool,
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
            filled: false,
            region,
            tumor_reads: reads(&calls[0]),
            reads: calls.iter().map(reads).sum(),
            tumor_depth: calls[0].depth.max(0) as u64,
            start: row.start,
            ref_allele: row.ref_allele.as_bytes(),
            alt_allele: row.alt_allele.as_bytes(),
        }
    }

    /// This rank for a record filled in for a position without one.
    pub fn filled(self) -> Self {
        Rank {
            filled: true,
            ..self
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
        })
    }

    /// Adds a record's VCF line under the call it makes, ranked against the call's other records.
    pub fn push(&mut self, call: SortKey, rank: Rank, line: &[u8]) -> Result<(), String> {
        self.key.clear();
        key::push_uint(&mut self.key, call.rid);
        key::push_uint(&mut self.key, call.pos);
        key::push_bytes(&mut self.key, call.ref_allele);
        key::push_bytes(&mut self.key, call.alt_allele);
        key::push_uint(&mut self.key, call.end);
        let call_length = u32::try_from(self.key.len()).map_err(|_| describe_length())?;
        key::push_uint(&mut self.key, u8::from(rank.filled));
        key::push_uint(&mut self.key, rank.region);
        key::push_uint(&mut self.key, u64::MAX - rank.tumor_reads);
        key::push_uint(&mut self.key, u64::MAX - rank.reads);
        key::push_uint(&mut self.key, u64::MAX - rank.tumor_depth);
        key::push_uint(&mut self.key, rank.start);
        key::push_bytes(&mut self.key, rank.ref_allele);
        key::push_bytes(&mut self.key, rank.alt_allele);
        self.value.clear();
        self.value.extend_from_slice(&call_length.to_le_bytes());
        self.value.extend_from_slice(line);
        self.sorter.push(&self.key, &self.value).map_err(describe)
    }

    /// Writes the line of the first record of each call in order, checking `stop` between records;
    /// returns how many lines were written, or `None` if `stop` returned true.
    pub fn write_to(
        self,
        out: &mut impl Write,
        stop: &dyn Fn() -> bool,
    ) -> Result<Option<u64>, String> {
        let mut merged: Merged = self.sorter.finish().map_err(describe)?;
        let mut previous: Vec<u8> = Vec::new();
        let mut written = 0;
        while let Some((key, value)) = merged.next_record().map_err(describe)? {
            if stop() {
                return Ok(None);
            }
            let (call_length, line) = value.split_at(4);
            let call_length = u32::from_le_bytes(call_length.try_into().unwrap()) as usize;
            let call = &key[..call_length];
            if previous == call {
                continue;
            }
            previous.clear();
            previous.extend_from_slice(call);
            out.write_all(line)
                .map_err(|error| format!("Could not write the output: {error}!"))?;
            written += 1;
        }
        Ok(Some(written))
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
            filled: false,
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

    fn ranked(records: &[(SortKey, Rank, &str)]) -> String {
        let dir = TempDir::new().unwrap();
        let mut sorter = RecordSorter::new(dir.path()).unwrap();
        for (key, rank, line) in records {
            sorter.push(*key, *rank, line.as_bytes()).unwrap();
        }
        let mut out = Vec::new();
        sorter.write_to(&mut out, &|| false).unwrap();
        String::from_utf8(out).unwrap()
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
    fn test_writes_a_record_vardict_wrote_before_a_filled_one_of_its_call() {
        let call = key(0, 100, "A", "");
        let lines = ranked(&[
            (call, rank(1, 0, 900, 100).filled(), "filled\n"),
            (call, rank(2, 0, 10, 100), "vardict\n"),
        ]);
        assert_eq!(lines, "vardict\n");
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
        let mut out = Vec::new();
        assert_eq!(sorter.write_to(&mut out, &|| false).unwrap(), Some(2048));
        assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn test_stops_when_asked() {
        let dir = TempDir::new().unwrap();
        let mut sorter = RecordSorter::new(dir.path()).unwrap();
        sorter
            .push(key(0, 1, "A", "C"), rank(0, 0, 0, 0), b"line\n")
            .unwrap();
        let mut out = Vec::new();
        assert_eq!(sorter.write_to(&mut out, &|| true).unwrap(), None);
        assert!(out.is_empty());
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
