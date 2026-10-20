//! Sorts output records by coordinate in bounded memory with `spillsort`, spilling to disk.
use std::io::Write;
use std::path::{Path, PathBuf};

use spillsort::{Config, Merged, Sorter, key};

/// The most memory the records being sorted are held in before they spill to the temporary
/// directory; `spillsort` peaks at about a quarter more.
pub const MAX_BYTES: usize = 64 << 20;

/// Records sorted by contig, in the header's order, then POS, REF, ALT and END. Records with equal
/// keys keep the order they were added in, and only the first of them is written.
pub struct RecordSorter {
    sorter: Sorter,
    key: Vec<u8>,
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

impl RecordSorter {
    /// Makes a sorter that spills to a private directory inside `tmp_dir`.
    pub fn new(tmp_dir: &Path) -> Result<Self, String> {
        let config = Config::default().max_bytes(MAX_BYTES).tmp_dir(tmp_dir);
        let sorter = Sorter::new(config).map_err(describe)?;
        Ok(RecordSorter {
            sorter,
            key: Vec::new(),
        })
    }

    /// Adds a record's VCF line under its key.
    pub fn push(&mut self, key: SortKey, line: &[u8]) -> Result<(), String> {
        self.encode(key);
        self.sorter.push(&self.key, line).map_err(describe)
    }

    /// Writes each record's line in order, skipping a record whose key an earlier record had, and
    /// checking `stop` between records; returns how many lines were written, or `None` if `stop`
    /// returned true.
    pub fn write_to(
        self,
        out: &mut impl Write,
        stop: &dyn Fn() -> bool,
    ) -> Result<Option<u64>, String> {
        let mut merged: Merged = self.sorter.finish().map_err(describe)?;
        let mut previous: Vec<u8> = Vec::new();
        let mut written = 0;
        while let Some((key, line)) = merged.next_record().map_err(describe)? {
            if stop() {
                return Ok(None);
            }
            if previous == key {
                continue;
            }
            previous.clear();
            previous.extend_from_slice(key);
            out.write_all(line)
                .map_err(|error| format!("Could not write the output: {error}!"))?;
            written += 1;
        }
        Ok(Some(written))
    }

    fn encode(&mut self, record: SortKey) {
        self.key.clear();
        key::push_uint(&mut self.key, record.rid);
        key::push_uint(&mut self.key, record.pos);
        key::push_bytes(&mut self.key, record.ref_allele);
        key::push_bytes(&mut self.key, record.alt_allele);
        key::push_uint(&mut self.key, record.end);
    }
}

/// Describes an error from sorting the output records.
fn describe(error: spillsort::Error) -> String {
    format!("Could not sort the output records: {error}!")
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

    fn sorted(records: &[(SortKey, &str)]) -> String {
        let dir = TempDir::new().unwrap();
        let mut sorter = RecordSorter::new(dir.path()).unwrap();
        for (key, line) in records {
            sorter.push(*key, line.as_bytes()).unwrap();
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
            sorter.push(key(0, pos, "A", "C"), &line).unwrap();
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
        sorter.push(key(0, 1, "A", "C"), b"line\n").unwrap();
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
