//! Reads bases from an indexed FASTA reference, holding one window of one contig at a time.
use std::fs::File;
use std::io::Read;
use std::path::Path;

use ahash::AHashMap;
use csv::ReaderBuilder;
use rust_htslib::faidx;

use crate::fai::{FaiRecord, fai_file};

/// The bases fetched on either side of a request that misses the window held.
const FLANK: u64 = 1 << 15;

/// An indexed FASTA reference whose bases are read uppercase.
pub struct Reference {
    reader: faidx::Reader,
    lengths: AHashMap<String, u64>,
    contig: String,
    start: u64,
    bases: Vec<u8>,
}

impl Reference {
    /// Opens an indexed FASTA file, uncompressed or bgzipped, after checking that its index can be
    /// read and that the file holds every base the index promises.
    pub fn from_path<P: AsRef<Path>>(fasta: P) -> Result<Self, String> {
        let fasta = fasta.as_ref();
        let fai = fai_file(fasta);
        let records = read_fai(&fai)?;
        check_fasta(fasta, &records)?;
        let reader = faidx::Reader::from_path(fasta)
            .map_err(|error| format!("Could not open the reference {fasta:?}: {error}!"))?;
        Ok(Reference {
            reader,
            lengths: records
                .into_iter()
                .map(|record| (record.name, record.length))
                .collect(),
            contig: String::new(),
            start: 0,
            bases: Vec::new(),
        })
    }

    /// Returns the uppercase bases of a contig over the 0-based, half-open range `start..end`.
    pub fn bases(&mut self, contig: &str, start: u64, end: u64) -> Result<&[u8], String> {
        if start >= end {
            self.length(contig)?;
            return Ok(&[]);
        }
        let held = self.start..self.start + self.bases.len() as u64;
        if self.contig != contig || start < held.start || end > held.end {
            self.fetch(contig, start, end)?;
        }
        let offset = (start - self.start) as usize;
        Ok(&self.bases[offset..offset + (end - start) as usize])
    }

    /// Returns the uppercase base of a contig at a 0-based position.
    pub fn base(&mut self, contig: &str, position: u64) -> Result<u8, String> {
        Ok(self.bases(contig, position, position + 1)?[0])
    }

    /// Returns the length of a contig, or an error naming a contig the reference lacks.
    pub fn length(&self, contig: &str) -> Result<u64, String> {
        self.lengths.get(contig).copied().ok_or_else(|| {
            format!("Contig '{contig}' is not in the reference: give --reference the FASTA VarDict was run with!")
        })
    }

    /// Holds the window of a contig around `start..end`.
    fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<(), String> {
        let length = self.length(contig)?;
        if end > length {
            return Err(format!(
                "Position {contig}:{end} is beyond the end of {contig}, which has {length} bases in the reference: give --reference the FASTA VarDict was run with!"
            ));
        }
        let from = start.saturating_sub(FLANK);
        let to = (end + FLANK).min(length);
        let bases = self
            .reader
            .fetch_seq(contig, from as usize, (to - 1) as usize)
            .map_err(|error| {
                format!(
                    "Could not read {contig}:{}-{to} from the reference: {error}!",
                    from + 1
                )
            })?;
        if bases.len() as u64 != to - from {
            return Err(format!(
                "Could not read {contig}:{}-{to} from the reference!",
                from + 1
            ));
        }
        self.contig.clear();
        self.contig.push_str(contig);
        self.start = from;
        self.bases = bases;
        self.bases.make_ascii_uppercase();
        Ok(())
    }
}

/// A contig's entry in a FASTA index.
struct Entry {
    name: String,
    length: u64,
    end: u64,
}

/// Reads the contigs of a FASTA index with the byte offset just past each one's last base.
fn read_fai(fai: &Path) -> Result<Vec<Entry>, String> {
    let describe = |error: csv::Error| format!("Could not read the FASTA index {fai:?}: {error}!");
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(fai)
        .map_err(describe)?;
    let mut entries = Vec::new();
    for record in reader.records() {
        let record = record.map_err(describe)?;
        let fai: FaiRecord = record.deserialize(None).map_err(describe)?;
        let last = fai.length.saturating_sub(1);
        let line_bases = fai.linebases.max(1);
        entries.push(Entry {
            name: fai.name.to_string(),
            length: fai.length,
            end: fai.offset + (last / line_bases) * fai.linewidth + last % line_bases + 1,
        });
    }
    Ok(entries)
}

/// Refuses a FASTA file htslib could not read every indexed base from: an uncompressed file shorter
/// than its index says, or a compressed file that is not bgzipped or lacks its `.gzi` index.
fn check_fasta(fasta: &Path, entries: &[Entry]) -> Result<(), String> {
    let mut file = File::open(fasta)
        .map_err(|error| format!("Could not open the reference {fasta:?}: {error}!"))?;
    let mut magic = [0u8; 14];
    let read = file.read(&mut magic).unwrap_or(0);
    if read >= 2 && magic[..2] == [0x1f, 0x8b] {
        let mut gzi = fasta.as_os_str().to_owned();
        gzi.push(".gzi");
        return match (
            read == magic.len() && magic[12..14] == *b"BC",
            Path::new(&gzi).exists(),
        ) {
            (true, true) => Ok(()),
            (true, false) => Err(format!(
                "Could not find the .gzi index of the bgzipped reference {fasta:?}!"
            )),
            (false, _) => Err(format!(
                "The reference {fasta:?} is compressed but not with bgzip!"
            )),
        };
    }
    let size = file
        .metadata()
        .map_err(|error| format!("Could not open the reference {fasta:?}: {error}!"))?
        .len();
    match entries.iter().find(|entry| entry.end > size) {
        Some(entry) => Err(format!(
            "The reference {fasta:?} ends before contig '{}' does in its index: rebuild the index with samtools faidx!",
            entry.name
        )),
        None => Ok(()),
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::io::Write;

    use pretty_assertions::assert_eq;
    use std::path::PathBuf;

    use tempfile::TempDir;

    use super::*;

    /// Writes contigs to an indexed FASTA file in a new directory, wrapping lines at 7 bases.
    pub(crate) fn fasta(contigs: &[(&str, &str)]) -> (TempDir, PathBuf) {
        let dir = TempDir::new().expect("Cannot create a temporary directory!");
        let path = dir.path().join("reference.fa");
        let mut file = File::create(&path).expect("Cannot create the FASTA file!");
        for (name, sequence) in contigs {
            writeln!(file, ">{name}").unwrap();
            for line in sequence.as_bytes().chunks(7) {
                file.write_all(line).unwrap();
                writeln!(file).unwrap();
            }
        }
        drop(file);
        faidx::build(&path).expect("Cannot index the FASTA file!");
        (dir, path)
    }

    #[test]
    fn test_reads_bases_uppercase_across_lines_and_contigs() {
        let (_dir, path) = fasta(&[("c1", "ACGTacgtNNac"), ("c2", "ttggCCAA")]);
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(reference.bases("c1", 2, 10).unwrap(), b"GTACGTNN");
        assert_eq!(reference.base("c2", 0).unwrap(), b'T');
        assert_eq!(reference.bases("c1", 0, 12).unwrap(), b"ACGTACGTNNAC");
        assert_eq!(reference.length("c2").unwrap(), 8);
        assert_eq!(reference.bases("c1", 4, 4).unwrap(), b"");
        assert_eq!(reference.bases("c2", 4, 4).unwrap(), b"");
        assert!(reference.bases("c9", 4, 4).is_err());
    }

    #[test]
    fn test_refuses_a_contig_it_lacks_and_a_position_past_the_end() {
        let (_dir, path) = fasta(&[("c1", "ACGT")]);
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(
            reference.base("c9", 0).unwrap_err(),
            "Contig 'c9' is not in the reference: give --reference the FASTA VarDict was run with!"
        );
        assert_eq!(
            reference.bases("c1", 2, 5).unwrap_err(),
            "Position c1:5 is beyond the end of c1, which has 4 bases in the reference: give --reference the FASTA VarDict was run with!"
        );
    }

    #[test]
    fn test_refuses_a_reference_without_an_index() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("reference.fa");
        std::fs::write(&path, ">c1\nACGT\n").unwrap();
        assert!(
            Reference::from_path(&path)
                .err()
                .unwrap()
                .starts_with("Could not read the FASTA index")
        );
    }

    #[test]
    fn test_refuses_a_reference_missing_its_sequence() {
        let (_dir, path) = fasta(&[("c1", "ACGTACGTACGT")]);
        std::fs::write(&path, ">c1\nACGTACG\n").unwrap();
        assert!(Reference::from_path(&path).err().unwrap().ends_with(
            "ends before contig 'c1' does in its index: rebuild the index with samtools faidx!"
        ));
        std::fs::remove_file(&path).unwrap();
        assert!(
            Reference::from_path(&path)
                .err()
                .unwrap()
                .starts_with("Could not open the reference")
        );
    }

    #[test]
    fn test_reads_a_reference_without_a_final_newline() {
        let (_dir, path) = fasta(&[("c1", "ACGTACGTACGTAC")]);
        std::fs::write(&path, ">c1\nACGTACG\nTACGTAC").unwrap();
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(reference.bases("c1", 10, 14).unwrap(), b"GTAC");
    }

    #[test]
    fn test_refuses_a_gzipped_reference_without_bgzip_framing() {
        let (_dir, path) = fasta(&[("c1", "ACGT")]);
        std::fs::write(&path, [0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 3]).unwrap();
        assert!(
            Reference::from_path(&path)
                .err()
                .unwrap()
                .ends_with("is compressed but not with bgzip!")
        );
    }
}
