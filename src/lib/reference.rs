//! Reads bases from an indexed FASTA reference, holding one window of one contig at a time.
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use ahash::AHashMap;
use csv::ReaderBuilder;
use noodles::bgzf;

use crate::fai::{FaiRecord, fai_file};

/// The bases fetched on either side of a request that misses the window held.
const FLANK: u64 = 1 << 15;

/// An indexed FASTA reference whose bases are read uppercase.
pub struct Reference {
    reader: Source,
    entries: AHashMap<String, Entry>,
    contig: String,
    start: u64,
    bases: Vec<u8>,
}

/// A FASTA file, uncompressed or bgzipped and read through its `.gzi` index.
enum Source {
    Plain(BufReader<File>),
    Bgzf(bgzf::io::IndexedReader<File>),
}

impl Reference {
    /// Opens an indexed FASTA file, uncompressed or bgzipped, after checking that its index can be
    /// read and that the file holds every base the index promises. A contig the index lists twice
    /// is read as its first entry describes it.
    pub fn from_path<P: AsRef<Path>>(fasta: P) -> Result<Self, String> {
        let fasta = fasta.as_ref();
        let fai = fai_file(fasta);
        let records = read_fai(&fai)?;
        let bgzipped = check_fasta(fasta, &records)?;
        let describe =
            |error: io::Error| format!("Could not open the reference {fasta:?}: {error}!");
        let reader = if bgzipped {
            bgzf::io::indexed_reader::Builder::default()
                .build_from_path(fasta)
                .map(Source::Bgzf)
        } else {
            File::open(fasta).map(|file| Source::Plain(BufReader::new(file)))
        }
        .map_err(describe)?;
        let mut entries = AHashMap::with_capacity(records.len());
        for record in records {
            entries.entry(record.name.clone()).or_insert(record);
        }
        Ok(Reference {
            reader,
            entries,
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
        self.entry(contig).map(|entry| entry.length)
    }

    /// Returns the length of the longest contig.
    pub fn longest_contig_length(&self) -> u64 {
        self.entries
            .values()
            .map(|entry| entry.length)
            .max()
            .unwrap_or(0)
    }

    /// Returns the index entry of a contig, or an error naming a contig the reference lacks.
    fn entry(&self, contig: &str) -> Result<&Entry, String> {
        self.entries.get(contig).ok_or_else(|| {
            format!("Contig '{contig}' is not in the reference: give --reference the FASTA VarDict was run with!")
        })
    }

    /// Holds the window of a contig around `start..end`.
    fn fetch(&mut self, contig: &str, start: u64, end: u64) -> Result<(), String> {
        let entry = self.entry(contig)?;
        let length = entry.length;
        if end > length {
            return Err(format!(
                "Position {contig}:{end} is beyond the end of {contig}, which has {length} bases in the reference: give --reference the FASTA VarDict was run with!"
            ));
        }
        let from = start.saturating_sub(FLANK);
        let to = (end + FLANK).min(length);
        let describe = |error: io::Error| {
            format!(
                "Could not read {contig}:{}-{to} from the reference: {error}!",
                from + 1
            )
        };
        let position = entry.position(from).map_err(describe)?;
        self.bases.clear();
        match &mut self.reader {
            Source::Plain(reader) => read_bases(reader, position, to - from, &mut self.bases),
            Source::Bgzf(reader) => read_bases(reader, position, to - from, &mut self.bases),
        }
        .map_err(describe)?;
        if self.bases.len() as u64 != to - from {
            return Err(format!(
                "Could not read {contig}:{}-{to} from the reference!",
                from + 1
            ));
        }
        self.contig.clear();
        self.contig.push_str(contig);
        self.start = from;
        self.bases.make_ascii_uppercase();
        Ok(())
    }
}

/// Reads up to `count` bases from a byte position of a FASTA file into `bases`, skipping line
/// breaks and every other byte that is not a printable character, as `samtools faidx` does.
fn read_bases<R: BufRead + Seek>(
    reader: &mut R,
    position: u64,
    count: u64,
    bases: &mut Vec<u8>,
) -> io::Result<()> {
    reader.seek(SeekFrom::Start(position))?;
    let count = count as usize;
    while bases.len() < count {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            break;
        }
        let mut used = 0;
        while used < buffer.len() && bases.len() < count {
            if buffer[used].is_ascii_graphic() {
                bases.push(buffer[used]);
            }
            used += 1;
        }
        reader.consume(used);
    }
    Ok(())
}

/// A contig's entry in a FASTA index.
struct Entry {
    name: String,
    length: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
    end: u64,
}

impl Entry {
    /// The byte position in the FASTA file of the contig's base at a 0-based position.
    fn position(&self, base: u64) -> io::Result<u64> {
        if self.line_bases == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "the FASTA index gives 0 bases per line for contig '{}'",
                    self.name
                ),
            ));
        }
        Ok(self.offset + base / self.line_bases * self.line_width + base % self.line_bases)
    }
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
            offset: fai.offset,
            line_bases: fai.linebases,
            line_width: fai.linewidth,
            end: fai.offset + (last / line_bases) * fai.linewidth + last % line_bases + 1,
        });
    }
    Ok(entries)
}

/// Refuses a FASTA file whose indexed bases cannot all be read: an uncompressed file shorter than
/// its index says, or a compressed file that is not bgzipped or lacks its `.gzi` index. Returns
/// whether the file is bgzipped.
fn check_fasta(fasta: &Path, entries: &[Entry]) -> Result<bool, String> {
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
            (true, true) => Ok(true),
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
        None => Ok(false),
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
        let index = noodles::fasta::fs::index(&path).expect("Cannot index the FASTA file!");
        noodles::fasta::fai::fs::write(fai_file(&path), &index)
            .expect("Cannot index the FASTA file!");
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
        assert_eq!(reference.longest_contig_length(), 12);
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

    #[test]
    fn test_reads_a_bgzipped_reference_across_its_blocks() {
        let (dir, path) = fasta(&[("c1", "ACGTacgtNNac"), ("c2", "ttggCCAA")]);
        let bgzipped = dir.path().join("reference.fa.gz");
        let mut writer = bgzf::io::Writer::new(File::create(&bgzipped).unwrap());
        let mut blocks = Vec::new();
        let text = std::fs::read(&path).unwrap();
        for (index, chunk) in text.chunks(5).enumerate() {
            writer.write_all(chunk).unwrap();
            writer.flush().unwrap();
            blocks.push((writer.position(), 5 * (index as u64 + 1)));
        }
        blocks.pop();
        writer.finish().unwrap();
        let mut gzi = bgzipped.as_os_str().to_owned();
        gzi.push(".gzi");
        bgzf::gzi::fs::write(gzi, &bgzf::gzi::Index::from(blocks)).unwrap();
        std::fs::copy(fai_file(&path), fai_file(&bgzipped)).unwrap();
        let mut reference = Reference::from_path(&bgzipped).unwrap();
        assert_eq!(reference.bases("c1", 2, 10).unwrap(), b"GTACGTNN");
        assert_eq!(reference.bases("c2", 0, 8).unwrap(), b"TTGGCCAA");
    }

    #[test]
    fn test_skips_what_is_not_a_base_as_samtools_faidx_does() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("reference.fa");
        std::fs::write(&path, ">c1\nACGT  \nacgt\t \nAC\n").unwrap();
        std::fs::write(fai_file(&path), "c1\t10\t4\t4\t7\n").unwrap();
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(reference.bases("c1", 2, 10).unwrap(), b"GTACGTAC");
    }

    #[test]
    fn test_reads_a_window_that_starts_past_the_first_line() {
        let sequence: String = (0..40_000u64)
            .map(|i| char::from(b"ACGT"[((i * i + 7 * i) / 3 % 4) as usize]))
            .collect();
        let (_dir, path) = fasta(&[("c1", &sequence)]);
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(
            reference.bases("c1", 39_000, 39_020).unwrap(),
            &sequence.as_bytes()[39_000..39_020]
        );
    }

    #[test]
    fn test_reads_a_contig_the_index_lists_twice_by_its_first_entry() {
        let (_dir, path) = fasta(&[("c1", "ACGT"), ("c2", "TTGGCCAA")]);
        let fai = std::fs::read_to_string(fai_file(&path)).unwrap();
        let second = fai.lines().nth(1).unwrap().replacen("c2", "c1", 1);
        std::fs::write(fai_file(&path), format!("{fai}{second}\n")).unwrap();
        let mut reference = Reference::from_path(&path).unwrap();
        assert_eq!(reference.length("c1").unwrap(), 4);
        assert_eq!(reference.bases("c1", 0, 4).unwrap(), b"ACGT");
    }
}
