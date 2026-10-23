//! The VCF output: standard output, a plain file, or a bgzipped and indexed file for a `.gz` path.
use std::fs::File;
use std::io::{self, BufRead, BufWriter, ErrorKind, Read, Seek, SeekFrom, Stdout, Write};
use std::path::{Path, PathBuf};

use ahash::AHashMap;
use noodles::bgzf;
use noodles::core::Position;
use noodles::csi::binning_index::index::ReferenceSequence;
use noodles::csi::binning_index::index::header::Builder as IndexHeader;
use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles::csi::binning_index::index::reference_sequence::index::{BinnedIndex, LinearIndex};
use noodles::csi::binning_index::{self, Indexer, ReferenceSequence as _};
use noodles::csi::{self, BinningIndex};
use noodles::tabix;

use crate::io::has_gzip_ext;

/// The empty block every complete BGZF file ends with.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// The index written beside a bgzipped regular file, as `tabix -p vcf` writes it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Index {
    /// A tabix index, `<output>.tbi`.
    Tbi,
    /// A CSI index, `<output>.csi`, for a reference with a contig too long for a tabix index.
    Csi,
}

impl Index {
    /// The longest contig a tabix index can address.
    pub const MAX_TBI_CONTIG_LENGTH: u64 = 1 << 29;

    /// The width of the smallest bin of either index, as a power of 2.
    const MIN_SHIFT: u8 = 14;

    /// The levels of bins in a tabix index.
    const TBI_DEPTH: u8 = 5;

    /// The fewest levels of bins in a CSI index.
    const MIN_CSI_DEPTH: u8 = 6;

    /// The index that can address every contig of a reference whose longest contig has this length.
    pub fn for_longest_contig(length: u64) -> Self {
        if length > Self::MAX_TBI_CONTIG_LENGTH {
            Index::Csi
        } else {
            Index::Tbi
        }
    }

    /// The path of this index beside a bgzipped file.
    fn path(self, output: &Path) -> PathBuf {
        let mut path = output.as_os_str().to_owned();
        path.push(match self {
            Index::Tbi => ".tbi",
            Index::Csi => ".csi",
        });
        PathBuf::from(path)
    }

    /// Writes this index of a bgzipped VCF file beside it and syncs it to disk.
    fn write(self, output: &Path) -> Result<(), String> {
        let path = self.path(output);
        index_vcf(output, self)
            .and_then(|index| match self {
                Index::Tbi => tabix::fs::write(&path, &index),
                Index::Csi => csi::fs::write(&path, &to_csi(&index)),
            })
            .and_then(|()| File::open(&path).and_then(|file| file.sync_all()))
            .map_err(|error| format!("Could not write the index {path:?}: {error}!"))
    }

    /// The levels of bins this index has for a VCF whose `##contig` lines give this longest length,
    /// as `tabix` chooses them: a CSI index adds levels until its bins span the contig and 256 bases.
    fn depth(self, longest_contig: u64) -> u8 {
        match self {
            Index::Tbi => Self::TBI_DEPTH,
            Index::Csi => {
                let mut depth = Self::MIN_CSI_DEPTH;
                while longest_contig.saturating_add(256) > span(depth) {
                    depth += 1;
                }
                depth
            }
        }
    }
}

/// The bases the bins of an index of this depth span.
fn span(depth: u8) -> u64 {
    1u64.checked_shl(u32::from(Index::MIN_SHIFT + 3 * depth))
        .unwrap_or(u64::MAX)
}

/// Indexes a bgzipped VCF file as `tabix -p vcf` does, with bins of the depth the kind of index
/// has; [`to_csi`] makes a CSI index of it. Each line not starting with `#` is a record on its
/// CHROM from its POS to the end of its REF allele, or to its INFO END when that is later; a line
/// without a numeric POS is skipped. A CSI depth comes from the longest `##contig` length, or from
/// 100 Gbp without one. Refuses records of a contig that are not together, or not sorted by POS,
/// and records beyond what the index can address.
fn index_vcf(path: &Path, kind: Index) -> io::Result<tabix::Index> {
    let mut reader = bgzf::io::Reader::new(File::open(path)?);
    let mut line = Vec::new();
    let mut longest_contig = 0;
    let mut indexer: Option<(Indexer<LinearIndex>, u64)> = None;
    let mut names: Vec<String> = Vec::new();
    let mut ids: AHashMap<String, usize> = AHashMap::new();
    let mut last: Option<(usize, u64)> = None;
    let mut chunk_start = reader.virtual_position();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        let chunk_end = reader.virtual_position();
        if line.ends_with(b"\n") {
            line.pop();
        }
        if line.ends_with(b"\r") {
            line.pop();
        }
        if line.first() == Some(&b'#') {
            if indexer.is_none() {
                longest_contig = longest_contig.max(contig_length(&line));
                chunk_start = chunk_end;
            }
            continue;
        }
        let (indexer, max_position) = indexer.get_or_insert_with(|| {
            let longest = if longest_contig == 0 {
                100 << 30
            } else {
                longest_contig
            };
            let depth = kind.depth(longest);
            (Indexer::new(Index::MIN_SHIFT, depth), span(depth))
        });
        let Some((name, start, end)) = parse_interval(&line) else {
            continue;
        };
        if start > *max_position || end > *max_position {
            return Err(io::Error::other(format!(
                "{name}:{} ends at {end}, beyond the {max_position} bases this index can address",
                start + 1
            )));
        }
        let id = match ids.get(name) {
            Some(&id) => id,
            None => {
                names.push(name.to_string());
                ids.insert(name.to_string(), names.len() - 1);
                names.len() - 1
            }
        };
        match last {
            Some((last_id, _)) if last_id != id && id < last_id => {
                return Err(io::Error::other(format!(
                    "the records of {name} are not together, as {} comes between them",
                    names[last_id]
                )));
            }
            Some((last_id, last_start)) if last_id == id && last_start > start => {
                return Err(io::Error::other(format!(
                    "the records of {name} are not sorted: {} comes before {}",
                    last_start + 1,
                    start + 1
                )));
            }
            _ => {}
        }
        if end < start {
            return Err(io::Error::other(format!(
                "{name}:{} ends at {end}, before it starts",
                start + 1
            )));
        }
        last = Some((id, start));
        let position =
            |position: u64| Position::try_from(position as usize).map_err(io::Error::other);
        indexer.add_record(
            Some((id, position(start + 1)?, position(end)?, true)),
            Chunk::new(chunk_start, chunk_end),
        )?;
        chunk_start = chunk_end;
    }
    let indexer = match indexer {
        Some((indexer, _)) => indexer,
        None => Indexer::new(Index::MIN_SHIFT, kind.depth(0)),
    };
    let count = names.len();
    let header = IndexHeader::vcf()
        .set_reference_sequence_names(names.into_iter().map(Into::into).collect())
        .build();
    Ok(indexer.set_header(header).build(count))
}

/// The CSI index of an index built with a CSI depth. Each bin's offset is the linear index's at the
/// bin's first window, as htslib sets it: a reader skips the chunks that end before the offset of
/// the bin a query starts in, so the offset must not pass a longer record, of a larger bin, that
/// overlaps the window.
fn to_csi(index: &tabix::Index) -> csi::Index {
    let depth = index.depth();
    let reference_sequences = index
        .reference_sequences()
        .iter()
        .map(|reference_sequence| {
            let windows = reference_sequence.index();
            let offsets: BinnedIndex = reference_sequence
                .bins()
                .keys()
                .map(|&bin| {
                    let offset = windows.get(first_window(bin, depth)).copied();
                    (bin, offset.unwrap_or_default())
                })
                .collect();
            ReferenceSequence::new(
                reference_sequence.bins().clone(),
                offsets,
                reference_sequence.metadata().cloned(),
            )
        })
        .collect();
    let mut builder = binning_index::Index::builder()
        .set_min_shift(index.min_shift())
        .set_depth(depth)
        .set_reference_sequences(reference_sequences);
    if let Some(header) = index.header() {
        builder = builder.set_header(header.clone());
    }
    if let Some(count) = index.unplaced_unmapped_record_count() {
        builder = builder.set_unplaced_unmapped_record_count(count);
    }
    builder.build()
}

/// The first window of the smallest bins that a bin of an index of this depth spans.
fn first_window(bin: usize, depth: u8) -> usize {
    let mut level = 0;
    let mut ancestor = bin;
    while ancestor > 0 {
        ancestor = (ancestor - 1) >> 3;
        level += 1;
    }
    let first_of_level = ((1 << (3 * level)) - 1) / 7;
    (bin - first_of_level) << (3 * (usize::from(depth) - level))
}

/// The length a `##contig` line gives, or 0 for any other line.
fn contig_length(line: &[u8]) -> u64 {
    let Some(rest) = line.strip_prefix(b"##contig") else {
        return 0;
    };
    let Some(at) = rest.windows(6).position(|window| window == b"length") else {
        return 0;
    };
    let value = &rest[at + 6..];
    let value = &value[value
        .iter()
        .take_while(|&&byte| byte == b' ' || byte == b'=')
        .count()..];
    leading_integer(value).map_or(0, |length| length.max(0) as u64)
}

/// The CHROM of a record line and the 0-based, half-open interval tabix indexes it over.
fn parse_interval(line: &[u8]) -> Option<(&str, u64, u64)> {
    let mut fields = line.split(|&byte| byte == b'\t');
    let name = std::str::from_utf8(fields.next()?).ok()?;
    let start = leading_integer(fields.next()?)?.saturating_sub(1).max(0) as u64;
    let mut end = 1;
    if let Some(ref_allele) = fields.nth(1).filter(|field| !field.is_empty()) {
        end = start + ref_allele.len() as u64;
    }
    if let Some(info_end) = fields.nth(3).and_then(info_end)
        && info_end > start as i64
    {
        end = info_end as u64;
    }
    Some((name, start, end))
}

/// The value of an INFO field's END key, read as tabix reads it: 0 when it holds no number.
fn info_end(info: &[u8]) -> Option<i64> {
    let value = if let Some(value) = info.strip_prefix(b"END=") {
        value
    } else {
        let at = info.windows(5).position(|window| window == b";END=")?;
        &info[at + 5..]
    };
    if value.first() == Some(&b'.') {
        return None;
    }
    Some(leading_integer(value).unwrap_or(0))
}

/// The decimal integer a field starts with, with an optional sign, saturating at the limits of an
/// `i64`; `None` when the field does not start with a digit.
fn leading_integer(field: &[u8]) -> Option<i64> {
    let (negative, digits) = match field.first() {
        Some(b'-') => (true, &field[1..]),
        Some(b'+') => (false, &field[1..]),
        _ => (false, field),
    };
    let digits = &digits[..digits
        .iter()
        .take_while(|byte| byte.is_ascii_digit())
        .count()];
    if digits.is_empty() {
        return None;
    }
    let magnitude = digits.iter().fold(0i64, |value, digit| {
        value
            .saturating_mul(10)
            .saturating_add(i64::from(digit - b'0'))
    });
    Some(if negative { -magnitude } else { magnitude })
}

/// Removes any tabix or CSI index beside a bgzipped file.
fn remove_indexes(output: &Path) -> Result<(), String> {
    for index in [Index::Tbi, Index::Csi] {
        let path = index.path(output);
        match std::fs::remove_file(&path) {
            Err(error) if error.kind() != ErrorKind::NotFound => {
                return Err(format!(
                    "Could not remove the stale index {path:?}: {error}!"
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

/// Where the VCF text goes. A file is created at once and, if it is a regular file, removed again
/// unless the run finishes; a device, pipe or FIFO given as the output is only ever closed. A
/// bgzipped regular file is indexed once it is complete, and its index removed along with it.
pub enum Output {
    /// Standard output, for no path or `-`, which notes when its reader has gone.
    Stdout {
        /// The buffered standard output.
        writer: BufWriter<Stdout>,
        /// Whether a write found the pipe closed.
        closed: bool,
    },
    /// A plain text file, and whether it is a regular file.
    Plain(PathBuf, BufWriter<File>, bool),
    /// A bgzipped file, and whether it is a regular file.
    Bgzf(PathBuf, bgzf::io::Writer<File>, bool),
}

impl Output {
    /// Creates the output at a path, bgzipped when it ends in `.gz`, or standard output.
    pub fn create(path: Option<&Path>) -> Result<Self, String> {
        let describe = |error: &dyn std::fmt::Display| {
            format!("Could not create the output {path:?}: {error}!")
        };
        match path {
            Some(path) if path.as_os_str() != "-" && has_gzip_ext(&path) => {
                let writer = File::create(path)
                    .map(bgzf::io::Writer::new)
                    .map_err(|error| describe(&error))?;
                let regular = std::fs::metadata(path).is_ok_and(|metadata| metadata.is_file());
                if regular && let Err(error) = remove_indexes(path) {
                    drop(writer);
                    let _ = std::fs::remove_file(path);
                    return Err(error);
                }
                Ok(Output::Bgzf(path.to_path_buf(), writer, regular))
            }
            Some(path) if path.as_os_str() != "-" => {
                let file = File::create(path).map_err(|error| describe(&error))?;
                let regular = file.metadata().is_ok_and(|metadata| metadata.is_file());
                Ok(Output::Plain(
                    path.to_path_buf(),
                    BufWriter::new(file),
                    regular,
                ))
            }
            _ => Ok(Output::Stdout {
                writer: BufWriter::new(std::io::stdout()),
                closed: false,
            }),
        }
    }

    /// Whether standard output's reader closed it, as `head` does once it has read enough.
    pub fn is_closed_pipe(&self) -> bool {
        matches!(self, Output::Stdout { closed: true, .. })
    }

    /// Flushes everything written and closes the output, checking each step, and removes a regular
    /// file that could not be completed. A bgzipped regular file is read back once closed to check
    /// that it ends with the BGZF end-of-file block, and is then indexed. Returns false, after
    /// removing the file and its index, when `stop` asks the run to stop once the index is written.
    pub fn finish(self, index: Index, stop: &dyn Fn() -> bool) -> Result<bool, String> {
        let describe =
            |error: &dyn std::fmt::Display| format!("Could not write the output: {error}!");
        match self {
            Output::Stdout { mut writer, .. } => match writer.flush() {
                Err(error) if error.kind() != ErrorKind::BrokenPipe => Err(describe(&error)),
                _ => Ok(true),
            },
            Output::Plain(path, writer, regular) => {
                let result = writer
                    .into_inner()
                    .map_err(|error| error.into_error())
                    .and_then(|file| if regular { file.sync_all() } else { Ok(()) });
                if result.is_err() && regular {
                    let _ = std::fs::remove_file(path);
                }
                result.map(|()| true).map_err(|error| describe(&error))
            }
            Output::Bgzf(path, writer, regular) => {
                let flushed = writer.finish().map(drop);
                if !regular {
                    return flushed.map(|()| true).map_err(|error| describe(&error));
                }
                let result = flushed
                    .and_then(|()| check_bgzf_end(&path))
                    .map_err(|error| describe(&error))
                    .and_then(|()| index.write(&path));
                let stopped = result.is_ok() && stop();
                if result.is_err() || stopped {
                    let _ = std::fs::remove_file(&path);
                    let _ = remove_indexes(&path);
                }
                result.map(|()| !stopped)
            }
        }
    }

    /// Closes the output without writing what it still buffers, and removes it if it is a regular
    /// file.
    pub fn discard(self) {
        let path = match self {
            Output::Stdout { writer, .. } => {
                let _ = writer.into_parts();
                return;
            }
            Output::Plain(path, writer, regular) => {
                let _ = writer.into_parts();
                regular.then_some(path)
            }
            Output::Bgzf(path, writer, regular) => {
                drop(writer);
                regular.then_some(path)
            }
        };
        if let Some(path) = path {
            let _ = std::fs::remove_file(path);
        }
    }
}

/// Checks that a closed bgzipped file ends with the BGZF end-of-file block and syncs it to disk.
fn check_bgzf_end(path: &Path) -> std::io::Result<()> {
    let mut file = File::open(path)?;
    let mut end = [0u8; BGZF_EOF.len()];
    let complete = file.seek(SeekFrom::End(-(BGZF_EOF.len() as i64))).is_ok()
        && file.read_exact(&mut end).is_ok()
        && end == BGZF_EOF;
    if !complete {
        return Err(std::io::Error::other(
            "the bgzipped file is incomplete, as the disk or a size limit may be full",
        ));
    }
    file.sync_all()
}

impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Output::Stdout { writer, closed } => writer.write(buf).inspect_err(|error| {
                *closed |= error.kind() == ErrorKind::BrokenPipe;
            }),
            Output::Plain(_, writer, _) => writer.write(buf),
            Output::Bgzf(_, writer, _) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Output::Stdout { writer, closed } => writer.flush().inspect_err(|error| {
                *closed |= error.kind() == ErrorKind::BrokenPipe;
            }),
            Output::Plain(_, writer, _) => writer.flush(),
            Output::Bgzf(_, writer, _) => writer.flush(),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use tempfile::TempDir;

    use super::*;

    /// A reader of the text of a bgzipped file.
    fn bgzf_reader(path: impl AsRef<Path>) -> bgzf::io::Reader<File> {
        File::open(path).map(bgzf::io::Reader::new).unwrap()
    }

    /// Whether a file starts with a BGZF block header, as htslib's `bgzf_is_bgzf` checks.
    fn is_bgzip(path: &Path) -> io::Result<bool> {
        let mut start = [0u8; 16];
        File::open(path)?.read_exact(&mut start)?;
        Ok(start == BGZF_EOF[..16])
    }

    /// A VCF of one contig `c` of this length, with a record at each 1-based position in turn.
    fn vcf(length: u64, positions: &[u64]) -> String {
        let mut text = format!("##fileformat=VCFv4.2\n##contig=<ID=c,length={length}>\n");
        text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        for position in positions {
            text.push_str(&format!("c\t{position}\t.\tA\tC\t.\t.\t.\n"));
        }
        text
    }

    /// Writes text to an output at a path and finishes it.
    fn write(
        path: &Path,
        text: &str,
        index: Index,
        stop: &dyn Fn() -> bool,
    ) -> Result<bool, String> {
        let mut output = Output::create(Some(path)).unwrap();
        output.write_all(text.as_bytes()).unwrap();
        output.finish(index, stop)
    }

    /// The positions of the records an index of a bgzipped VCF finds over a 1-based region of `c`.
    fn query(path: &Path, start: u64, end: u64) -> Vec<String> {
        let csi = Index::Csi.path(path);
        let chunks = if csi.exists() {
            bin_chunks(&csi::fs::read(csi).unwrap(), start, end)
        } else {
            bin_chunks(&tabix::fs::read(Index::Tbi.path(path)).unwrap(), start, end)
        };
        let mut reader = bgzf_reader(path);
        let mut positions = Vec::new();
        for chunk in chunks {
            reader.seek(chunk.start()).unwrap();
            while reader.virtual_position() < chunk.end() {
                let mut line = String::new();
                reader.read_line(&mut line).unwrap();
                let fields: Vec<&str> = line.split('\t').collect();
                let position: u64 = fields[1].parse().unwrap();
                if position <= end && position + fields[3].len() as u64 > start {
                    positions.push(fields[1].to_string());
                }
            }
        }
        positions
    }

    /// The chunks of the bins of the first contig that overlap a 1-based region, in file order;
    /// noodles refuses to query the last base a tabix index can address, so the bins are found here.
    fn bin_chunks<I: noodles::csi::binning_index::index::reference_sequence::Index>(
        index: &binning_index::Index<I>,
        start: u64,
        end: u64,
    ) -> Vec<Chunk> {
        let (min_shift, depth) = (u32::from(index.min_shift()), u32::from(index.depth()));
        let (start, end) = ((start - 1) as usize, (end - 1) as usize);
        let bins = index.reference_sequences()[0].bins();
        let mut chunks: Vec<Chunk> = (0..=depth)
            .flat_map(|level| {
                let shift = min_shift + 3 * (depth - level);
                let first = ((1usize << (3 * level)) - 1) / 7;
                (first + (start >> shift))..=(first + (end >> shift))
            })
            .filter_map(|id| bins.get(&id))
            .flat_map(|bin| bin.chunks().iter().copied())
            .collect();
        chunks.sort_by_key(|chunk| chunk.start());
        chunks
    }

    /// The names of the files in a directory, sorted.
    fn files(dir: &TempDir) -> Vec<String> {
        let mut names: Vec<String> = std::fs::read_dir(dir.path())
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        names.sort();
        names
    }

    #[test]
    fn test_writes_a_plain_file() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf");
        assert!(write(&path, "line\n", Index::Tbi, &|| false).unwrap());
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "line\n");
        assert_eq!(files(&dir), ["out.vcf"]);
    }

    #[test]
    fn test_writes_a_bgzipped_file_for_a_gz_path() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let text = vcf(100, &[10]);
        assert!(write(&path, &text, Index::Tbi, &|| false).unwrap());
        assert!(is_bgzip(&path).unwrap());
        let mut written = String::new();
        bgzf_reader(&path).read_to_string(&mut written).unwrap();
        assert_eq!(written, text);
    }

    #[rstest]
    #[case(Index::Tbi, "out.vcf.gz.tbi", b"TBI\x01")]
    #[case(Index::Csi, "out.vcf.gz.csi", b"CSI\x01")]
    fn test_indexes_a_bgzipped_file_in_place_of_any_stale_index(
        #[case] index: Index,
        #[case] name: &str,
        #[case] magic: &[u8; 4],
    ) {
        let dir = TempDir::new().unwrap();
        for stale in ["out.vcf.gz.tbi", "out.vcf.gz.csi"] {
            std::fs::write(dir.path().join(stale), "stale").unwrap();
        }
        let path = dir.path().join("out.vcf.gz");
        assert!(write(&path, &vcf(100, &[10, 20, 30]), index, &|| false).unwrap());
        assert_eq!(files(&dir), ["out.vcf.gz", name]);
        let mut start = [0u8; 4];
        bgzf_reader(dir.path().join(name))
            .read_exact(&mut start)
            .unwrap();
        assert_eq!(&start, magic);
        assert_eq!(query(&path, 15, 30), ["20", "30"]);
    }

    #[rstest]
    #[case(1 << 29, Index::Tbi)]
    #[case((1 << 29) + 1, Index::Csi)]
    fn test_indexes_the_last_base_of_the_longest_contig_it_chooses_for(
        #[case] length: u64,
        #[case] index: Index,
    ) {
        assert_eq!(Index::for_longest_contig(length), index);
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let text = vcf(length, &[length]);
        assert!(write(&path, &text, index, &|| false).unwrap());
        assert_eq!(query(&path, length, length), [length.to_string()]);
        if index == Index::Csi {
            assert!(write(&path, &text, Index::Tbi, &|| false).is_err());
        }
    }

    #[test]
    fn test_removes_the_file_and_its_index_when_stopped_once_indexed() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let index = dir.path().join("out.vcf.gz.tbi");
        assert!(!write(&path, &vcf(100, &[10]), Index::Tbi, &|| index.exists()).unwrap());
        assert!(files(&dir).is_empty());
    }

    #[test]
    fn test_removes_the_file_when_it_cannot_be_indexed() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let error = write(&path, &vcf(100, &[20, 10]), Index::Tbi, &|| false).unwrap_err();
        assert!(error.starts_with("Could not write the index"), "{error}");
        assert!(files(&dir).is_empty());
    }

    #[test]
    fn test_refuses_a_stale_index_it_cannot_remove() {
        let dir = TempDir::new().unwrap();
        let stale = dir.path().join("out.vcf.gz.csi");
        std::fs::create_dir(&stale).unwrap();
        std::fs::write(stale.join("file"), "").unwrap();
        let error = Output::create(Some(&dir.path().join("out.vcf.gz")))
            .err()
            .unwrap();
        assert!(
            error.starts_with("Could not remove the stale index"),
            "{error}"
        );
        assert_eq!(files(&dir), ["out.vcf.gz.csi"]);
    }

    #[test]
    fn test_refuses_a_bgzipped_file_without_its_end_of_file_block() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        std::fs::write(&path, &BGZF_EOF[..20]).unwrap();
        assert!(check_bgzf_end(&path).is_err());
        std::fs::write(&path, BGZF_EOF).unwrap();
        assert!(check_bgzf_end(&path).is_ok());
    }

    #[test]
    fn test_discarding_removes_the_file() {
        let dir = TempDir::new().unwrap();
        for name in ["out.vcf", "out.vcf.gz"] {
            let path = dir.path().join(name);
            let mut output = Output::create(Some(&path)).unwrap();
            output.write_all(b"line\n").unwrap();
            output.discard();
            assert!(!path.exists());
        }
    }

    #[test]
    fn test_refuses_a_path_it_cannot_create() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("missing").join("out.vcf");
        assert!(
            Output::create(Some(&path))
                .err()
                .unwrap()
                .starts_with("Could not create the output")
        );
    }

    #[rstest]
    #[case(Index::Tbi, 0, 5)]
    #[case(Index::Tbi, 1 << 40, 5)]
    #[case(Index::Csi, 0, 6)]
    #[case(Index::Csi, (1 << 32) - 256, 6)]
    #[case(Index::Csi, (1 << 32) - 255, 7)]
    #[case(Index::Csi, 100 << 30, 8)]
    fn test_index_depth_spans_the_longest_contig_as_tabix_chooses_it(
        #[case] index: Index,
        #[case] longest_contig: u64,
        #[case] depth: u8,
    ) {
        assert_eq!(index.depth(longest_contig), depth);
    }

    #[rstest]
    #[case(b"c\t10\t.\tA\tC\t.\t.\t.", Some(("c", 9, 10)))]
    #[case(b"c\t10\t.\tACGT\tA\t.\t.\t.", Some(("c", 9, 13)))]
    #[case(b"c\t10\t.\tA\t<DEL>\t.\t.\tTYPE=DEL;END=40", Some(("c", 9, 40)))]
    #[case(b"c\t10\t.\tA\t<DEL>\t.\t.\tEND=40;TYPE=DEL", Some(("c", 9, 40)))]
    #[case(b"c\t10\t.\tA\t<DEL>\t.\t.\tEND=5", Some(("c", 9, 10)))]
    #[case(b"c\t10\t.\tA\t<DEL>\t.\t.\tEND=.", Some(("c", 9, 10)))]
    #[case(b"c\t10\t.\tA\t<DEL>\t.\t.\tSVEND=40", Some(("c", 9, 10)))]
    #[case(b"c\t0\t.\tA\tC\t.\t.\t.", Some(("c", 0, 1)))]
    #[case(b"c\t10\t.\t\tC\t.\t.\t.", Some(("c", 9, 1)))]
    #[case(b"c\tx\t.\tA\tC\t.\t.\t.", None)]
    #[case(b"c", None)]
    fn test_parse_interval_reads_a_record_as_tabix_does(
        #[case] line: &[u8],
        #[case] expected: Option<(&str, u64, u64)>,
    ) {
        assert_eq!(parse_interval(line), expected);
    }

    #[rstest]
    #[case(
        "c1\t20\t.\tA\tC\t.\t.\t.\nc1\t10\t.\tA\tC\t.\t.\t.\n",
        "the records of c1 are not sorted: 20 comes before 10"
    )]
    #[case(
        "c1\t10\t.\tA\tC\t.\t.\t.\nc2\t10\t.\tA\tC\t.\t.\t.\nc1\t20\t.\tA\tC\t.\t.\t.\n",
        "the records of c1 are not together, as c2 comes between them"
    )]
    #[case("c1\t10\t.\t\tC\t.\t.\t.\n", "c1:10 ends at 1, before it starts")]
    fn test_refuses_records_tabix_cannot_index(#[case] records: &str, #[case] message: &str) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let text = format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{records}"
        );
        let error = write(&path, &text, Index::Tbi, &|| false).unwrap_err();
        assert_eq!(
            error,
            format!(
                "Could not write the index {:?}: {message}!",
                Index::Tbi.path(&path)
            )
        );
    }

    #[test]
    fn test_a_csi_bin_offset_does_not_pass_a_longer_record_overlapping_the_bin() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let mut text = vcf(300_000, &[]);
        text.push_str("c\t10\t.\tA\t<DEL>\t.\t.\tEND=300000\nc\t20000\t.\tA\tC\t.\t.\t.\n");
        assert!(write(&path, &text, Index::Csi, &|| false).unwrap());
        let index = csi::fs::read(Index::Csi.path(&path)).unwrap();
        let reference_sequence = &index.reference_sequences()[0];
        let first_record = reference_sequence.metadata().unwrap().start_position();
        let depth = u32::from(index.depth());
        let bin = ((1usize << (3 * depth)) - 1) / 7 + (20000 - 1) / (1 << 14);
        assert!(reference_sequence.bins().contains_key(&bin));
        assert_eq!(reference_sequence.index().get(&bin), Some(&first_record));
    }

    #[test]
    fn test_a_csi_bin_above_the_leaves_takes_the_offset_of_its_first_window() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let mut text = vcf(300_000, &[10, 20_000, 135_000]);
        text.push_str("c\t140000\t.\tA\t<DEL>\t.\t.\tEND=150000\nc\t160000\t.\tA\tC\t.\t.\t.\n");
        assert!(write(&path, &text, Index::Csi, &|| false).unwrap());
        let index = csi::fs::read(Index::Csi.path(&path)).unwrap();
        let level = u32::from(index.depth()) - 1;
        let bin = ((1usize << (3 * level)) - 1) / 7 + (140_000 - 1) / (1 << (14 + 3));
        // The text fits one BGZF block, so a record's virtual position is its offset in the text.
        let first_record_in_window = text.find("c\t135000\t").unwrap() as u64;
        let reference_sequence = &index.reference_sequences()[0];
        assert!(reference_sequence.bins().contains_key(&bin));
        assert_eq!(
            reference_sequence.index().get(&bin),
            Some(&bgzf::VirtualPosition::from(first_record_in_window))
        );
    }

    #[rstest]
    #[case(Index::Tbi)]
    #[case(Index::Csi)]
    fn test_an_index_finds_a_record_through_its_info_end(#[case] kind: Index) {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let text = [
            "##fileformat=VCFv4.2",
            r#"##INFO=<ID=END,Number=1,Type=Integer,Description="End">"#,
            "##contig=<ID=c,length=100000>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "c\t10\t.\tA\t<DEL>\t.\t.\tEND=40000",
            "c\t30000\t.\tA\tC\t.\t.\t.",
            "c\t50000\t.\tA\tC\t.\t.\t.\n",
        ]
        .join("\n");
        assert!(write(&path, &text, kind, &|| false).unwrap());
        let file = File::open(&path).unwrap();
        let mut reader = match kind {
            Index::Tbi => noodles::vcf::io::IndexedReader::new(
                file,
                tabix::fs::read(kind.path(&path)).unwrap(),
            ),
            Index::Csi => {
                noodles::vcf::io::IndexedReader::new(file, csi::fs::read(kind.path(&path)).unwrap())
            }
        };
        let header = reader.read_header().unwrap();
        let region = "c:35000-35000".parse().unwrap();
        let positions: Vec<usize> = reader
            .query(&header, &region)
            .unwrap()
            .records()
            .map(|record| usize::from(record.unwrap().variant_start().unwrap().unwrap()))
            .collect();
        assert_eq!(positions, [10]);
    }
}
