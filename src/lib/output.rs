//! The VCF output: standard output, a plain file, or a bgzipped and indexed file for a `.gz` path.
use std::fs::File;
use std::io::{BufWriter, ErrorKind, Read, Seek, SeekFrom, Stdout, Write};
use std::path::{Path, PathBuf};

use rust_htslib::bcf::index::{Type, build};
use rust_htslib::bgzf;

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
        let kind = match self {
            Index::Tbi => Type::Tbx,
            Index::Csi => Type::Csi(14),
        };
        build(output, Some(path.as_path()), 0, kind)
            .map_err(|error| error.msg)
            .and_then(|()| {
                File::open(&path)
                    .and_then(|file| file.sync_all())
                    .map_err(|error| error.to_string())
            })
            .map_err(|error| format!("Could not write the index {path:?}: {error}!"))
    }
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
    Bgzf(PathBuf, bgzf::Writer, bool),
}

impl Output {
    /// Creates the output at a path, bgzipped when it ends in `.gz`, or standard output.
    pub fn create(path: Option<&Path>) -> Result<Self, String> {
        let describe = |error: &dyn std::fmt::Display| {
            format!("Could not create the output {path:?}: {error}!")
        };
        match path {
            Some(path) if path.as_os_str() != "-" && has_gzip_ext(&path) => {
                let writer = bgzf::Writer::from_path(path).map_err(|error| describe(&error))?;
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
    /// file that could not be completed. A bgzipped regular file is closed by htslib, which does not
    /// report a failure, so it is read back to check that it ends with the BGZF end-of-file block,
    /// and is then indexed. Returns false, after removing the file and its index, when `stop` asks
    /// the run to stop once the index is written.
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
            Output::Bgzf(path, mut writer, regular) => {
                let flushed = writer.flush();
                drop(writer);
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
    use rust_htslib::bgzf::Reader;
    use rust_htslib::tbx::{self, Read as TbxRead};
    use tempfile::TempDir;

    use super::*;

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
        let mut reader = tbx::Reader::from_path(path).unwrap();
        let tid = reader.tid("c").unwrap();
        reader.fetch(tid, start - 1, end).unwrap();
        reader
            .records()
            .map(|record| {
                let record = String::from_utf8(record.unwrap()).unwrap();
                record.split('\t').nth(1).unwrap().to_string()
            })
            .collect()
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
        assert!(bgzf::is_bgzip(&path).unwrap());
        let mut written = String::new();
        Reader::from_path(&path)
            .unwrap()
            .read_to_string(&mut written)
            .unwrap();
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
        Reader::from_path(dir.path().join(name))
            .unwrap()
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
}
