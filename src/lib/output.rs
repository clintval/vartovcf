//! The VCF output: standard output, a plain file, or a bgzipped file for a `.gz` path.
use std::fs::File;
use std::io::{BufWriter, ErrorKind, Read, Seek, SeekFrom, Stdout, Write};
use std::path::{Path, PathBuf};

use rust_htslib::bgzf;

use crate::io::has_gzip_ext;

/// The empty block every complete BGZF file ends with.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Where the VCF text goes. A file is created at once and, if it is a regular file, removed again
/// unless the run finishes; a device, pipe or FIFO given as the output is only ever closed.
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
    /// report a failure, so it is read back to check that it ends with the BGZF end-of-file block.
    pub fn finish(self) -> Result<(), String> {
        let describe =
            |error: &dyn std::fmt::Display| format!("Could not write the output: {error}!");
        let (path, result) = match self {
            Output::Stdout { mut writer, .. } => {
                return match writer.flush() {
                    Err(error) if error.kind() != ErrorKind::BrokenPipe => Err(describe(&error)),
                    _ => Ok(()),
                };
            }
            Output::Plain(path, writer, regular) => {
                let result = writer
                    .into_inner()
                    .map_err(|error| error.into_error())
                    .and_then(|file| if regular { file.sync_all() } else { Ok(()) });
                (regular.then_some(path), result)
            }
            Output::Bgzf(path, mut writer, regular) => {
                let flushed = writer.flush();
                drop(writer);
                let result = flushed.and_then(|()| {
                    if regular {
                        check_bgzf_end(&path)
                    } else {
                        Ok(())
                    }
                });
                (regular.then_some(path), result)
            }
        };
        result.map_err(|error| {
            if let Some(path) = path {
                let _ = std::fs::remove_file(path);
            }
            describe(&error)
        })
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
    use rust_htslib::bgzf::Reader;
    use tempfile::TempDir;

    use super::*;

    #[test]
    fn test_writes_a_plain_file() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf");
        let mut output = Output::create(Some(&path)).unwrap();
        output.write_all(b"line\n").unwrap();
        output.finish().unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "line\n");
    }

    #[test]
    fn test_writes_a_bgzipped_file_for_a_gz_path() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("out.vcf.gz");
        let mut output = Output::create(Some(&path)).unwrap();
        output.write_all(b"line\n").unwrap();
        output.finish().unwrap();
        assert!(bgzf::is_bgzip(&path).unwrap());
        let mut text = String::new();
        Reader::from_path(&path)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        assert_eq!(text, "line\n");
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
