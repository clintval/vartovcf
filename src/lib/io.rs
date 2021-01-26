//! A module for commonly used IO functions.
use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use crate::path::GZIP_EXTENSION;

/// A path to this system's standard input file descriptor.
pub fn stdin() -> PathBuf { PathBuf::from_str("/dev/stdin").unwrap() }

/// A path to this system's standard output file descriptor.
pub fn stdout() -> PathBuf { PathBuf::from_str("/dev/stdout").unwrap() }

/// Test if a path has a gzip extension or not.
pub fn has_gzip_ext<P>(path: &P) -> bool where P: AsRef<Path> {
    let ext = path.as_ref().extension().and_then(OsStr::to_str);
    matches!(ext, Some(GZIP_EXTENSION))
}
