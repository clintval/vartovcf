//! A module for commonly used IO functions.
use std::path::PathBuf;
use std::str::FromStr;

/// A path to this system's standard input file descriptor.
pub fn stdin() -> PathBuf { PathBuf::from_str("/dev/stdin").unwrap() }

/// A path to this system's standard output file descriptor.
pub fn stdout() -> PathBuf { PathBuf::from_str("/dev/stdout").unwrap() }
