//! A module for commonly used functions in Rust.

/// Namespace for common IO utilities.
pub mod io {
    use std::path::PathBuf;
    use std::str::FromStr;

    /// The default log level for the `vartovcf` tool.
    pub const DEFAULT_LOG_LEVEL: &str = "info";

    /// A path to this system's standard input file descriptor.
    pub fn stdin() -> PathBuf { PathBuf::from_str("/dev/stdin").unwrap() }

    /// A path to this system's standard output file descriptor.
    pub fn stdout() -> PathBuf { PathBuf::from_str("/dev/stdout").unwrap() }
}
