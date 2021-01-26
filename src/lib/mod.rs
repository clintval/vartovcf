//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]
pub mod io;
pub mod fai;
pub mod record;
pub mod vartovcf;

/// The default log level for the `vartovcf` tool.
pub const DEFAULT_LOG_LEVEL: &str = "info";

/// Namespace for path parts and extensions.
mod path {

    /// The Gzip extension.
    pub const GZIP_EXTENSION: &str = "gz";
}
