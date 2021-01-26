//! A library for working with VarDict/VarDictJava output.
#![warn(missing_docs)]
#![warn(missing_doc_code_examples)]
pub mod io;
pub mod records;
pub mod vartovcf;

/// The default log level for the `vartovcf` tool.
pub const DEFAULT_LOG_LEVEL: &str = "info";
