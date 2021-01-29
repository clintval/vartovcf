//! A module for commonly used IO functions.
use std::ffi::OsStr;
use std::path::Path;

use crate::path::GZIP_EXTENSION;

/// Test if a path has a gzip extension or not.
pub fn has_gzip_ext<P>(path: &P) -> bool
where
    P: AsRef<Path>,
{
    let ext = path.as_ref().extension().and_then(OsStr::to_str);
    matches!(ext, Some(GZIP_EXTENSION))
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_has_gzip_ext() {
        assert!(has_gzip_ext(&PathBuf::from("prefix.gz")));
        assert!(has_gzip_ext(&PathBuf::from("prefix.infix.gz")));
        assert!(!has_gzip_ext(&PathBuf::from("gz")));
        assert!(!has_gzip_ext(&PathBuf::from("prefix.vcf")));
        assert!(!has_gzip_ext(&PathBuf::from("prefix.zip")));
    }
}
