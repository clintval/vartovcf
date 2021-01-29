#[cfg(test)]
mod tests {
    use anyhow::Result;
    use assert_cmd::prelude::*;
    use std::process::Command;
    use tempfile::NamedTempFile;

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_success() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file.");
        let output = &output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--input").arg("tests/nras.var")
            .arg("--output").arg(output)
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .unwrap().assert().success();

        Ok(())
    }
}
