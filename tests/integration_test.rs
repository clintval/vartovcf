#[cfg(test)]
mod tests {
    use std::fs::read_to_string;

    use anyhow::Result;
    use assert_cmd::cmd::Command;
    use assert_cmd::prelude::*;
    use file_diff::diff;
    use rstest::rstest;
    use tempfile::NamedTempFile;

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_success() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_success_on_g_vcf() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.g.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.g.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_on_skippable_records() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.skippable.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.skippable.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_drops_only_exact_duplicate_calls() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.duplicates.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.duplicates.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_keeps_full_precision_and_clamps_af() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.deep.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.deep.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_labels_calls_with_the_given_filters() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.filters.var")
            .arg("--output").arg(output)
            .arg("--filter-near-read-end").arg("8")
            .arg("--filter-low-mean-mapq").arg("10")
            .arg("--filter-homopolymer-indel").arg("13")
            .arg("--filter-homopolymer-indel-max-af").arg("0.275")
            .arg("--filter-tandem-repeat-indel").arg("13")
            .arg("--filter-tandem-repeat-indel-max-af").arg("0.2")
            .arg("--filter-high-mean-mismatches").arg("5.25")
            .arg("--filter-same-read-position").arg("0.35")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.filters.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_on_complex_calls() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.complex.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.complex.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_success_streaming_io() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .pipe_stdin("tests/calls.var")?
            .assert()
            .stdout(read_to_string("tests/calls.vcf")?);
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_success_streaming_io_with_dash() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("-")
            .arg("--output").arg("-")
            .pipe_stdin("tests/calls.var")?
            .assert()
            .stdout(read_to_string("tests/calls.vcf")?);
        Ok(())
    }

    #[rstest]
    #[case("tests/calls.no-fisher.reference.var")]
    #[case("tests/calls.no-fisher.variant.var")]
    #[rustfmt::skip]
    fn run_end_to_end_refuses_rows_without_fisher_columns(#[case] input: &str) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .pipe_stdin(input)?
            .assert()
            .code(1)
            .stdout("");

        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains("run VarDictJava with --fisher!"), "{stderr}");
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_fails_cleanly_on_malformed_row() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .pipe_stdin("tests/calls.malformed.var")?
            .assert()
            .code(1);

        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains("Could not parse column 8 on line 2"), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_fails_cleanly_on_an_sv_without_a_length() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .pipe_stdin("tests/calls.sv-without-length.var")?
            .assert()
            .code(1);

        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains("chr13:24684729: cannot read the DEL length from the genotype column 'G/G'"), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_with_skip_non_variants() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.non-variants-skipped.var")
            .arg("--output").arg(output)
            .arg("--skip-non-variants")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.non-variants-skipped.vcf"));
        Ok(())
    }
}
