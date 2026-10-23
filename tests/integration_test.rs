#[cfg(test)]
mod tests {
    use std::fs::read_to_string;

    use anyhow::Result;
    use assert_cmd::cmd::Command;
    use assert_cmd::prelude::*;
    use file_diff::diff;
    use rstest::rstest;
    use tempfile::{NamedTempFile, TempDir};

    use noodles::sam::alignment::record::data::field::Tag;
    use std::path::{Path, PathBuf};
    use streampile::testing::{Frag, Pair, SamBuilder, Strand};

    #[cfg(unix)]
    use std::os::unix::fs::FileTypeExt;

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
            .arg("--filter-strand-bias").arg("0.01")
            .arg("--filter-strand-bias-min-odds-ratio").arg("5")
            .arg("--filter-strand-bias-max-af").arg("0.25")
            .arg("--filter-low-af").arg("0.0001")
            .arg("--filter-low-qmean").arg("30")
            .arg("--filter-low-dp").arg("3")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.filters.vcf"));
        Ok(())
    }

    #[rstest]
    #[case("--filter-low-hicnt")]
    #[case("--filter-low-hicnt-fraction")]
    fn run_end_to_end_refuses_the_removed_hicnt_filters(
        #[case] option: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args([
                "--reference",
                "tests/reference.fa",
                "--sample",
                "dna00001",
                option,
                "1",
            ])
            .assert()
            .code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(
            stderr.contains(&format!("unexpected argument '{option}'")),
            "{stderr}"
        );
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_takes_the_sample_name_from_the_input() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--input").arg("tests/calls.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.vcf"));
        Ok(())
    }

    #[test]
    fn run_end_to_end_refuses_a_normal_sample_for_tumor_only_input()
    -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args([
                "--reference",
                "tests/reference.fa",
                "--input",
                "tests/calls.var",
            ])
            .args(["--normal-sample", "normal"])
            .assert()
            .code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(
            stderr.contains("--normal-sample was given but the input is tumor-only!"),
            "{stderr}"
        );
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_on_tumor_normal_complex_and_sv_calls() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/tumor-normal.fa")
            .arg("--input").arg("tests/calls.tumor-normal.edge.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.tumor-normal.edge.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_recovers_tumor_normal_hicnt() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/tumor-normal.fa")
            .arg("--input").arg("tests/calls.tumor-normal.hicnt.var")
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.tumor-normal.hicnt.vcf"));
        Ok(())
    }

    #[test]
    fn run_end_to_end_names_the_samples_as_given() -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string("tests/calls.tumor-normal.var")?.replace("T|N\t", "T.bam|N\t");
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args(["--reference", "tests/tumor-normal.fa"])
            .args(["--sample", "T", "--normal-sample", "N"])
            .write_stdin(rows)
            .assert()
            .success();
        let stdout = String::from_utf8(assert.get_output().stdout.clone())?;
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stdout.contains("\tFORMAT\tT\tN\n"), "{stdout}");
        assert!(
            stderr.contains("The input names a sample 'T.bam'; writing it as 'T' as given."),
            "{stderr}"
        );
        Ok(())
    }

    /// The tumor-normal fixture's rows as VarDict writes them when it names only the tumor.
    fn tumor_normal_rows_naming_only_the_tumor() -> String {
        read_to_string("tests/calls.tumor-normal.var")
            .unwrap()
            .replace("T|N\t", "T\t")
    }

    #[rstest]
    #[case(
        None,
        Err("The input is tumor-normal but no --normal-sample was given!")
    )]
    #[case(Some("N"), Ok(()))]
    fn run_end_to_end_needs_the_normal_sample_when_the_input_names_only_the_tumor(
        #[case] normal_sample: Option<&str>,
        #[case] expected: Result<(), &str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args(["--reference", "tests/tumor-normal.fa"])
            .write_stdin(tumor_normal_rows_naming_only_the_tumor());
        if let Some(normal_sample) = normal_sample {
            cmd.arg("--normal-sample").arg(normal_sample);
        }
        match expected {
            Ok(()) => {
                let assert = cmd.assert().success();
                let stdout = String::from_utf8(assert.get_output().stdout.clone())?;
                assert!(stdout.contains("\tFORMAT\tT\tN\n"), "{stdout}");
            }
            Err(message) => {
                let assert = cmd.assert().code(1);
                let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
                assert!(stderr.contains(message), "{stderr}");
            }
        }
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_on_tumor_normal_calls() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/tumor-normal.fa")
            .arg("--input").arg("tests/calls.tumor-normal.var")
            .arg("--output").arg(output)
            .arg("--filter-low-mean-mapq").arg("30")
            .arg("--filter-low-af").arg("0.05")
            .arg("--filter-low-dp").arg("10")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.tumor-normal.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_folds_leaky_mnv_partials() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.leaky-mnv.var")
            .arg("--output").arg(output)
            .arg("--fix-leaky-mnv")
            .arg("--filter-strand-bias").arg("0.01")
            .arg("--filter-low-af").arg("0.02")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.leaky-mnv.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_folds_tumor_normal_leaky_mnv_partials() -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/tumor-normal.fa")
            .arg("--normal-sample").arg("N")
            .arg("--input").arg("tests/calls.tumor-normal.leaky-mnv.var")
            .arg("--output").arg(output)
            .arg("--fix-leaky-mnv")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.tumor-normal.leaky-mnv.vcf"));
        Ok(())
    }

    #[rstest]
    #[case(false, "LOW_AF")]
    #[case(true, "PASS")]
    #[rustfmt::skip]
    fn run_end_to_end_filters_a_folded_mnv_on_its_merged_counts(
        #[case] fix_leaky_mnv: bool,
        #[case] expected: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/tumor-normal.fa")
            .arg("--normal-sample").arg("N")
            .arg("--input").arg("tests/calls.tumor-normal.leaky-mnv.var")
            .arg("--output").arg(output)
            .arg("--filter-low-af").arg("0.51");
        if fix_leaky_mnv {
            cmd.arg("--fix-leaky-mnv");
        }
        cmd.unwrap().assert().success();

        let vcf = std::fs::read_to_string(output)?;
        let filter = vcf
            .lines()
            .find(|line| line.starts_with("cG\t200\t.\tCC\tTT\t"))
            .and_then(|line| line.split('\t').nth(6))
            .expect("the CC>TT record");
        assert_eq!(filter, expected);
        Ok(())
    }

    #[rstest]
    #[case(
        "tests/reference.fa",
        "tests/calls.normalize.var",
        "tests/calls.normalize.vcf"
    )]
    #[case(
        "tests/tumor-normal.fa",
        "tests/calls.tumor-normal.normalize.var",
        "tests/calls.tumor-normal.normalize.vcf"
    )]
    fn run_end_to_end_normalizes_each_record(
        #[case] reference: &str,
        #[case] input: &str,
        #[case] expected: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args([
            "--reference",
            reference,
            "--input",
            input,
            "--output",
            output,
        ])
        .unwrap()
        .assert()
        .success();

        assert!(diff(output, expected));
        Ok(())
    }

    #[rstest]
    #[case(true, "tests/calls.leaky-mnv.overlapping.vcf")]
    #[case(false, "")]
    fn run_end_to_end_writes_a_call_of_overlapping_regions_once(
        #[case] fix_leaky_mnv: bool,
        #[case] expected: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args(["--reference", "tests/reference.fa", "--output", output])
            .args(["--input", "tests/calls.leaky-mnv.overlapping.var"]);
        if fix_leaky_mnv {
            cmd.arg("--fix-leaky-mnv");
        }
        cmd.unwrap().assert().success();

        let vcf = read_to_string(output)?;
        let records: Vec<String> = vcf
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').take(5).collect::<Vec<_>>().join(" "))
            .collect();
        if fix_leaky_mnv {
            assert!(diff(output, expected));
            assert_eq!(records, vec!["chr1 1100 . GG CA"]);
        } else {
            assert_eq!(
                records,
                vec!["chr1 1100 . G C", "chr1 1100 . GG CA", "chr1 1101 . G A"]
            );
        }
        Ok(())
    }

    /// The CHROM, POS, REF, ALT and the first sample's AD of each record converted from some rows.
    fn records(rows: &str, args: &[&str]) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        let reference = if rows.starts_with("T|N") {
            "tests/tumor-normal.fa"
        } else {
            "tests/reference.fa"
        };
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args(["--reference", reference])
            .args(args)
            .write_stdin(rows)
            .assert()
            .success();
        let vcf = String::from_utf8(assert.get_output().stdout.clone())?;
        Ok(vcf
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                let ad = fields[9].split(':').nth(1).unwrap_or_default();
                format!(
                    "{} {} {} {} {ad}",
                    fields[0], fields[1], fields[3], fields[4]
                )
            })
            .collect())
    }

    #[rstest]
    #[case::partial_first(false, true, vec!["chr1 1101 G A .,105"])]
    #[case::mnv_first(true, true, vec!["chr1 1101 G A .,105"])]
    #[case::partial_first_without_the_fold(false, false, vec!["chr1 1101 G A .,100"])]
    #[case::mnv_first_without_the_fold(true, false, vec!["chr1 1101 G A .,100"])]
    fn run_end_to_end_writes_rows_that_normalize_to_one_call_once(
        #[case] reversed: bool,
        #[case] fix_leaky_mnv: bool,
        #[case] expected: Vec<&str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string("tests/calls.leaky-mnv.same-call.var")?;
        let rows: String = if reversed {
            rows.lines().rev().map(|line| format!("{line}\n")).collect()
        } else {
            rows
        };
        let args: &[&str] = if fix_leaky_mnv {
            &["--fix-leaky-mnv"]
        } else {
            &[]
        };
        assert_eq!(records(&rows, args)?, expected);
        Ok(())
    }

    #[rstest]
    #[case::the_most_depth_wins_a_tie_of_reads(
        "tests/calls.same-call.tie.var",
        vec!["chr1 1153 C CA 300,20"]
    )]
    #[case::the_tumor_reads_decide_before_the_normal_reads(
        "tests/calls.tumor-normal.same-call.var",
        vec!["cL 21 CA C 270,50"]
    )]
    fn run_end_to_end_chooses_among_rows_of_one_call_whatever_their_order(
        #[case] input: &str,
        #[case] expected: Vec<&str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string(input)?;
        let reversed: String = rows.lines().rev().map(|line| format!("{line}\n")).collect();
        assert_eq!(records(&rows, &[])?, expected);
        assert_eq!(records(&reversed, &[])?, expected);
        Ok(())
    }

    #[rstest]
    #[case::the_default_rate(&[], vec!["chr1 1100 GG CA .,108"])]
    #[case::a_stricter_rate(
        &["--leaky-mnv-rate", "0.03"],
        vec!["chr1 1100 G C 196,8", "chr1 1100 GG CA .,100"]
    )]
    fn run_end_to_end_folds_only_what_the_leak_rate_explains(
        #[case] rate: &[&str],
        #[case] expected: Vec<&str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string("tests/calls.leaky-mnv.rate.var")?;
        let args: Vec<&str> = ["--fix-leaky-mnv"]
            .into_iter()
            .chain(rate.iter().copied())
            .collect();
        assert_eq!(records(&rows, &args)?, expected);
        Ok(())
    }

    #[rstest]
    #[case::partial_written_by_an_earlier_region(
        "tests/calls.leaky-mnv.partial-first.var",
        vec!["chr1 1100 G C 300,5", "chr1 1100 GG CA .,101"]
    )]
    #[case::partial_taken_by_an_mnv_an_earlier_region_wrote_is_dropped(
        "tests/calls.leaky-mnv.three-regions.var",
        vec!["chr1 1100 GG CA .,100"]
    )]
    fn run_end_to_end_never_counts_a_partial_twice(
        #[case] input: &str,
        #[case] expected: Vec<&str>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string(input)?;
        assert_eq!(records(&rows, &["--fix-leaky-mnv"])?, expected);
        Ok(())
    }

    #[test]
    fn run_end_to_end_keeps_a_reference_only_row_in_mixed_case()
    -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string("tests/calls.mixed-case.var")?;
        assert_eq!(records(&rows, &[])?, vec!["chr1 1100 G . 300"]);
        Ok(())
    }

    #[rstest]
    #[case(
        "tests/tumor-normal.fa",
        "tests/calls.var",
        "Contig 'chr12' is not in the reference: give --reference the FASTA VarDict was run with!"
    )]
    #[case(
        "tests/reference.fa",
        "tests/calls.tumor-normal.var",
        "Contig 'cA' is not in the reference: give --reference the FASTA VarDict was run with!"
    )]
    #[case(
        "tests/missing.fa",
        "tests/calls.var",
        "Could not read the FASTA index \"tests/missing.fa.fai\""
    )]
    fn run_end_to_end_refuses_a_reference_vardict_did_not_use(
        #[case] reference: &str,
        #[case] input: &str,
        #[case] message: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args(["--reference", reference, "--input", input])
            .assert()
            .code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains(message), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        Ok(())
    }

    #[test]
    fn run_end_to_end_refuses_a_ref_allele_the_reference_lacks()
    -> Result<(), Box<dyn std::error::Error>> {
        let rows =
            read_to_string("tests/calls.var")?.replace("\t880\t880\tT\tA\t", "\t880\t880\tG\tA\t");
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
            .write_stdin(rows)
            .assert()
            .code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(
            stderr.contains("The REF allele 'G' at chr1:880 is 'T' in the reference: give --reference the FASTA VarDict was run with!"),
            "{stderr}"
        );
        Ok(())
    }

    #[test]
    fn run_end_to_end_refuses_the_removed_mode_option() -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args(["--reference", "tests/reference.fa", "--mode", "TumorOnly"])
            .assert()
            .code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains("unexpected argument '--mode'"), "{stderr}");
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
    #[case(
        "tests/calls.no-fisher.var",
        &["--reference", "tests/reference.fa", "--sample", "dna00001"],
        "tests/calls.vcf"
    )]
    #[case(
        "tests/calls.tumor-normal.no-fisher.var",
        &[
            "--reference", "tests/tumor-normal.fa",
            "--filter-low-mean-mapq", "30",
            "--filter-low-af", "0.05",
            "--filter-low-dp", "10",
        ],
        "tests/calls.tumor-normal.vcf"
    )]
    #[rustfmt::skip]
    fn run_end_to_end_reads_rows_without_fisher_columns(
        #[case] input: &str,
        #[case] args: &[&str],
        #[case] expected: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .args(args)
            .arg("--input").arg(input)
            .arg("--output").arg(output)
            .unwrap().assert().success();

        assert!(diff(output, expected));
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
        assert!(stderr.contains("chr13:1141: cannot read the DEL length from the genotype column 'G/G'"), "{stderr}");
        assert!(!stderr.contains("panicked"), "{stderr}");
        Ok(())
    }

    #[test]
    fn run_end_to_end_writes_an_ambiguity_code_in_ref_as_its_first_base()
    -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let fasta = dir.path().join("reference.fa");
        let mut bases = std::fs::read("tests/reference.fa")?;
        let offset = 6 + 882 / 60 * 61 + 882 % 60;
        assert_eq!(bases[offset].to_ascii_uppercase(), b'G');
        bases[offset] = b'R';
        std::fs::write(&fasta, bases)?;
        std::fs::copy(
            "tests/reference.fa.fai",
            dir.path().join("reference.fa.fai"),
        )?;
        let rows: String = read_to_string("tests/calls.var")?
            .lines()
            .filter(|row| row.contains("\tchr1\t883\t"))
            .map(|row| {
                let mut fields: Vec<&str> = row.split('\t').collect();
                let gt = format!("R/{}", fields[6]);
                fields[5] = "R";
                fields[13] = &gt;
                format!("{}\n", fields.join("\t"))
            })
            .collect();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .arg("--reference")
            .arg(&fasta)
            .args(["--sample", "dna00001"])
            .write_stdin(rows)
            .assert()
            .success();
        let stdout = String::from_utf8(assert.get_output().stdout.clone())?;
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        let records: Vec<String> = stdout
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').take(5).collect::<Vec<_>>().join(" "))
            .collect();
        assert_eq!(records, vec!["chr1 883 . A T"]);
        assert!(
            stderr.contains("Skipping chr1:883: its REF R is written as A, the same as its ALT A."),
            "{stderr}"
        );
        Ok(())
    }

    #[test]
    fn run_end_to_end_warns_of_an_sv_that_ends_before_its_pos()
    -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string("tests/calls.var")?;
        let mut fields: Vec<&str> = rows.lines().next().unwrap().split('\t').collect();
        fields[4] = "5";
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
            .write_stdin(format!("{}\n", fields.join("\t")))
            .assert()
            .success();
        let stdout = String::from_utf8(assert.get_output().stdout.clone())?;
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stdout.contains("\tTYPE=INV;END=5;"), "{stdout}");
        assert!(
            stderr.contains("The INV at chr12:447 ends at 5, before its POS;"),
            "{stderr}"
        );
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
            .arg("--input").arg("tests/calls.g.var")
            .arg("--output").arg(output)
            .arg("--skip-non-variants")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.non-variants-skipped.vcf"));
        Ok(())
    }

    #[test]
    fn run_end_to_end_sorts_the_records() -> Result<(), Box<dyn std::error::Error>> {
        let rows = read_to_string("tests/calls.var")?;
        let reversed: String = rows.lines().rev().map(|line| format!("{line}\n")).collect();
        let tmp_dir = TempDir::new()?;
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
            .arg("--tmp-dir")
            .arg(tmp_dir.path())
            .write_stdin(reversed)
            .assert()
            .success()
            .stdout(read_to_string("tests/calls.vcf")?);
        assert_eq!(std::fs::read_dir(tmp_dir.path())?.count(), 0);
        Ok(())
    }

    #[rstest]
    #[case::the_option_wins_over_tmpdir(Some("tests"), "missing", true)]
    #[case::tmpdir_when_no_option(None, "tests", true)]
    #[case::an_empty_tmpdir_is_unset(None, "", true)]
    #[case::a_missing_option(Some("missing"), "tests", false)]
    #[case::a_missing_tmpdir(None, "missing", false)]
    #[case::a_file(Some("tests/calls.var"), "", false)]
    fn run_end_to_end_resolves_the_temporary_directory(
        #[case] option: Option<&str>,
        #[case] tmpdir: &str,
        #[case] works: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args([
            "--reference",
            "tests/reference.fa",
            "--input",
            "tests/calls.var",
        ])
        .env("TMPDIR", tmpdir);
        if let Some(option) = option {
            cmd.args(["--tmp-dir", option]);
        }
        if works {
            cmd.assert()
                .success()
                .stdout(read_to_string("tests/calls.vcf")?);
        } else {
            let assert = cmd.assert().code(1).stdout("");
            let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
            let named = option.unwrap_or(tmpdir);
            assert!(
                stderr.contains(&format!(
                    "Could not use the temporary directory \"{named}\""
                )),
                "{stderr}"
            );
        }
        Ok(())
    }

    /// The CHROM to ALT fields of the records a bgzipped VCF's index finds over a 1-based region.
    fn query(
        path: &std::path::Path,
        contig: &str,
        start: u64,
        end: u64,
    ) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        use noodles::core::{Position, Region};
        use noodles::{csi, tabix, vcf};

        let file = std::fs::File::open(path)?;
        let csi = format!("{}.csi", path.display());
        let mut reader = if std::path::Path::new(&csi).exists() {
            vcf::io::IndexedReader::new(file, csi::fs::read(csi)?)
        } else {
            vcf::io::IndexedReader::new(file, tabix::fs::read(format!("{}.tbi", path.display()))?)
        };
        let header = reader.read_header()?;
        let region = Region::new(
            contig,
            Position::try_from(start as usize)?..=Position::try_from(end as usize)?,
        );
        let mut records = Vec::new();
        for record in reader.query(&header, &region)?.records() {
            let record = record?;
            let fields = [
                record.reference_sequence_name().to_string(),
                record.variant_start().transpose()?.unwrap().to_string(),
                Some(record.ids().as_ref())
                    .filter(|ids| !ids.is_empty())
                    .unwrap_or(".")
                    .to_string(),
                record.reference_bases().to_string(),
                record.alternate_bases().as_ref().to_string(),
            ];
            records.push(fields.join(" "));
        }
        Ok(records)
    }

    #[test]
    fn run_end_to_end_indexes_a_bgzipped_output() -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let output = dir.path().join("out.vcf.gz");
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
            .args(["--input", "tests/calls.var"])
            .arg("--output")
            .arg(&output)
            .assert()
            .success();
        assert!(dir.path().join("out.vcf.gz.tbi").exists());
        assert!(!dir.path().join("out.vcf.gz.csi").exists());
        assert_eq!(
            query(&output, "chr1", 881, 883)?,
            ["chr1 883 . G A", "chr1 883 . G T"]
        );
        assert_eq!(
            query(&output, "chr13", 2000, 2000)?,
            ["chr13 1141 . G <DEL>"]
        );
        assert!(query(&output, "chr12", 1, 446)?.is_empty());
        Ok(())
    }

    #[rstest]
    #[case(Some("out.vcf"), &["out.vcf"])]
    #[case(Some("-"), &[])]
    #[case(None, &[])]
    fn run_end_to_end_indexes_no_plain_output(
        #[case] output: Option<&str>,
        #[case] expected: &[&str],
    ) -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let tests = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests");
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.current_dir(dir.path())
            .arg("--reference")
            .arg(tests.join("reference.fa"))
            .arg("--input")
            .arg(tests.join("calls.var"));
        if let Some(output) = output {
            cmd.args(["--output", output]);
        }
        cmd.assert().success();
        let mut files: Vec<String> = std::fs::read_dir(dir.path())?
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        files.sort();
        assert_eq!(files, expected);
        Ok(())
    }

    #[test]
    fn run_end_to_end_writes_a_csi_index_for_a_contig_too_long_for_tabix()
    -> Result<(), Box<dyn std::error::Error>> {
        use std::io::{Seek, SeekFrom, Write};

        let dir = TempDir::new()?;
        let length: u64 = (1 << 29) + 1000;
        let position = length - 500;
        let header = b">chrBig\n";
        let fasta = dir.path().join("reference.fa");
        let mut file = std::fs::File::create(&fasta)?;
        file.write_all(header)?;
        // A sparse file: only the bases around the call are written.
        file.set_len(header.len() as u64 + length + 1)?;
        file.seek(SeekFrom::End(-50_001))?;
        file.write_all(&[b'A'; 50_000])?;
        file.write_all(b"\n")?;
        drop(file);
        std::fs::write(
            dir.path().join("reference.fa.fai"),
            format!(
                "chrBig\t{length}\t{}\t{length}\t{}\n",
                header.len(),
                length + 1
            ),
        )?;

        let row = read_to_string("tests/calls.var")?;
        let mut fields: Vec<String> = row
            .lines()
            .nth(2)
            .unwrap()
            .split('\t')
            .map(String::from)
            .collect();
        fields[2] = "chrBig".to_string();
        fields[3] = position.to_string();
        fields[4] = position.to_string();
        fields[5] = "A".to_string();
        fields[6] = "C".to_string();
        fields[13] = "A/C".to_string();
        fields[34] = format!("chrBig:{}-{}", position - 100, position + 100);

        let output = dir.path().join("out.vcf.gz");
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.arg("--reference")
            .arg(&fasta)
            .args(["--sample", "dna00001"])
            .arg("--output")
            .arg(&output)
            .write_stdin(format!("{}\n", fields.join("\t")))
            .assert()
            .success();
        assert!(dir.path().join("out.vcf.gz.csi").exists());
        assert!(!dir.path().join("out.vcf.gz.tbi").exists());
        assert_eq!(
            query(&output, "chrBig", position, position)?,
            [format!("chrBig {position} . A C")]
        );
        Ok(())
    }

    #[test]
    fn run_end_to_end_removes_the_output_after_an_error() -> Result<(), Box<dyn std::error::Error>>
    {
        let dir = TempDir::new()?;
        for index in ["out.vcf.gz.tbi", "out.vcf.gz.csi"] {
            std::fs::write(dir.path().join(index), "stale")?;
        }
        for name in ["out.vcf", "out.vcf.gz"] {
            let output = dir.path().join(name);
            let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
            cmd.args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
                .args(["--input", "tests/calls.malformed.var"])
                .arg("--output")
                .arg(&output)
                .arg("--tmp-dir")
                .arg(dir.path())
                .assert()
                .code(1);
            assert!(!output.exists());
        }
        assert_eq!(std::fs::read_dir(dir.path())?.count(), 0);
        Ok(())
    }

    #[cfg(unix)]
    #[rstest]
    #[case("INT", 130)]
    #[case("TERM", 143)]
    fn run_end_to_end_cleans_up_when_a_signal_stops_it(
        #[case] signal: &str,
        #[case] code: i32,
        #[values("out.vcf", "out.vcf.gz")] name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Write;
        use std::process::Stdio;
        use std::time::{Duration, Instant};

        let dir = TempDir::new()?;
        let output = dir.path().join(name);
        let mut child =
            std::process::Command::new(assert_cmd::cargo::cargo_bin(env!("CARGO_PKG_NAME")))
                .args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
                .arg("--output")
                .arg(&output)
                .arg("--tmp-dir")
                .arg(dir.path())
                .stdin(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()?;
        let rows = read_to_string("tests/calls.var")?;
        let mut stdin = child.stdin.take().unwrap();
        stdin.write_all(rows.as_bytes())?;
        stdin.flush()?;
        let started = Instant::now();
        while !output.exists() {
            assert!(
                started.elapsed() < Duration::from_secs(30),
                "the output never appeared"
            );
            std::thread::sleep(Duration::from_millis(10));
        }
        let pid = child.id().to_string();
        let status = std::process::Command::new("kill")
            .args([&format!("-{signal}"), &pid])
            .status()?;
        assert!(status.success());
        // The run may already have stopped on a row still in the pipe, which closes it.
        match stdin.write_all(rows.as_bytes()) {
            Err(error) if error.kind() == std::io::ErrorKind::BrokenPipe => {}
            result => result?,
        }
        drop(stdin);
        assert_eq!(child.wait()?.code(), Some(code));
        assert!(!output.exists());
        assert_eq!(std::fs::read_dir(dir.path())?.count(), 0);
        Ok(())
    }

    /// A reference-only row for every base of the test reference's chr1, enough output to fill
    /// several blocks of a bgzipped file.
    fn reference_only_rows() -> String {
        let template = read_to_string("tests/calls.g.var").unwrap();
        let template: Vec<&str> = template.lines().nth(1).unwrap().split('\t').collect();
        let fasta = read_to_string("tests/reference.fa").unwrap();
        let chr1: String = fasta
            .split('>')
            .find(|contig| contig.starts_with("chr1\n"))
            .unwrap()
            .lines()
            .skip(1)
            .collect();
        let mut rows = String::new();
        for (index, base) in chr1.to_uppercase().chars().enumerate() {
            let mut fields = template.clone();
            let position = (index + 1).to_string();
            let base = base.to_string();
            let genotype = format!("{base}/{base}");
            fields[3] = &position;
            fields[4] = &position;
            fields[5] = &base;
            fields[6] = &base;
            fields[13] = &genotype;
            fields[34] = "chr1:1-4000";
            rows.push_str(&fields.join("\t"));
            rows.push('\n');
        }
        rows
    }

    #[cfg(unix)]
    #[rstest]
    #[case("out.vcf")]
    #[case("out.vcf.gz")]
    fn run_end_to_end_removes_an_output_it_cannot_complete(
        #[case] name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Write;
        use std::process::Stdio;

        let dir = TempDir::new()?;
        let rows = reference_only_rows();
        let run = |output: &std::path::Path,
                   limit: Option<u64>|
         -> Result<Option<i32>, Box<dyn std::error::Error>> {
            let script = match limit {
                Some(blocks) => format!("trap '' XFSZ; ulimit -f {blocks}; exec \"$0\" \"$@\""),
                None => "exec \"$0\" \"$@\"".to_string(),
            };
            let mut command = std::process::Command::new("bash");
            command
                .args(["-c", &script])
                .arg(assert_cmd::cargo::cargo_bin(env!("CARGO_PKG_NAME")))
                .args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
                .arg("--output")
                .arg(output)
                .arg("--tmp-dir")
                .arg(dir.path())
                .stdin(Stdio::piped())
                .stderr(Stdio::null());
            if limit.is_some() {
                // The file size limit would also cut short a coverage build's profile of this run.
                command.env("LLVM_PROFILE_FILE", "/dev/null");
            }
            let mut child = command.spawn()?;
            child.stdin.take().unwrap().write_all(rows.as_bytes())?;
            Ok(child.wait()?.code())
        };
        let full = dir.path().join(format!("full.{name}"));
        assert_eq!(run(&full, None)?, Some(0));
        let size = std::fs::metadata(&full)?.len();
        let output = dir.path().join(name);
        assert_eq!(run(&output, Some((size - 1) / 1024))?, Some(1));
        assert!(!output.exists());
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn run_end_to_end_ends_quietly_when_standard_output_closes()
    -> Result<(), Box<dyn std::error::Error>> {
        use std::io::{Read, Write};
        use std::process::Stdio;

        let mut child =
            std::process::Command::new(assert_cmd::cargo::cargo_bin(env!("CARGO_PKG_NAME")))
                .args(["--reference", "tests/reference.fa", "--sample", "dna00001"])
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()?;
        child
            .stdin
            .take()
            .unwrap()
            .write_all(reference_only_rows().as_bytes())?;
        let mut start = [0u8; 16];
        child.stdout.take().unwrap().read_exact(&mut start)?;
        let output = child.wait_with_output()?;
        let stderr = String::from_utf8(output.stderr)?;
        assert_eq!(output.status.code(), Some(0), "{stderr}");
        assert!(!stderr.contains("ERROR"), "{stderr}");
        Ok(())
    }

    #[test]
    fn run_end_to_end_refuses_an_unwritable_output_before_reading_input()
    -> Result<(), Box<dyn std::error::Error>> {
        use std::process::Stdio;
        use std::time::{Duration, Instant};

        let mut child =
            std::process::Command::new(assert_cmd::cargo::cargo_bin(env!("CARGO_PKG_NAME")))
                .args([
                    "--reference",
                    "tests/reference.fa",
                    "--output",
                    "missing/out.vcf",
                ])
                .stdin(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()?;
        let started = Instant::now();
        let status = loop {
            if let Some(status) = child.try_wait()? {
                break status;
            }
            if started.elapsed() > Duration::from_secs(30) {
                child.kill()?;
                panic!("vartovcf waited for input before refusing its output");
            }
            std::thread::sleep(Duration::from_millis(10));
        };
        assert_eq!(status.code(), Some(1));
        Ok(())
    }

    #[cfg(unix)]
    #[rstest]
    #[case("out.vcf")]
    #[case("out.vcf.gz")]
    fn run_end_to_end_writes_to_a_fifo_and_leaves_it_in_place(
        #[case] name: &str,
    ) -> Result<(), Box<dyn std::error::Error>> {
        use std::io::Read;

        let dir = TempDir::new()?;
        let fifo = dir.path().join(name);
        assert!(
            std::process::Command::new("mkfifo")
                .arg(&fifo)
                .status()?
                .success()
        );
        let reader = {
            let fifo = fifo.clone();
            std::thread::spawn(move || std::fs::read(fifo))
        };
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd.args([
            "--reference",
            "tests/reference.fa",
            "--input",
            "tests/calls.var",
        ])
        .arg("--output")
        .arg(&fifo)
        .arg("--tmp-dir")
        .arg(dir.path())
        .assert()
        .success();
        let bytes = reader.join().unwrap()?;
        let text = if name.ends_with(".gz") {
            let copy = dir.path().join("copy.vcf.gz");
            std::fs::write(&copy, &bytes)?;
            let mut text = String::new();
            noodles::bgzf::io::Reader::new(std::fs::File::open(&copy)?)
                .read_to_string(&mut text)?;
            text
        } else {
            String::from_utf8(bytes)?
        };
        assert_eq!(text, read_to_string("tests/calls.vcf")?);
        assert!(std::fs::symlink_metadata(&fifo)?.file_type().is_fifo());
        Ok(())
    }

    #[cfg(unix)]
    #[rstest]
    #[case("/dev/null", false)]
    #[case("/dev/stdout", true)]
    fn run_end_to_end_writes_to_a_device(
        #[case] output: &str,
        #[case] captured: bool,
    ) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        let assert = cmd
            .args([
                "--reference",
                "tests/reference.fa",
                "--input",
                "tests/calls.var",
            ])
            .args(["--output", output])
            .assert()
            .success();
        if captured {
            assert.stdout(read_to_string("tests/calls.vcf")?);
        }
        assert!(std::path::Path::new(output).exists());
        Ok(())
    }

    /// The reference of `tests/reference.fa` over chr1:1081-1120, which holds the fixtures' MNV
    /// GG>CA at 1100-1101.
    const REFERENCE_1081_1120: &str = "TCAAACCTGCCAGGCGGTCGGGTCGCGGACCTCGGTCGAA";

    /// A 40-base read over chr1:1081-1120 with these substitutions at 1-based positions.
    fn read(edits: &[(usize, u8)]) -> String {
        let mut bases = REFERENCE_1081_1120.as_bytes().to_vec();
        for &(position, base) in edits {
            bases[position - 1081] = base;
        }
        String::from_utf8(bases).unwrap()
    }

    /// Base qualities of 40, except at these 1-based positions.
    fn quals(edits: &[(usize, u8)]) -> Vec<u8> {
        let mut quals = vec![40u8; 40];
        for &(position, quality) in edits {
            quals[position - 1081] = quality;
        }
        quals
    }

    /// A fragment over chr1:1081-1120 with these substitutions and qualities, and NM set to the
    /// number of substitutions.
    fn frag(name: &str, edits: &[(usize, u8)], low: &[(usize, u8)]) -> Frag {
        Frag::at(1081)
            .name(name)
            .bases(read(edits))
            .quals(quals(low))
            .attr(Tag::EDIT_DISTANCE, edits.len() as i32)
    }

    /// Writes reads to an indexed BAM in `dir` and returns its path.
    fn write_bam(dir: &Path, name: &str, builder: &SamBuilder) -> PathBuf {
        let path = dir.join(format!("{name}.bam"));
        builder.write_bam(&path).unwrap();
        let index = noodles::bam::fs::index(&path).unwrap();
        noodles::bam::bai::fs::write(dir.join(format!("{name}.bam.bai")), &index).unwrap();
        path
    }

    /// The tumor-only reads of `tests/calls.reads.var`, which is VarDictJava 1.8.4's output for
    /// them with `-q 30 -Q 20 -F 0x504 -z 1 -P 0 -K`: 10 reads of the MNV GG>CA at 1100-1101; G>C
    /// at 1100 on a read whose mate proves the MNV molecule, one with a no-call at 1101 and one with
    /// the reference there; G>A at 1101 on two reads with a no-call at 1100, which VarDict folds
    /// into the MNV itself; G>T at 1100 on two reads; and five reference reads.
    fn tumor_only_reads() -> SamBuilder {
        let mut builder = SamBuilder::new().read_length(40);
        for index in 0..10 {
            builder.add_frag(frag(
                &format!("m{index}"),
                &[(1100, b'C'), (1101, b'A')],
                &[],
            ));
        }
        builder.add_pair(
            Pair::at(1081, 1081)
                .name("p1")
                .bases1(read(&[(1100, b'C'), (1101, b'A')]))
                .quals1(quals(&[(1101, 20)]))
                .bases2(read(&[(1100, b'C'), (1101, b'A')]))
                .quals2(quals(&[]))
                .strand2(Strand::Minus)
                .attr(Tag::EDIT_DISTANCE, 2),
        );
        builder.add_frag(frag("p2", &[(1100, b'C'), (1101, b'N')], &[]));
        builder.add_frag(frag("p3", &[(1100, b'C')], &[]));
        for index in 0..2 {
            builder.add_frag(frag(
                &format!("a{index}"),
                &[(1100, b'N'), (1101, b'A')],
                &[],
            ));
            builder.add_frag(frag(&format!("t{index}"), &[(1100, b'T')], &[]));
        }
        for index in 0..5 {
            builder.add_frag(frag(&format!("r{index}"), &[], &[]));
        }
        builder
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_folds_leaky_mnv_partials_by_their_reads() -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let bam = write_bam(dir.path(), "tumor", &tumor_only_reads());
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--sample").arg("dna00001")
            .arg("--input").arg("tests/calls.reads.var")
            .arg("--output").arg(output)
            .arg("--fix-leaky-mnv")
            .arg("--reads").arg(&bam)
            .arg("--reads-base-quality").arg("30")
            .arg("--filter-low-af").arg("0.1")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.reads.vcf"));
        Ok(())
    }

    /// The tumor and normal reads of `tests/calls.tumor-normal.reads.var`, which is VarDictJava
    /// 1.8.4's output for them with `-q 30 -Q 20 -F 0x504 -z 1 -P 0 -K`, over the three-base MNV
    /// GGG>CAC at 1100-1102: both samples carry it on 10 reads. The tumor has G>C at 1100 on a read
    /// with the reference at 1101 and 1102 and one with a no-call at 1101, and GG>CA at 1100-1101 on
    /// two reads whose mates prove the MNV molecule; the normal has G>C on a read with a no-call at
    /// 1101 and GG>CA on a read with a no-call at 1102.
    fn tumor_normal_reads() -> (SamBuilder, SamBuilder) {
        let mut tumor = SamBuilder::new().read_length(40);
        let mut normal = SamBuilder::new().read_length(40);
        let mnv = [(1100, b'C'), (1101, b'A'), (1102, b'C')];
        for builder in [&mut tumor, &mut normal] {
            for index in 0..10 {
                builder.add_frag(frag(&format!("m{index}"), &mnv, &[]));
            }
            for index in 0..5 {
                builder.add_frag(frag(&format!("r{index}"), &[], &[]));
            }
        }
        tumor.add_frag(frag("s1", &[(1100, b'C')], &[]));
        tumor.add_frag(frag("s2", &[(1100, b'C'), (1101, b'N'), (1102, b'C')], &[]));
        for index in 0..2 {
            tumor.add_pair(
                Pair::at(1081, 1081)
                    .name(format!("p{index}"))
                    .bases1(read(&mnv))
                    .quals1(quals(&[(1102, 20)]))
                    .bases2(read(&mnv))
                    .quals2(quals(&[]))
                    .strand2(Strand::Minus)
                    .attr(Tag::EDIT_DISTANCE, 3),
            );
        }
        normal.add_frag(frag("n1", &[(1100, b'C'), (1101, b'N'), (1102, b'C')], &[]));
        normal.add_frag(frag("n2", &[(1100, b'C'), (1101, b'A'), (1102, b'N')], &[]));
        (tumor, normal)
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_folds_tumor_normal_leaky_mnv_partials_by_their_reads() -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let (tumor, normal) = tumor_normal_reads();
        let tumor = write_bam(dir.path(), "tumor", &tumor);
        let normal = write_bam(dir.path(), "normal", &normal);
        let output = NamedTempFile::new().expect("Cannot create temporary file!");
        let output = output.path().to_str().unwrap();
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--normal-sample").arg("N")
            .arg("--input").arg("tests/calls.tumor-normal.reads.var")
            .arg("--output").arg(output)
            .arg("--fix-leaky-mnv")
            .arg("--reads").arg(format!("{}|{}", tumor.display(), normal.display()))
            .arg("--reads-base-quality").arg("30")
            .unwrap().assert().success();

        assert!(diff(output, "tests/calls.tumor-normal.reads.vcf"));
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_refuses_one_bam_for_two_samples() -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let (tumor, _) = tumor_normal_reads();
        let tumor = write_bam(dir.path(), "tumor", &tumor);
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--normal-sample").arg("N")
            .arg("--input").arg("tests/calls.tumor-normal.reads.var")
            .arg("--output").arg(dir.path().join("out.vcf"))
            .arg("--fix-leaky-mnv")
            .arg("--reads").arg(&tumor);
        let assert = cmd.assert().failure().code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains("--reads names 1 BAM(s) for 2 sample(s)"), "{stderr}");
        Ok(())
    }

    #[test]
    #[rustfmt::skip]
    fn run_end_to_end_refuses_reads_of_another_reference() -> Result<(), Box<dyn std::error::Error>> {
        let dir = TempDir::new()?;
        let mut builder = SamBuilder::new().read_length(40);
        builder.add_frag(frag("r0", &[], &[]).contig(1));
        let bam = write_bam(dir.path(), "tumor", &builder);
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/tumor-normal.fa")
            .arg("--normal-sample").arg("N")
            .arg("--input").arg("tests/calls.tumor-normal.leaky-mnv.var")
            .arg("--output").arg(dir.path().join("out.vcf"))
            .arg("--fix-leaky-mnv")
            .arg("--reads").arg(format!("{0}|{0}", bam.display()));
        let assert = cmd.assert().failure().code(1);
        let stderr = String::from_utf8(assert.get_output().stderr.clone())?;
        assert!(stderr.contains("Contig 'cG' is not in the reads"), "{stderr}");
        Ok(())
    }

    #[rstest]
    #[case::without_the_fix(&["--reads", "x.bam"])]
    #[case::with_a_rate(&["--fix-leaky-mnv", "--reads", "x.bam", "--leaky-mnv-rate", "0.1"])]
    #[case::read_options_without_reads(&["--fix-leaky-mnv", "--reads-base-quality", "30"])]
    #[rustfmt::skip]
    fn run_end_to_end_refuses_read_options_out_of_place(#[case] args: &[&str]) -> Result<(), Box<dyn std::error::Error>> {
        let mut cmd = Command::cargo_bin(env!("CARGO_PKG_NAME"))?;
        cmd
            .arg("--reference").arg("tests/reference.fa")
            .arg("--input").arg("tests/calls.reads.var")
            .args(args);
        cmd.assert().failure().code(1);
        Ok(())
    }
}
