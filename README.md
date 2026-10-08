# vartovcf

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/vartovcf/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/vartovcf/badges/version.svg)](http://bioconda.github.io/recipes/vartovcf/README.html)
[![Build Status](https://github.com/clintval/vartovcf/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/clintval/vartovcf/actions/workflows/rust.yml)
[![Coverage Status](https://coveralls.io/repos/github/clintval/vartovcf/badge.svg?branch=main)](https://coveralls.io/github/clintval/vartovcf?branch=main)
[![Language](https://img.shields.io/badge/language-rust-a72144.svg)](https://www.rust-lang.org/)

Convert variants from VarDict/VarDictJava into VCF v4.2 format.

![The Pacific Northwest - Fish Lake](.github/img/cover.jpg)

Install with the Conda or Mamba package manager after setting your [Bioconda channels](https://bioconda.github.io/#usage):

```bash
❯ mamba install vartovcf
```

Or build from source with Rust 1.88 or newer, a C toolchain, and libclang (used by `bindgen` to generate the htslib bindings):

```bash
❯ cargo install --locked --git https://github.com/clintval/vartovcf
```

### Features

- Unlike the Perl script bundled with VarDict, this tool streams record-by-record
- This tool is kept lean on purpose: it converts formats, and applies FILTER labels only when asked
- The output is compliant with the VCF v4.2 and v4.3 specifications
- Output VCF records are unsorted and a call to `bcftools sort` is recommended
- Tumor-only (`var2vcf_valid.pl`) and tumor-normal (`var2vcf_paired.pl`) output are both supported, told apart from the first row
- VarDictJava must be run with `--fisher`

### Example Usage

Replace the call to `var2vcf_valid.pl` with `vartovcf` in a typical VarDictJava stream like the one below. The `--filter-*` values reproduce the thresholds `var2vcf_valid.pl` applies by default; leave out any you don't want, and the label is not applied.

```bash
❯ vardict-java \
    -b input.bam \
    -G hg38.fa \
    -N dna00001 \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | vartovcf --reference hg38.fa --sample dna00001 \
      --filter-near-read-end 8 \
      --filter-low-mean-mapq 10 \
      --filter-homopolymer-indel 13 --filter-homopolymer-indel-max-af 0.275 \
      --filter-tandem-repeat-indel 13 --filter-tandem-repeat-indel-max-af 0.2 \
      --filter-high-mean-mismatches 5.25 \
      --filter-same-read-position 0.35 \
      --filter-strand-bias 0.01 --filter-strand-bias-min-odds-ratio 5 --filter-strand-bias-max-af 0.25 \
  | bcftools sort -Oz > variants.vcf.gz
```

In pileup mode (`-p`) VarDict keeps every candidate and switches off its own call rules. Add `LOW_AF`, `LOW_QMEAN` and `LOW_DP` to mark the candidates its `-f` and `-q` rules would have dropped, here with VarDict's defaults (`-f` as given, `-q 22.5`):

```bash
❯ vardict-java -p ... --fisher calling-intervals.bed \
  | vartovcf --reference hg38.fa --sample dna00001 --skip-non-variants \
      --filter-near-read-end 8 \
      --filter-low-mean-mapq 10 \
      --filter-homopolymer-indel 13 --filter-homopolymer-indel-max-af 0.275 \
      --filter-tandem-repeat-indel 13 --filter-tandem-repeat-indel-max-af 0.2 \
      --filter-high-mean-mismatches 5.25 \
      --filter-same-read-position 0.35 \
      --filter-strand-bias 0.01 --filter-strand-bias-min-odds-ratio 5 --filter-strand-bias-max-af 0.25 \
      --filter-low-af 0.05 \
      --filter-low-qmean 22.5 \
      --filter-low-dp 3 \
  | bcftools sort -Oz > candidates.vcf.gz
```

### Tumor-normal

Replace `var2vcf_paired.pl` the same way. VarDictJava writes the sample names into its output: both as `tumor|normal` when given `-N "tumor|normal"` with a BED file of regions, only the tumor's with `-R`, and without `-N` a name derived from the BAM paths, such as `T.bam|N` from relative paths or just `T` from absolute ones. `vartovcf` reads them from there, but names given with `--sample` and `--normal-sample`, or `--sample "tumor|normal"`, take their place with a warning when they differ; give them so empty input still gets a header for both samples. The `--filter-*` options are the same and test the tumor sample.

```bash
❯ vardict-java \
    -b "tumor.bam|normal.bam" \
    -G hg38.fa \
    -N "dna00001|dna00002" \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | vartovcf --reference hg38.fa --sample dna00001 --normal-sample dna00002 \
      --filter-near-read-end 8 \
      --filter-low-mean-mapq 10 \
      --filter-homopolymer-indel 13 --filter-homopolymer-indel-max-af 0.275 \
      --filter-tandem-repeat-indel 13 --filter-tandem-repeat-indel-max-af 0.2 \
      --filter-high-mean-mismatches 5.25 \
      --filter-same-read-position 0.35 \
      --filter-strand-bias 0.01 --filter-strand-bias-min-odds-ratio 5 --filter-strand-bias-max-af 0.25 \
  | bcftools sort -Oz > somatic.vcf.gz
```

Each record has a tumor and a normal sample column with the tumor-only FORMAT fields. VarDict does not write `HICNT` for a pair, so `vartovcf` recovers each sample's count exactly from its `AD` and the ratio of high- to low-quality ALT reads that VarDict does write, and leaves it `.` in the rare cases that ratio cannot settle, which the header describes. Because `AD` counts ALT reads at any base quality, `HICNT` is the ALT count to compare with a depth filtered by base quality. A sample with no reads of the allele has `.` for its ALT-read statistics and a `0/0` genotype, or `./.` when it has no depth either. Two INFO fields describe the pair:

- `VARDICT_STATUS`: VarDict's label for the allele, such as `StrongSomatic` or `Germline`, from where the allele is found, its AF and VarDict's own call rules rather than a statistical test; its header description gives each label's rule.
- `TUMOR_NORMAL_FISHER_P`: the one-sided Fisher exact p-value that the allele makes up more of the tumor's reads than of the normal's.

FILTER describes the tumor's evidence, not whether a call is somatic, so select somatic calls with these, for example `bcftools view -f PASS -i 'INFO/VARDICT_STATUS ~ "Somatic" && INFO/TUMOR_NORMAL_FISHER_P < 0.05'`.

With `-p`, VarDictJava writes no reference-only rows for a pair.

> [!WARNING]
> **Known VarDictJava issue:** in tumor-normal mode VarDict walks the tumor's alleles at a position from the best-supported down and stops at the first one that fails its call rules, so it never writes the alleles ranked below that one, even ones that pass them and that tumor-only mode would write; and it looks at alleles only the normal carries only when the top-ranked tumor allele fails. VarDictJava 1.6.0 fixed the same stop-at-the-first-failure bug in its other branches but not in this loop (`SomaticPostProcessModule.java:157-160` and `:221` in 1.8.4, and `vardict.pl` before it), so these calls never reach `vartovcf`. On 63 tumor-normal pairs from a targeted panel it dropped 31 of 7,955 rows (0.4%): 3 tumor alleles, all homopolymer stutter, and 28 low-level alleles only the normal carried. A fix is proposed in [joachimwolff/VarDictJava#4](https://github.com/joachimwolff/VarDictJava/pull/4).

### Filters

No FILTER label is applied unless its option is given; with none given, the FILTER column is `.`. Once any label is applied, a call that fails none of them is `PASS`, and a record stays `.` when it has no ALT allele or its tumor (or only) sample has no reads of it. A label is not applied when the value it tests is missing, and for tumor-normal input every label tests the tumor sample. Each label's threshold is written into its `##FILTER` description.

| Label | Option | Applied when |
|---|---|---|
| `NEAR_READ_END` | `--filter-near-read-end <MIN_MEAN_DIST>` | `MEAN_DIST_TO_READ_END` is below the threshold |
| `LOW_MEAN_MAPQ` | `--filter-low-mean-mapq <MIN_MEAN_MAPQ>` | `MEAN_MAPQ` is below the threshold |
| `HOMOPOLYMER_INDEL` | `--filter-homopolymer-indel <MIN_COPIES>` and optionally `--filter-homopolymer-indel-max-af <MAX_AF>` | a one-base insertion or deletion in a homopolymer of at least `MIN_COPIES` copies, at an `AF` below `MAX_AF` if given |
| `TANDEM_REPEAT_INDEL` | `--filter-tandem-repeat-indel <MIN_COPIES>` and optionally `--filter-tandem-repeat-indel-max-af <MAX_AF>` | a one-unit insertion or deletion in a 2-6 bp tandem repeat of at least `MIN_COPIES` copies, at an `AF` below `MAX_AF` if given |
| `HIGH_MEAN_MISMATCHES` | `--filter-high-mean-mismatches <MAX_MEAN_MISMATCHES>` | `MEAN_MISMATCHES` is above the threshold |
| `SAME_READ_POSITION` | `--filter-same-read-position <MAX_AF>` | 2 or more ALT reads all place the variant at the same read position (`ALT_READ_POS_VARIES` 0) and `AF` is below the threshold |
| `STRAND_BIAS` | `--filter-strand-bias <MAX_P>`, optionally `--filter-strand-bias-min-odds-ratio <MIN_ODDS_RATIO>` and `--filter-strand-bias-max-af <MAX_AF>` | `STRAND_BIAS_FISHER_P` is below `MAX_P`, the folded odds ratio of the REF and ALT strand counts is above `MIN_ODDS_RATIO` or the table has an empty cell (if given), and `AF` is below `MAX_AF` (if given) |
| `LOW_AF` | `--filter-low-af <MIN_AF>` | `AF` is below the threshold, useful when VarDict runs with `-p`, which ignores its own `-f` |
| `LOW_QMEAN` | `--filter-low-qmean <MIN_QMEAN>` | `QMEAN` is below the threshold, useful when VarDict runs with `-p`, which keeps calls failing its own `-q` |
| `LOW_DP` | `--filter-low-dp <MIN_DP>` | `DP` is below the threshold |

### Benchmarks

Measured on an Apple M3 Max running macOS 26.6.2, with `vartovcf` 3.0.0 (pre-release build of 5008e5e), VarDictJava 1.8.4 and the `var2vcf_valid.pl` it ships, run by Perl 5.34.1. The input is 60,501 rows of VarDictJava `-p --fisher` output for part of a targeted panel. `vartovcf` streams its input, while `var2vcf_valid.pl` holds all of it in memory to sort it, so its peak memory is 10 times higher (184 MB against 18 MB).

```bash
❯ vartovcf --reference hg19.fa --sample dna00001 < test.var > /dev/null
[2026-10-05T01:04:45Z INFO  vartovcf] Input stream: STDIN
[2026-10-05T01:04:45Z INFO  vartovcf] Output stream: STDOUT
[2026-10-05T01:04:46Z INFO  proglog] [main] Processed 60498 variant records

❯ hyperfine --warmup 5 'vartovcf -r hg19.fa -s dna00001 < test.var > /dev/null'
Benchmark 1: vartovcf -r hg19.fa -s dna00001 < test.var > /dev/null
  Time (mean ± σ):     153.6 ms ±   3.3 ms    [User: 147.6 ms, System: 4.7 ms]
  Range (min … max):   148.5 ms … 160.8 ms    18 runs

❯ hyperfine --warmup 5 'var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null'
Benchmark 1: var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null
  Time (mean ± σ):     386.8 ms ±   9.3 ms    [User: 361.7 ms, System: 18.1 ms]
  Range (min … max):   376.7 ms … 404.1 ms    10 runs

❯ /usr/bin/time -l vartovcf -r hg19.fa -s dna00001 < test.var > /dev/null
...
            18399232  maximum resident set size

❯ /usr/bin/time -l var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null
...
           184303616  maximum resident set size
```

### Development

The checks CI runs are cargo aliases, so they can be run locally before pushing:

```bash
❯ cargo ci-fmt
❯ cargo ci-lint
❯ cargo ci-test
```
