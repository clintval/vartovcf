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
- At this time, only tumor-only mode (`var2vcf_valid.pl`) is supported
- VarDictJava must be run with `--fisher`; rows without the strand-bias p-value and odds-ratio columns are refused

### Filters

No FILTER label is applied unless its option is given; with none given, the FILTER column is `.`. Once any label is applied, a call that fails none of them is `PASS`, and a record without an ALT allele stays `.`. Each label's threshold is written into its `##FILTER` description.

| Label | Option | Applied when |
|---|---|---|
| `NEAR_READ_END` | `--filter-near-read-end <MIN_MEAN_DIST>` | `MEAN_DIST_TO_READ_END` is below the threshold |

### Example Usage

Replace to call to `var2vcf_valid.pl` with `vartovcf` in a typical VarDictJava stream like:

```bash
❯ vardict-java \
    -b input.bam \
    -G hg38.fa \
    -N dna00001 \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | vartovcf --reference hg38.fa --sample dna00001 \
  | bcftools sort -Oz > variants.vcf.gz
```

### Benchmarks

```bash
❯ vartovcf --reference hs38DH.fa --sample dna00001 < test.var > /dev/null
[2025-10-21T01:16:49Z INFO  vartovcf] Input stream: STDIN
[2025-10-21T01:16:49Z INFO  vartovcf] Output stream: STDOUT
[2025-10-21T01:16:49Z INFO  proglog] [main] Processed 60181 variant records

❯ hyperfine --warmup 5 'vartovcf -r hs38DH.fa -s dna00001 < test.var > /dev/null'
Benchmark #1: vartovcf -r /references/hs38DH.fa -s dna00001 < test.var > /dev/null
  Time (mean ± σ):     174.7 ms ±   1.5 ms    [User: 165.7 ms, System: 8.0 ms]
  Range (min … max):   173.2 ms … 178.2 ms    16 runs

❯ hyperfine --warmup 5 'var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null'
Benchmark #1: var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null
  Time (mean ± σ):     359.4 ms ±   2.5 ms    [User: 329.2 ms, System: 25.8 ms]
  Range (min … max):   356.1 ms … 363.6 ms    10 runs
```

### Development

The checks CI runs are cargo aliases, so they can be run locally before pushing:

```bash
❯ cargo ci-fmt
❯ cargo ci-lint
❯ cargo ci-test
```
