# vartovcf

[![Build Status](https://github.com/clintval/vartovcf/workflows/CI/badge.svg)](https://github.com/clintval/vartovcf/actions)
[![Coverage Status](https://coveralls.io/repos/github/clintval/vartovcf/badge.svg?branch=main)](https://coveralls.io/github/clintval/vartovcf?branch=main)
[![Language](https://img.shields.io/badge/language-rust-a72144.svg)](https://www.rust-lang.org/)

Convert variants from VarDict/VarDictJava into VCF v4.2 format.

![The Pacific Northwest - Fish Lake](.github/img/cover.jpg)

```bash
❯ VarDictJava/build/install/VarDict/bin/VarDict \
    -b input.bam \
    -G /references/hg38.fa \
    -N dna00001 \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | vartovcf -r /references/hg38.fa -s dna00001 \
  | bcftools sort -Ou \
  | bcftools norm -Ou --remove-duplicates \
  | bgzip -dc > variants.vcf.gz
```

#### Features

- Unlike the Perl script bundled with VarDict, this tool streams record-by-record
- This tool is kept lean on purpose and is solely responsible for fast format conversion 
- The output is compliant with the VCF v4.2 and v4.3 specifications
- Output VCF records are unsorted and a call to `bcftools sort` is advised
- Output VCF records may also have duplicates so a call to `bcftools norm -D` is advised
- At this time, only tumor-only mode (`var2vcf_valid.pl`) is supported

#### Benchmarks

```bash
❯ vartovcf -r /references/hs38DH.fa -s dna00001 < test.var > /dev/null
[2021-01-31T03:45:12Z INFO  vartovcf] Input stream: STDIN
[2021-01-31T03:45:12Z INFO  vartovcf] Output stream: STDOUT
[2021-01-31T03:45:12Z INFO  vartovcflib::progress_logger] Logging records started at: 2021-01-30 22:45:12.
[2021-01-31T03:45:17Z INFO  vartovcflib::progress_logger] Processed 100000 variant records. Elapsed time: 00:00:05. Time of last 100000 variant records: 00:00:05
[2021-01-31T03:45:21Z INFO  vartovcflib::progress_logger] Processed 156471 variant records. Elapsed time: 00:00:08. Time of last 56471 variant records: 00:00:03

❯ hyperfine --warmup 5 'vartovcf -r /references/hs38DH.fa -s dna00001 < test.var > /dev/null'
Benchmark #1: vartovcf -r /references/hs38DH.fa -s dna00001 < test.var > /dev/null
  Time (mean ± σ):      1.565 s ±  0.043 s    [User: 1.543 s, System: 0.019 s]
  Range (min … max):    1.528 s …  1.679 s    10 runs

❯ hyperfine --warmup 5 'var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null'
Benchmark #1: var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null
  Time (mean ± σ):      1.815 s ±  0.021 s    [User: 1.680 s, System: 0.128 s]
  Range (min … max):    1.789 s …  1.848 s    10 runs
```
