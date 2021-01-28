# vartovcf

[![Build Status](https://github.com/clintval/vartovcf/workflows/CI/badge.svg)](https://github.com/clintval/vartovcf/actions)
[![Coverage Status](https://coveralls.io/repos/github/clintval/vartovcf/badge.svg?branch=main)](https://coveralls.io/github/clintval/vartovcf?branch=main)
[![Language](https://img.shields.io/badge/language-rust-a72144.svg)](https://www.rust-lang.org/)

Convert variants from VarDict/VarDictJava into VCF format, fast.

![The Pacific Northwest - Fish Lake](docs/cover.jpg)

```bash
❯ VarDictJava/build/install/VarDict/bin/VarDict \
    -b input.bam \
    -G /references/hg38.fa \
    -N dna00001 \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | vartovcf -r /references/hg38.fa -s dna00001 \
  | picard SortVcf -I:/dev/stdin -O:variants.vcf.gz
```

#### Features

- Unlike the Perl script bundled with VarDict, this tool streams record-by-record
- Only the conversion of variant records generated with tumor-only mode is supported
- Output VCF records are unsorted and a call to `picard SortVcf` is advised
- The output is compliant with the VCF v4.2 and v4.3 specifications

#### Benchmarks

```bash
❯ vartovcf -i test.var -o /dev/null -r /references/hs38DH.fa -s dna00001
[2021-01-28T06:50:44Z INFO  vartovcflib::vartovcf] Input file: "/Users/clintval/test.var"
[2021-01-28T06:50:44Z INFO  vartovcflib::vartovcf] Output file: "/dev/null"
[2021-01-28T06:50:45Z INFO  vartovcflib::vartovcf] Processed 22353 variant records

❯ hyperfine --warmup 5 'vartovcf -i test.var -o /dev/null -r /references/hs38DH.fa' -s dna00001
Benchmark #1: vartovcf -i test.var -o /dev/null -r /references/hs38DH.fa -s dna00001
  Time (mean ± σ):     210.2 ms ±   2.5 ms    [User: 204.4 ms, System: 4.2 ms]
  Range (min … max):   206.9 ms … 218.0 ms    14 runs

❯ hyperfine --warmup 5 'var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null'
Benchmark #1: var2vcf_valid.pl -N dna00001 -f 0.0 -E < test.var > /dev/null
  Time (mean ± σ):     443.2 ms ±  17.0 ms    [User: 416.7 ms, System: 29.1 ms]
  Range (min … max):   419.6 ms … 469.5 ms    10 runs
```
