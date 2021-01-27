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

Features:

- Unlike the Perl script bundled with VarDict, this tool streams record-by-record
- Output VCF records are unsorted and a call to `picard SortVcf` is advised
- The output is compatible with the VCF v4.2 specification

**Currently under Development**
