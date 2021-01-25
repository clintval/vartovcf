# vartovcf

**Currently under Development**

Convert variants in VarDict/VarDictJava custom format into VCF format.

![The Pacific Northwest - Fish Lake](docs/cover.jpg)

```bash
❯ VarDictJava/build/install/VarDict/bin/VarDict \
    -G /path/to/hg19.fa \
    -N sample1 \
    -b input.bam \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | vartovcf \
  | picard SortVcf -I:/dev/stdin -O:variants.vcf.gz
```

Features:

- Unlike the Perl script bundled with VarDict, this tool streams record-by-record
- Variants are output unsorted and a subsequent call to picard SortVcf is recommended
- Compatibility with the VCF specification v4.2
