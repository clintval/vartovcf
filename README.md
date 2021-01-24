# var-to-vcf

Convert variants in VarDict/VarDictJava custom format into VCF format

```bash
❯ VarDictJava/build/install/VarDict/bin/VarDict \
    -G /path/to/hg19.fa \
    -N sample1 \
    -b input.bam \
    -c1 -S2 -E3 -g4 -f0.05 \
    --fisher \
    calling-intervals.bed \
  | var-to-vcf -N sample1 -A -E -f0.05 \
  | picard SortVcf -I:/dev/stdin -O:/dev/stdout
  > variants.vcf
```

Features:

- Unlike the Perl script bundled with VarDict, this tool streams
- Variants are output unsorted in the VCF and picard SortVcf is recommended
- Faster and more intuitive in its use than the bundled VarDict Perl script
- Compatibility with the VCF specification v4.2
