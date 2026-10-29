# vartovcf

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/vartovcf/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/vartovcf/badges/version.svg)](http://bioconda.github.io/recipes/vartovcf/README.html)
[![Build Status](https://github.com/clintval/vartovcf/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/clintval/vartovcf/actions/workflows/rust.yml)
[![Coverage Status](https://coveralls.io/repos/github/clintval/vartovcf/badge.svg?branch=main)](https://coveralls.io/github/clintval/vartovcf?branch=main)
[![Language](https://img.shields.io/badge/language-rust-a72144.svg)](https://www.rust-lang.org/)

Convert variants from VarDict/VarDictJava into VCF v4.2 format.

![The Pacific Northwest - Fish Lake](.github/img/cover.jpg)

Install with the Conda or Mamba package manager after setting your [Bioconda channels](https://bioconda.github.io/#usage):

```console
mamba install vartovcf
```

## Features

- Output is always coordinate sorted and normalized, and indexed when bgzipped
- Reads tumor-only and tumor-normal output, see [Tumor-Normal](#tumor-normal)
- Applies FILTER labels only when asked, see [Filters](#filters)
- Optionally folds an MNV's partial rows back into it, see [Leaky MNVs](#leaky-mnvs)
- Writes gVCFs with a record and depth at every position, see [gVCFs](#gvcfs)

## Quick Start

Replace the call to `var2vcf_valid.pl` with `vartovcf`:

```console
vardict-java \
    -b input.bam \
    -G hg38.fa \
    -N dna00001 \
    -c1 -S2 -E3 -g4 -f0.05 \
    calling-intervals.bed \
  | vartovcf \
      --reference hg38.fa \
      --sample dna00001 \
      --output variants.vcf.gz
```

The reference must be the one VarDict ran with.
An output ending in `.gz` is bgzipped and indexed; without `--output`, plain VCF goes to standard output.

## Filters

No FILTER label is applied unless its option is given.
These are the thresholds `var2vcf_valid.pl` applies by default, though not every rule around them is reproduced:

```console
vartovcf \
    --reference hg38.fa \
    --sample dna00001 \
    --filter-near-read-end 8 \
    --filter-low-mean-mapq 10 \
    --filter-homopolymer-indel 13 --filter-homopolymer-indel-max-af 0.275 \
    --filter-tandem-repeat-indel 13 --filter-tandem-repeat-indel-max-af 0.2 \
    --filter-high-mean-mismatches 5.25 \
    --filter-same-read-position 0.35 \
    --filter-strand-bias 0.01 --filter-strand-bias-min-odds-ratio 5 --filter-strand-bias-max-af 0.25 \
    --output variants.vcf.gz
```

In pileup mode (`-p`), VarDict keeps every candidate, so add these to label the ones its `-f`, `-q` and `-r`, and `var2vcf_valid.pl`'s `-d`, would drop:

```console
    --filter-low-af 0.05 \
    --filter-low-qmean 22.5 \
    --filter-low-hicnt 2 \
    --filter-low-dp 3
```

| Label | Applied when |
|---|---|
| `NEAR_READ_END` | `MEAN_DIST_TO_READ_END` is below the threshold |
| `LOW_MEAN_MAPQ` | `MEAN_MAPQ` is below the threshold |
| `HOMOPOLYMER_INDEL` | An indel or Complex call one base longer or shorter than REF sits in a homopolymer of at least the threshold's copies |
| `TANDEM_REPEAT_INDEL` | An indel or Complex call one unit longer or shorter than REF sits in a 2-6 bp tandem repeat of at least the threshold's copies |
| `HIGH_MEAN_MISMATCHES` | `MEAN_MISMATCHES` is above the threshold |
| `SAME_READ_POSITION` | 2 or more ALT reads all place the variant at one read position, at an `AF` below the threshold |
| `STRAND_BIAS` | `STRAND_BIAS_FISHER_P` is below the threshold |
| `LOW_AF` | `AF` is below the threshold |
| `LOW_QMEAN` | `QMEAN` is below the threshold |
| `LOW_HICNT` | `HICNT` is below the threshold |
| `LOW_DP` | `DP` is below the threshold |

The `-max-af` and `-min-odds-ratio` options narrow a label to calls below that `AF` or above that odds ratio.
Once any label is applied, a call that fails none is `PASS`.

## Tumor-Normal

Replace `var2vcf_paired.pl` the same way, naming the normal sample:

```console
vardict-java \
    -b "tumor.bam|normal.bam" \
    -G hg38.fa \
    -N "dna00001|dna00002" \
    -c1 -S2 -E3 -g4 -f0.05 \
    calling-intervals.bed \
  | vartovcf \
      --reference hg38.fa \
      --sample dna00001 \
      --normal-sample dna00002 \
      --output somatic.vcf.gz
```

Each record gets INFO `VARDICT_STATUS`, VarDict's label such as `StrongSomatic` or `Germline`, and `TUMOR_NORMAL_FISHER_P`, the one-sided Fisher p-value that the allele is more of the tumor's reads than the normal's.
FILTER labels test the tumor only, so select somatic calls with these:

```console
bcftools view -f PASS -i 'INFO/VARDICT_STATUS ~ "Somatic" && INFO/TUMOR_NORMAL_FISHER_P < 0.05' somatic.vcf.gz
```

> [!WARNING]
> VarDictJava's tumor-normal mode stops at the first allele at a position that fails its call rules, so it never writes the alleles ranked below it, about 0.4% of rows.
> A fix is proposed in [joachimwolff/VarDictJava#4](https://github.com/joachimwolff/VarDictJava/pull/4).

## Leaky MNVs

VarDict builds a multi-nucleotide variant (MNV) from a read's nearby mismatches.
A read that ends inside the MNV, or has a no-call or low-quality base there, counts toward a shorter allele instead.
Here reads 3 and 4 leave `G>C` and `G>A` beside the MNV `GG>CA`:

```text
Reference  GTCGGGT
Read 1     GTCCAGT
Read 2     GTCCAGT
Read 3     GTCCNGT
Read 4     GTCNAGT
```

The examples run from a clone of this repository on its test data:

```console
vartovcf --reference tests/reference.fa --input tests/calls.leaky-mnv.var --output calls.vcf.gz
bcftools query -r chr1:1100-1101 -f '%POS\t%REF\t%ALT[\t%AD]\n' calls.vcf.gz
```

```text
1100	G	C	196,6
1100	G	T	196,30
1100	GG	CA	.,200
1101	G	A	196,4
```

With `--fix-leaky-mnv`, the 6-read and 4-read partials fold back into the MNV, and `G>T`, another allele, stays:

```console
vartovcf --reference tests/reference.fa --input tests/calls.leaky-mnv.var --fix-leaky-mnv --output folded.vcf.gz
bcftools query -r chr1:1100-1101 -f '%POS\t%REF\t%ALT\t%MNV_PARTIALS[\t%AD\t%MNV_PARTIAL_READS]\n' folded.vcf.gz
```

```text
1100	G	T	.	196,30	.
1100	GG	CA	2	.,210	10
```

A partial folds when it carries the MNV's ALT bases and has few enough reads to be leakage at `--leaky-mnv-rate`, 0.05 by default.
At that rate, an MNV folds a 1-read partial from 2 reads up, a 2-read one from 6, and a 10-read one from 101.
A folded record sums the ALT reads, recomputes `AF`, `GT`, FILTER and the Fisher p-values, and counts what it took in INFO `MNV_PARTIALS` and FORMAT `MNV_PARTIAL_READS`.

> [!TIP]
> To estimate your library's rate, run without `--fix-leaky-mnv` and, for MNVs with a few hundred reads, take each partial's `AD[1]` over the sum of its and the MNV's.
> Set the rate at a high quantile of those, such as the 99th.

> [!WARNING]
> VarDictJava folds a partial starting after the MNV's first base whenever it has fewer reads, even an SNV of another haplotype, and writes no row for it, so `--fix-leaky-mnv` cannot undo it.

### Deciding Partials by Their Reads

With `--reads`, a partial is decided by its reads in the BAM VarDict read, instead of by their count:

```console
vardict-java ... \
  | vartovcf \
      --reference hg38.fa \
      --sample dna00001 \
      --fix-leaky-mnv \
      --reads input.bam \
      --output variants.vcf.gz
```

Each template of a partial is read at the MNV's other changed positions:

| There, the template has | It is | Its reads |
|---|---|---|
| The MNV's ALT | The MNV's molecule | Move to the MNV |
| Another base | Another allele | Stay |
| A no-call, a low-quality base, or no read | Unknown | Follow the row's other reads |

A partial with reads left keeps only them, and FORMAT `MNV_MOVED_READS` counts the reads that left.
Give tumor-normal BAMs as VarDict takes them, `tumor.bam|normal.bam`, and set these to the run's VarDict options:

| Option | VarDict's | Default |
|---|---|---|
| `--reads-base-quality` | `-q` | 22.5 |
| `--reads-min-mapq` | `-Q` | 0 |
| `--reads-exclude-flags` | `-F` | `0x504` |
| `--reads-max-mismatches` | `-m` | 8 |
| `--reads-mnv-gap` | `-X` | 2 |

## gVCFs

In pileup mode (`-p`), VarDict leaves some positions of its regions without a record.
`--fill-gvcf` writes a reference-only record at each one, so every position has a record, and `--context-3mer` adds each record's reference 3-mer:

```console
vardict-java -p ... calling-intervals.bed \
  | vartovcf \
      --reference hg38.fa \
      --sample dna00001 \
      --fill-gvcf \
      --fill-gvcf-bed calling-intervals.bed \
      --context-3mer \
      --output dna00001.g.vcf.gz
```

A filled record's INFO `GVCF_FILL` says where its depth came from:

| `GVCF_FILL` | The position |
|---|---|
| `ROWS` | Has reads VarDict counted in rows written elsewhere, such as a deletion it writes at the base before |
| `SPANNED` | Lies inside a deletion, MNV or Complex call, whose reads are its depth |
| `NO_READS` | Has no reads |

On the test data, a 3-base deletion at 1106 leaves positions to fill inside it, and no reads reach past 1111:

```console
vartovcf --reference tests/reference.fa --input tests/calls.fill-gvcf.var --fill-gvcf --context-3mer --output calls.g.vcf.gz
bcftools query -r chr1:1106-1112 -f '%POS\t%REF\t%ALT\t%CONTEXT_3mer\t%GVCF_FILL[\t%AD\t%SITE_DP]\n' calls.g.vcf.gz
```

```text
1106	C	.	GCG	.	280	280
1106	CGGA	C	GCG	.	0,280	280
1107	G	.	CGG	ROWS	0	280
1108	G	.	GGA	SPANNED	0	280
1109	A	.	GAC	SPANNED	0	280
1110	C	.	ACC	.	250	250
1111	C	.	CCT	.	240	240
1112	T	.	CTC	NO_READS	0	0
```

`--fill-gvcf-bed` also fills the regions of VarDict's BED that no read covers.
Set `--fill-gvcf-bed-extend` to VarDict's `-x`, and add `--fill-gvcf-bed-one-based` unless VarDict ran with `-z 1` or read a 4-column BED without `-c`; a mismatch fails the run, naming the options that fix it.
VarDict must run without `-3`, which moves indels from where it counted them.
Tumor-normal input is refused, as VarDict writes no reference-only rows for a pair.

## Site Depth

Every record carries FORMAT `SITE_DP`, the depth at its POS: the reads VarDict counted there as REF or any allele, plus those of deletions, MNVs and Complex calls spanning it.
It is the same on every record at a position, so `uniq` leaves one line per position, and it is complete only for pileup (`-p`) output.
Unlike `DP`, it never counts the `N` bases VarDict's `-K` adds.

Take depths from the whole file, before leaving out any record:

```console
# Depth at each position
bcftools query -f '%CHROM\t%POS[\t%SITE_DP]\n' calls.g.vcf.gz | uniq

# Total depth
bcftools query -f '%CHROM\t%POS[\t%SITE_DP]\n' calls.g.vcf.gz | uniq | awk '{n += $3} END {print n}'

# Depth by reference 3-mer, a mutation spectrum's denominators
bcftools query -f '%CHROM\t%POS\t%CONTEXT_3mer[\t%SITE_DP]\n' calls.g.vcf.gz | uniq \
  | awk '{n[$3] += $4} END {for (c in n) print c, n[c]}'

# Mutations with their context and ALT reads
bcftools query -i 'ALT!="."' -f '%CHROM\t%POS\t%REF\t%ALT\t%CONTEXT_3mer[\t%AD{1}]\n' calls.g.vcf.gz
```

MF is the number of mutations (MF min) or their ALT reads (MF max) over the total depth.
A deletion's `AD` is counted at its first deleted base, not POS, so its `AD[1]` over `SITE_DP` can top 1.

## Development and Testing

The checks CI runs are cargo aliases, so they can be run locally before pushing:

```console
cargo ci-fmt
cargo ci-lint
cargo ci-test
```
