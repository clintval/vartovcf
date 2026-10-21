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

Or build from source with Rust 1.89 or newer, a C toolchain, and libclang (used by `bindgen` to generate the htslib bindings):

```bash
❯ cargo install --locked --git https://github.com/clintval/vartovcf
```

### Features

- Sorts in bounded memory, where VarDict's Perl script holds every record, see [Sorting](#sorting)
- Indexes a bgzipped output file as `tabix -p vcf` does, see [Output](#output)
- Stays lean on purpose: converts formats and applies FILTER labels only when asked
- Writes VCF compliant with the v4.2 and v4.3 specifications
- Normalizes like `bcftools norm -f` and writes uppercase alleles, see [Normalization](#normalization)
- Reads tumor-only and tumor-normal output, detecting which from the first row
- Optionally folds an MNV's partial rows back into it, see [Leaky MNVs](#leaky-mnvs)
- Runs the strand-bias Fisher test itself, so VarDictJava's `-fisher` is optional

### Example Usage

Replace the call to `var2vcf_valid.pl` with `vartovcf` in a typical VarDictJava stream like the one below. The `--filter-*` values reproduce the thresholds `var2vcf_valid.pl` applies by default; leave out any you don't want, and the label is not applied.

```bash
❯ vardict-java \
    -b input.bam \
    -G hg38.fa \
    -N dna00001 \
    -c1 -S2 -E3 -g4 -f0.05 \
    calling-intervals.bed \
  | vartovcf --reference hg38.fa --sample dna00001 \
      --filter-near-read-end 8 \
      --filter-low-mean-mapq 10 \
      --filter-homopolymer-indel 13 --filter-homopolymer-indel-max-af 0.275 \
      --filter-tandem-repeat-indel 13 --filter-tandem-repeat-indel-max-af 0.2 \
      --filter-high-mean-mismatches 5.25 \
      --filter-same-read-position 0.35 \
      --filter-strand-bias 0.01 --filter-strand-bias-min-odds-ratio 5 --filter-strand-bias-max-af 0.25 \
      --output variants.vcf.gz
```

In pileup mode (`-p`) VarDict keeps every candidate and switches off its own call rules. Add `LOW_AF`, `LOW_QMEAN` and `LOW_DP` to mark the candidates its `-f` and `-q` rules would have dropped, here with VarDict's defaults (`-f` as given, `-q 22.5`):

```bash
❯ vardict-java -p ... calling-intervals.bed \
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
      --output candidates.vcf.gz
```

### Tumor-normal

Replace `var2vcf_paired.pl` the same way. VarDictJava writes the sample names into its output: both as `tumor|normal` when given `-N "tumor|normal"` with a BED file of regions, only the tumor's with `-R`, and without `-N` a name derived from the BAM paths, such as `T.bam|N` from relative paths or just `T` from absolute ones. `vartovcf` reads them from there, but names given with `--sample` and `--normal-sample`, or `--sample "tumor|normal"`, take their place with a warning when they differ; give them so empty input still gets a header for both samples. The `--filter-*` options are the same and test the tumor sample.

```bash
❯ vardict-java \
    -b "tumor.bam|normal.bam" \
    -G hg38.fa \
    -N "dna00001|dna00002" \
    -c1 -S2 -E3 -g4 -f0.05 \
    calling-intervals.bed \
  | vartovcf --reference hg38.fa --sample dna00001 --normal-sample dna00002 \
      --filter-near-read-end 8 \
      --filter-low-mean-mapq 10 \
      --filter-homopolymer-indel 13 --filter-homopolymer-indel-max-af 0.275 \
      --filter-tandem-repeat-indel 13 --filter-tandem-repeat-indel-max-af 0.2 \
      --filter-high-mean-mismatches 5.25 \
      --filter-same-read-position 0.35 \
      --filter-strand-bias 0.01 --filter-strand-bias-min-odds-ratio 5 --filter-strand-bias-max-af 0.25 \
      --output somatic.vcf.gz
```

Each record has a tumor and a normal sample column with the tumor-only FORMAT fields. VarDict does not write `HICNT` for a pair, so `vartovcf` recovers each sample's count exactly from its `AD` and the ratio of high- to low-quality ALT reads that VarDict does write, and leaves it `.` in the rare cases that ratio cannot settle, which the header describes. Because `AD` counts ALT reads at any base quality, `HICNT` is the ALT count to compare with a depth filtered by base quality. A sample with no reads of the allele has `.` for its ALT-read statistics and a `0/0` genotype, or `./.` when it has no depth either. Two INFO fields describe the pair:

- `VARDICT_STATUS`: VarDict's label for the allele, such as `StrongSomatic` or `Germline`, from where the allele is found, its AF and VarDict's own call rules rather than a statistical test; its header description gives each label's rule.
- `TUMOR_NORMAL_FISHER_P`: the one-sided Fisher exact p-value that the allele makes up more of the tumor's reads than of the normal's.

FILTER describes the tumor's evidence, not whether a call is somatic, so select somatic calls with these, for example `bcftools view -f PASS -i 'INFO/VARDICT_STATUS ~ "Somatic" && INFO/TUMOR_NORMAL_FISHER_P < 0.05'`.

With `-p`, VarDictJava writes no reference-only rows for a pair.

> [!WARNING]
> **Known VarDictJava issue:** in tumor-normal mode VarDict walks the tumor's alleles at a position from the best-supported down and stops at the first one that fails its call rules, so it never writes the alleles ranked below that one, even ones that pass them and that tumor-only mode would write; and it looks at alleles only the normal carries only when the top-ranked tumor allele fails. VarDictJava 1.6.0 fixed the same stop-at-the-first-failure bug in its other branches but not in this loop (`SomaticPostProcessModule.java:157-160` and `:221` in 1.8.4, and `vardict.pl` before it), so these calls never reach `vartovcf`. In practice it drops about 0.4% of rows, most of them low-level alleles only the normal carries. A fix is proposed in [joachimwolff/VarDictJava#4](https://github.com/joachimwolff/VarDictJava/pull/4).

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

### Leaky MNVs

VarDict builds a multi-nucleotide variant (MNV) from a read's mismatches up to `-X` (default 2) matching bases apart. A read that ends inside the MNV, or has a no-call (N) or a low-quality base there, counts toward a shorter allele instead: an SNV or a shorter MNV carrying the MNV's ALT bases. Here the fourth read's N leaves `TG>AA` beside `TGG>AAA`:

```text
Reference  CGTGGTT
Read 1     CGAAATT
Read 2     CGAAATT
Read 3     CGAAATT
Read 4     CGAANTT
```

VarDict folds some of these partials back itself and writes the rest as calls of their own. `--fix-leaky-mnv` folds them into their MNV with the directional adjacency rule of [UMI-tools](https://doi.org/10.1101/gr.209601.116) (Smith, Heger and Sudbery 2017):

- A row is a partial of an MNV when it lies within the MNV's span and carries the MNV's ALT bases over its own span.
- An MNV claims a partial when it has at least `2 * child - 1` ALT reads in every sample, so an SNV on another haplotype, with about as many reads, stays a record of its own.
- Rows are visited from the most supported down. When two MNVs can claim the same partial, it goes to the one where fewer of the partial's ends fall inside the MNV, since each such end needs a read end or a no-call to explain it; a tie goes to the MNV with more reads.

In the example, `TGG>AAA` is written with `AD` 4 and `DP` 4. A folded record gains INFO `MNV_PARTIALS` (rows folded in) and FORMAT `MNV_PARTIAL_READS` (ALT reads added), and its other fields are combined as VarDict merges an allele's reads:

| Field | On a folded record |
|---|---|
| `AD[1]`, `ADF[1]`, `ADR[1]`, `HICNT` | The sums over the MNV and its partials |
| `DP` | The MNV's, plus the reads of partials starting after `POS` |
| `AF`, `GT`, `FILTER`, `STRAND_BIAS_FISHER_P`, `TUMOR_NORMAL_FISHER_P` | Recomputed from the sums |
| `QMEAN`, `MEAN_MAPQ`, `MEAN_MISMATCHES`, `MEAN_DIST_TO_READ_END` | Means over all the ALT reads, each partial's reads measured at its own bases |
| `ALT_READ_POS_VARIES` | 1 when any row's is, or the rows' `MEAN_DIST_TO_READ_END` differ |
| `REALIGNED_FRAC_OF_DP` | Counts the partials' reads as reassigned to the MNV |
| `REF`, `AD[0]`, `TYPE`, `REPEAT_UNIT_*`, `VARDICT_STATUS` | The MNV's own |

Folding never pushes `AD[1]` above `DP`, though a VarDict row already above it stays so. Only substitutions within one VarDict region are folded. When regions overlap, a partial folded in one is not written again from another, and a partial another region already wrote is not folded, so no read is counted twice. A partial only a later region writes, beside an MNV an earlier region already wrote, is dropped, so its reads are counted nowhere. With few reads a partial cannot be told from an SNV on another haplotype, so a 2-read SNV under a 4-read MNV is folded.

> [!TIP]
> On duplex consensus reads, folded partials carry about 1% of their MNV's reads, while SNVs on another haplotype carry about as many as the MNV and stay apart. Some partials are `PASS` SNVs, often pieces of germline MNVs, that would otherwise read as mutations of their own; folding them typically raises the MNV's `AF` by about 0.006.

> [!WARNING]
> **Known VarDictJava issue:** VarDict's own fold takes a partial starting after the MNV's first base whenever it has fewer reads than the MNV (`VariationRealigner.java:335` in 1.8.4). An SNV on another haplotype with slightly fewer reads is folded in too, which can nearly double the MNV's `AF`, and no row is written for the SNV, so `--fix-leaky-mnv` cannot undo it.

### Normalization

Every record is normalized against `--reference`, which must be the FASTA VarDict was run with: an indel moves to its left-most position and alleles are trimmed to one shared base, as `bcftools norm -f` does, so `bcftools norm` is not needed. A REF allele that differs from the reference is an error. Alleles are written uppercase even where the reference is soft-masked, unlike `bcftools norm` ([samtools/bcftools#2606](https://github.com/samtools/bcftools/issues/2606)). Reference-only rows and symbolic structural variants are not moved.

Each call is written once. Rows of one VarDict region that normalize to one call, such as one indel VarDict placed at two positions in a repeat, hold different reads of it, and the row with the most ALT reads in the tumor (or only) sample is written, then the most in all samples, then the most depth, so the choice never depends on row order. Overlapping BED regions repeat calls, and the first region's record is written. With `--fix-leaky-mnv`, a region's rows are held until it ends, so memory grows with the largest VarDict region.

### Sorting

Records are sorted by contig, in the reference's order, then by POS, REF and ALT, as `bcftools sort` does, so `bcftools sort` is not needed, and they are written only once the input ends. Up to 64 MiB of records are held in memory, plus one VarDict region's rows with `--fix-leaky-mnv`, and the rest spill to `--tmp-dir`, `$TMPDIR` by default, which is cleaned up when `vartovcf` exits or fails. A SIGINT or SIGTERM stops it at its next input row or the next record written and removes the output file, though output already sent to standard output may end in a partial record; a second signal exits at once and leaves its temporary files and an incomplete output file.

### Output

Records go to standard output unless `--output` names a file, which is bgzipped when its name ends in `.gz`.
A bgzipped regular file is indexed once it is complete, as `tabix -p vcf` would index it, so `bcftools index` is not needed.
The index is a tabix index beside it, `<output>.tbi`, or a CSI index, `<output>.csi`, when a reference contig is longer than the 2^29 bases (512 Mbp) a tabix index can address.
Standard output, a plain `.vcf` file, a FIFO and a device such as `/dev/null` get no index.
Any index an earlier run left beside a bgzipped output file is removed when the run starts, so none outlives the file it described.
When the run fails, or a SIGINT or SIGTERM stops it, even while the index is being written, the output file and its index are both removed.

### Benchmarks

Measured on an Apple M3 Max running macOS 26.6.2, with `vartovcf` 3.0.0 (pre-release build of 5008e5e), VarDictJava 1.8.4 and the `var2vcf_valid.pl` it ships, run by Perl 5.34.1. The input is 60,501 rows of VarDictJava `-p -fisher` output for part of a targeted panel. `var2vcf_valid.pl` holds all of its input in memory to sort it, so its peak memory was 10 times that of `vartovcf` 3.0.0, which wrote its records unsorted (184 MB against 18 MB).

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
