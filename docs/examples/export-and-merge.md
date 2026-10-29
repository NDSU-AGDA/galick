---
title: Export genotype tables and merge VCF callsets in R
description: Export VCF calls to exchange formats and merge samples or markers with explicit allele reconciliation, phase boundaries and loss reporting.
---
# Export and merge VCFs

Export one VCF to several exchange formats, then combine callsets by
sample or by marker. Everything in this walkthrough is pure R; Beagle,
PLINK, and bcftools are not required for the plain-text examples. Exporting
`.vcf.gz` requires bgzip and produces indexable BGZF output.
Replacing a VCF with `overwrite = TRUE` removes old `.tbi`/`.csi` indexes;
index the new file again before random access.

## Build two small callsets

The two files describe the same three marker records but contain
different samples. At `rs1`, the second file has `ALT=T,G`, so its allele
index `2` refers to `G`; the first file uses `ALT=G`, where `G` is index
`1`. This is the allele-numbering difference that `merge_vcf()` must
reconcile.

```r
library(galick)

work_dir <- tempfile("galick_export_merge_")
dir.create(work_dir)

cohort_a <- file.path(work_dir, "cohort_a.vcf")
cohort_b <- file.path(work_dir, "cohort_b.vcf")

writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=1000>",
  "##contig=<ID=chr2,length=500>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tline_A",
  "chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1",
  "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t1|1",
  "chr2\t50\trs3\tT\tC\t.\tPASS\t.\tGT\t0/."
), cohort_a)

writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=1000>",
  "##contig=<ID=chr2,length=500>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tline_B",
  "chr1\t100\trs1\tA\tT,G\t.\tPASS\t.\tGT\t2/1",
  "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0|1",
  "chr2\t50\trs3\tT\tC\t.\tPASS\t.\tGT\t1/1"
), cohort_b)
```

Validate before transforming anything:

```r
validate_vcf(cohort_a, target = "general")$ok
validate_vcf(cohort_b, target = "general")$ok

detect_subgenome_labels(cohort_a)
detect_genotype_ploidy(cohort_a)
```

Subgenome labels describe observed chromosome names only; a chromosome
subset cannot establish organism ploidy. Genotype ploidy is the maximum
observed allele-slot count in GT and guards diploid-only exports.

## Export to exchange formats

`export_genotypes()` reads the VCF once per call and returns a
`galick_result` report. VCF retains FORMAT fields, phase sets and reference
metadata. CSV and TSV retain GT strings but not additional FORMAT fields.

```r
export_dir <- file.path(work_dir, "exports")
dir.create(export_dir)

export_reports <- list(
  vcf = export_genotypes(
    cohort_a,
    file.path(export_dir, "cohort_a.copy.vcf"),
    format = "vcf"
  ),
  hapmap = export_genotypes(
    cohort_a,
    file.path(export_dir, "cohort_a.hmp.txt"),
    format = "hapmap"
  ),
  csv = export_genotypes(
    cohort_a,
    file.path(export_dir, "cohort_a.csv"),
    format = "csv"
  ),
  tsv = export_genotypes(
    cohort_a,
    file.path(export_dir, "cohort_a.tsv"),
    format = "tsv"
  )
)

vapply(export_reports, function(x) x$summary$sites_out, integer(1))
```

HapMap uses standard IUPAC codes for diploid biallelic SNPs. Phase is lost,
and partially missing calls become `N`; both changes appear in the report.
CSV and TSV retain VCF genotype indices exactly.

```r
read.delim(
  file.path(export_dir, "cohort_a.hmp.txt"),
  check.names = FALSE
)[, c("rs#", "alleles", "line_A")]

read.csv(file.path(export_dir, "cohort_a.csv"))[, c("ID", "line_A")]
```

Numeric output writes samples as rows and markers as columns. The default
`field = "GT"` counts all ALT copies and is invariant to unphased allele order.
Use `field = "DS"` for expected dosage from an imputation program; missing DS
never falls back to GT. Both modes require biallelic markers.

```r
dosage_report <- export_genotypes(
  cohort_a,
  file.path(export_dir, "cohort_a.dosage.tsv"),
  format = "dosage"
)

read.delim(file.path(export_dir, "cohort_a.dosage.tsv"))
```

Pure-R PLINK export writes PED/MAP files. It does not create binary
BED/BIM/FAM files; use `vcf_to_plink()` when that external PLINK
workflow is required.

```r
plink_report <- export_genotypes(
  cohort_a,
  file.path(export_dir, "cohort_a"),
  format = "plink"
)

basename(plink_report$summary$output_paths)
plink_report$warnings$code
```

PED/MAP cannot retain the `|` phase separator. Partially missing genotypes
are written as `0 0`, since PLINK cannot load a half-missing PED call.

## Polyploid refusal is the default

PLINK and HapMap refuse genotypes with more than two allele slots. Numeric
GT counts retain all copies from a K-ploid genotype, with values from 0 to K.

```r
tetraploid <- file.path(work_dir, "tetraploid.vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=1000>",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tline_P",
  "chr1\t300\tpoly1\tA\tG\t.\tPASS\t.\tGT\t0/0/1/1"
), tetraploid)

detect_genotype_ploidy(tetraploid)

tryCatch(
  export_genotypes(
    tetraploid,
    file.path(export_dir, "polyploid"),
    format = "plink"
  ),
  error = function(e) class(e)[1]
)
```

Export the complete allele count without reducing genotype ploidy:

```r
poly_report <- export_genotypes(
  tetraploid,
  file.path(export_dir, "polyploid.tsv"),
  format = "dosage",
  field = "GT"
)

read.delim(file.path(export_dir, "polyploid.tsv"))
```

## Merge different samples

Use `mode = "by_sample"` when every input represents the same marker
records and each file contributes sample columns.

```r
sample_merged <- file.path(work_dir, "merged-by-sample.vcf")
sample_report <- merge_vcf(
  c(cohort_a, cohort_b),
  sample_merged,
  mode = "by_sample",
  sample_collision = "error"
)

unname(sample_report$summary$sites_in)
sample_report$summary$sites_out
sample_report$summary$samples_out
nrow(sample_report$summary$sites_dropped)
```

Inspecting the merged records shows the allele union and reindexing at
`rs1`. The output ALT list is `G,T`; the second file's `2/1` therefore
becomes `1/2`.

```r
merged_lines <- readLines(sample_merged)
cat(merged_lines[!startsWith(merged_lines, "##")], sep = "\n")

validate_vcf(sample_merged, target = "general")$ok
```

The merged VCF feeds directly into existing extraction functions:

```r
vcf_to_dosage(sample_merged)
```

Sample collisions default to an error. Use `sample_collision =
"suffix"` to keep both columns under unique names, or `"coalesce"` to
combine identical calls and fill fully missing calls. Coalescing refuses
two different called genotypes for the same sample, treating unphased allele
permutations as equal. Phased colliding samples require independent phase
reconciliation and cannot be coalesced by this function.

Sample merges recompute `AC` and `AN`, union failed filters, and remove QUAL
and other INFO values with a collected warning. Estimated `AF` is not replaced
with hard-call frequency. All merges retain GT and PS; omitted FORMAT fields
such as DS and GP are listed in `dropped_format_fields`. Separate input phase
blocks receive distinct PS identifiers. Use VCF export to preserve all FORMAT
fields without combining cohorts.

## Merge different markers

Use `mode = "by_marker"` when inputs contain the same samples and
different marker records. Sample columns may be in different orders;
the merge aligns them by name.

```r
marker_a <- file.path(work_dir, "markers_a.vcf")
marker_b <- file.path(work_dir, "markers_b.vcf")

writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=1000>",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tline_A\tline_B",
  "chr1\t300\trs4\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1"
), marker_a)

writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr2,length=500>",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tline_B\tline_A",
  "chr2\t25\trs5\tG\tC\t.\tPASS\t.\tGT\t1/1\t0/1"
), marker_b)

marker_merged <- file.path(work_dir, "merged-by-marker.vcf")
marker_report <- merge_vcf(
  c(marker_a, marker_b),
  marker_merged,
  mode = "by_marker"
)

unname(marker_report$summary$sites_in)
marker_report$summary$sites_out
marker_report$summary$samples_out

vcf_to_hardcalls(marker_merged, encoding = "012")
```

`by_marker` refuses duplicate canonical records across files; distinct SNPs
and indels sharing a coordinate are retained. `by_sample` requires equal
position multiplicities. Repeated coordinates align by compatible allele
identity or unique matching IDs; ambiguous mappings fail. Both modes refuse
conflicting reference dictionaries. Normalize variant representations before
merging when inputs were produced by different callers.

## Read the reports

Both functions return `galick_result` objects. The most useful audit
fields are:

| Field | Meaning |
|---|---|
| `sites_in`, `sites_out` | Record counts before and after the operation. |
| `sites_dropped` | Dropped sites and reasons such as REF mismatch. |
| `samples_in`, `samples_out` | Sample counts before and after the operation. |
| `collisions` | Sample collisions and their resolution. |
| `genotype_ploidy` | Maximum allele-slot count observed in GT. |
| `subgenome_labels` | Observed chromosome suffixes, without a ploidy inference. |
| `dropped_info_fields`, `dropped_format_fields` | Annotations deliberately omitted by a merge. |
| `output_paths` | Files actually written. |

Any site dropped because of a REF mismatch or unresolvable allele is
listed with its input path and original record index in `sites_dropped`.
If filtering removes all records, the result contains errors and no output
is written. Unknown reference-contig lengths remain unknown.

## Format and analysis boundaries

These operations keep variant representations explicit; they do not verify
REF against a genome FASTA, left-align indels, or establish cross-file phase.
Normalize independently produced callsets before merging. All inputs are held
in memory, so test resource use on representative data before a large run.

The format contracts follow the [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf),
[TASSEL HapMap documentation](https://maize-genetics.github.io/tassel/user_manual/load/#hapmap),
and [PLINK PED format](https://www.cog-genomics.org/plink/1.9/formats#ped).
The normalization boundary is motivated by
[Tan et al. (2015), Unified representation of genetic variants](https://doi.org/10.1093/bioinformatics/btv112).
