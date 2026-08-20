# Export and merge VCFs

Export one VCF to several exchange formats, then combine callsets by
sample or by marker. Everything in this walkthrough is pure R; Beagle,
PLINK, and bcftools are not required.

## Build two small callsets

The two files describe the same three marker records but contain
different samples. At `rs1`, the second file has `ALT=T,G`, so its allele
index `2` refers to `G`; the first file uses `ALT=G`, where `G` is index
`1`. This is the allele-numbering difference that `merge_vcf()` must
reconcile.

```r
library(beagler)

work_dir <- tempfile("beagler_export_merge_")
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
#> [1] TRUE
validate_vcf(cohort_b, target = "general")$ok
#> [1] TRUE

detect_chrom_ploidy(cohort_a)
#> [1] 2
detect_genotype_ploidy(cohort_a)
#> [1] 2
```

The two ploidy functions answer different questions. Chromosome ploidy
is inferred from chromosome names; genotype ploidy counts allele slots
in GT and is the value used to guard lossy exports.

## Export to exchange formats

`export_genotypes()` reads the VCF once per call and returns a
`beagler_result` report. VCF, HapMap, CSV, and TSV retain the GT allele
slots and phase separators.

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
#>    vcf hapmap    csv    tsv
#>      3      3      3      3
```

The HapMap export uses explicit allele strings such as `C|G`, instead
of collapsing phased calls to one IUPAC character. CSV and TSV retain
the VCF genotype indices exactly.

```r
read.delim(
  file.path(export_dir, "cohort_a.hmp.txt"),
  check.names = FALSE
)[, c("rs#", "alleles", "line_A")]
#>   rs# alleles line_A
#> 1 rs1     A/G     A/G
#> 2 rs2     C/G     G|G
#> 3 rs3     T/C     T/.

read.csv(file.path(export_dir, "cohort_a.csv"))[, c("ID", "line_A")]
#>    ID line_A
#> 1 rs1    0/1
#> 2 rs2    1|1
#> 3 rs3    0/.
```

Dosage output writes samples as rows and markers as columns:

```r
dosage_report <- export_genotypes(
  cohort_a,
  file.path(export_dir, "cohort_a.dosage.tsv"),
  format = "dosage"
)

read.delim(file.path(export_dir, "cohort_a.dosage.tsv"))
#>   sample rs1 rs2 rs3
#> 1 line_A   1   2  NA
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
#> [1] "cohort_a.ped" "cohort_a.map"
plink_report$warnings$code
#> [1] "export_plink_phase_dropped"
```

The warning is explicit because PED/MAP cannot retain the `|` phase
separator.

## Polyploid refusal is the default

PLINK and the package's dosage table are treated as two-slot outputs.
They refuse a genotype with more than two allele slots unless a
downcoding rule is named explicitly.

```r
tetraploid <- file.path(work_dir, "tetraploid.vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1,length=1000>",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tline_P",
  "chr1\t300\tpoly1\tA\tG\t.\tPASS\t.\tGT\t0/0/1/1"
), tetraploid)

detect_genotype_ploidy(tetraploid)
#> [1] 4

tryCatch(
  export_genotypes(
    tetraploid,
    file.path(export_dir, "polyploid"),
    format = "plink"
  ),
  error = function(e) class(e)[1]
)
#> [1] "beagler_error_export_polyploid_refused"
```

When losing the later allele slots is intentional, opt in:

```r
poly_report <- export_genotypes(
  tetraploid,
  file.path(export_dir, "polyploid"),
  format = "plink",
  downcode = "first_two"
)

poly_report$warnings$code
#> [1] "export_downcoded_genotypes"
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
#> [1] 3 3
sample_report$summary$sites_out
#> [1] 3
sample_report$summary$samples_out
#> [1] 2
nrow(sample_report$summary$sites_dropped)
#> [1] 0
```

Inspecting the merged records shows the allele union and reindexing at
`rs1`. The output ALT list is `G,T`; the second file's `2/1` therefore
becomes `1/2`.

```r
merged_lines <- readLines(sample_merged)
cat(merged_lines[!startsWith(merged_lines, "##")], sep = "\n")
#> #CHROM POS ID  REF ALT QUAL FILTER INFO FORMAT line_A line_B
#> chr1   100 rs1 A   G,T .    PASS   .    GT     0/1    1/2
#> chr1   200 rs2 C   G   .    PASS   .    GT     1|1    0|1
#> chr2   50  rs3 T   C   .    PASS   .    GT     0/.    1/1

validate_vcf(sample_merged, target = "general")$ok
#> [1] TRUE
```

The merged VCF feeds directly into existing extraction functions:

```r
vcf_to_dosage(sample_merged)
#>        rs1 rs2 rs3
#> line_A   1   2  NA
#> line_B   2   1   2
```

Sample collisions default to an error. Use `sample_collision =
"suffix"` to keep both columns under unique names, or `"coalesce"` to
combine identical calls and fill fully missing calls. Coalescing refuses
two different called genotypes for the same sample.

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
#> [1] 1 1
marker_report$summary$sites_out
#> [1] 2
marker_report$summary$samples_out
#> [1] 2

vcf_to_hardcalls(marker_merged, encoding = "012")
#>        rs4 rs5
#> line_A "0" "1"
#> line_B "1" "2"
```

`by_marker` refuses overlapping canonical marker records. `by_sample`
refuses different marker sets. These mode contracts prevent a partial
join from silently dropping or inventing data.

## Read the reports

Both functions return `beagler_result` objects. The most useful audit
fields are:

| Field | Meaning |
|---|---|
| `sites_in`, `sites_out` | Record counts before and after the operation. |
| `sites_dropped` | Dropped sites and reasons such as REF mismatch. |
| `samples_in`, `samples_out` | Sample counts before and after the operation. |
| `collisions` | Sample collisions and their resolution. |
| `genotype_ploidy` | Maximum allele-slot count observed in GT. |
| `chrom_ploidy` | Ploidy inferred independently from chromosome names. |
| `output_paths` | Files actually written. |

Any site dropped because of a REF mismatch or unresolvable allele is
listed in `sites_dropped`; it is never omitted silently.
