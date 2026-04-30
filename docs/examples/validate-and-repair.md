# Validate and repair

The `validate_vcf` → `repair_vcf` → `validate_vcf` workflow on a
deliberately broken fixture.

## A dirty starting point

`inst/extdata/tiny_dirty.vcf` ships nine records that exercise every
drop condition `repair_vcf()` knows about:

```r
library(beagler)

dirty <- system.file("extdata", "tiny_dirty.vcf", package = "beagler")
cat(readLines(dirty), sep = "\n")
```

```text
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID    REF ALT QUAL FILTER INFO FORMAT s1     s2
chr1   200 rs1   A   T   .    PASS   .    GT     0/0    1/1
chr1   100 .     C   G   .    PASS   .    GT     0/1    1/1
chr1   300 rs3   G   A,C .    PASS   .    GT     0/0    0/1
chr1   100 rs4   T   A   .    PASS   .    GT     1/1    0/0
chr2   -5  rs5   A   T   .    PASS   .    GT     0/0    0/1
chr2   50  rs6   G   G   .    PASS   .    GT     0/0    0/1
chr2   100 rs1   C   A   .    PASS   .    GT     0/0    1/1
chr2   200 rs8   T   G   .    PASS   .    GT     0      1
chr2   300 rs9   A   C   .    PASS   .    GT:DP  0/1:30 1/1:25
```

## Validate first

```r
v <- validate_vcf(dirty, target = "beagle")
v$ok
#> [1] FALSE

v$errors$code
#>  [1] "vcf_invalid_pos"            "vcf_duplicated_marker_ids"
#>  [3] "vcf_duplicated_chrom_pos"   "vcf_multiallelic"
#>  [5] "vcf_beagle_ref_eq_alt"      "vcf_beagle_non_diploid"
#>  [7] "vcf_unsorted_records"
v$warnings$code
#> [1] "vcf_missing_marker_ids" "vcf_beagle_not_bgzf"
```

## Repair

```r
clean <- tempfile(fileext = ".vcf")
result <- repair_vcf(dirty, clean)
result$ok
#> [1] TRUE
```

`$summary$step_summary` records what each step did:

```r
result$summary$step_summary
#> # A tibble: 5 × 4
#>   step                  n_in  n_out n_changed
#>   <chr>                <int> <int>     <int>
#> 1 drop_invalid_records     9     6         3   # bad pos, REF==ALT, non-diploid
#> 2 filter_biallelic         6     5         1   # multiallelic ALT
#> 3 fix_marker_ids           5     5         1   # filled '.' ID
#> 4 remove_duplicates        5     3         2   # dup ID + dup (chrom,pos)
#> 5 sort                     3     3         3   # ordered records
```

Step actions also surface as `$warnings`:

```r
result$warnings$code
#> [1] "repair_vcf_drop_invalid_records_applied"
#> [2] "repair_vcf_filter_biallelic_applied"
#> [3] "repair_vcf_fix_marker_ids_applied"
#> [4] "repair_vcf_remove_duplicates_applied"
#> [5] "repair_vcf_sort_applied"
```

## Re-validate

```r
validate_vcf(clean, target = "beagle")$ok
#> [1] TRUE
```

## What lands in the output

`repair_vcf()` writes a *minimal* VCF: just `CHROM`, `POS`, `ID`,
`REF`, `ALT`, and `FORMAT=GT`. The `GT:DP` cell from the dirty input
gets rewritten to GT-only:

```r
cat(readLines(clean)[6:9], sep = "\n")
#> #CHROM  POS  ID                REF  ALT  QUAL  FILTER  INFO  FORMAT  s1   s2
#> chr1    100  chr1_100_C_G      C    G    .     PASS    .     GT      0/1  1/1
#> chr1    200  rs1               A    T    .     PASS    .     GT      0/0  1/1
#> chr2    300  rs9               A    C    .     PASS    .     GT      0/1  1/1
```

## Step selection

You can opt out of steps. If the remaining steps can't make the
output pass Beagle validation, the function reports
`repair_vcf_remaining_*` errors and refuses to write:

```r
out <- tempfile(fileext = ".vcf")
result <- repair_vcf(dirty, out, steps = "drop_invalid_records")
result$ok
#> [1] FALSE
result$errors$code
#> [1] "repair_vcf_remaining_vcf_multiallelic"
#> [2] "repair_vcf_remaining_vcf_duplicated_marker_ids"
#> [3] "repair_vcf_remaining_vcf_duplicated_chrom_pos"
file.exists(out)
#> [1] FALSE
```

## Chromosome-style normalization

If your VCF mixes `chr1` and `1` for the same chromosome, the default
`chrom_style = "auto"` keeps them separate. Pass `"numeric"` or
`"chr_prefixed"` to normalize labels in the output **and** use the
normalized form for duplicate detection:

```r
# chr1 + 1 with the same POS → seen as duplicates under "numeric"
repair_vcf("messy.vcf", "out.vcf", chrom_style = "numeric")
```

## File-level blockers (no repair)

Some problems can't be safely repaired. They show up as `$errors`
and the function writes nothing:

- `repair_vcf_missing_chrom_header` — no `#CHROM` line.
- `repair_vcf_no_samples` — Beagle needs samples.
- `repair_vcf_duplicated_samples` — no obvious resolution.

Fix these by hand (or upstream) before retrying.
