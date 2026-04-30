# Quick tour

A 5-minute walkthrough that uses the bundled fixtures, skips Beagle
when it isn't available, and ends with a dosage matrix you can pass to
downstream tools.

## Bundled fixtures

beagler ships three tiny example files in `inst/extdata/`:

```r
library(beagler)

csv_path <- system.file("extdata", "tiny.csv", package = "beagler")
map_path <- system.file("extdata", "tiny_map.tsv", package = "beagler")
hmp_path <- system.file("extdata", "tiny.hmp.txt", package = "beagler")

read.csv(csv_path)[1:3, 1:5]
#>   sample_id rs1 rs2 rs3 rs4
#> 1        s1   B   A   B   H
#> 2        s2   A   B   H   B
#> 3        s3   H   B   A   B
```

Five samples × 20 markers, AHB encoded.

## CSV to VCF

```r
out <- tempfile(fileext = ".vcf")
result <- csv_to_vcf(csv_path, map_path, out)
result$ok
#> [1] TRUE
result$summary$inferred_encoding
#> [1] "AHB"
```

If anything goes wrong, no VCF is written and `result$errors` carries
classed codes you can dispatch on.

## Validate

`validate_vcf()` reports issues against three targets:

```r
validate_vcf(out, target = "general")$ok    # TRUE
validate_vcf(out, target = "beagle")$ok     # TRUE
validate_vcf(out, target = "plink")$ok      # TRUE
```

`target = "beagle"` adds Beagle-specific checks: diploid throughout,
no multiallelic sites, REF ≠ ALT, BGZF or warning.

## Repair (if needed)

```r
dirty <- system.file("extdata", "tiny_dirty.vcf", package = "beagler")
v <- validate_vcf(dirty, target = "beagle")
v$ok
#> [1] FALSE
v$errors$code
#>  [1] "vcf_invalid_pos"           "vcf_invalid_chrom"
#>  [3] "vcf_duplicated_marker_ids" "vcf_duplicated_chrom_pos"
#>  ...

clean <- tempfile(fileext = ".vcf")
repair_vcf(dirty, clean)$ok
#> [1] TRUE

validate_vcf(clean, target = "beagle")$ok
#> [1] TRUE
```

`repair_vcf()` runs five steps in order and re-validates the result
before writing — if the requested steps couldn't make the output
pass, nothing lands and the result reports `repair_vcf_remaining_*`
errors.

## Run Beagle

The `nzchar()` guard makes this chunk a no-op when Beagle isn't
configured:

```r
if (nzchar(Sys.getenv("BEAGLE_PATH"))) {
  out_dir <- tempfile("tour_")
  dir.create(out_dir)
  out_prefix <- file.path(out_dir, "phased")

  run <- run_beagle(out, out_prefix, mode = "phase", nthreads = 2L)
  run$status
  #> [1] 0
  run$output_paths
  #> [1] ".../phased.vcf.gz"
}
```

`run` is a `beagler_run` object with status, runtime, the full argv,
log paths, and a JSON run record (with SHA-256 hashes of the inputs
and outputs).

## Pull results back

```r
imputed <- if (nzchar(Sys.getenv("BEAGLE_PATH"))) run$output_paths[1] else out

dosage <- vcf_to_dosage(imputed)
dim(dosage)            # samples × markers
#> [1]  5 20
dosage[1:3, 1:4]
#>    rs1 rs2 rs3 rs4
#> s1   2   0   2   1
#> s2   0   2   1   2
#> s3   1   2   0   2

hard <- vcf_to_hardcalls(imputed, encoding = "AHB")
hard[1:3, 1:4]
#>    rs1 rs2 rs3 rs4
#> s1   B   A   B   H
#> s2   A   B   H   B
#> s3   H   B   A   B
```

## What's next

- The [result classes](result-classes.md) page explains the shape of
  `beagler_result` and `beagler_run`.
- [Reference](../reference/index.md) lists every export grouped by
  area.
- [Examples](../examples/index.md) walk through end-to-end workflows
  including chromosome-parallel imputation and accuracy benchmarking.
