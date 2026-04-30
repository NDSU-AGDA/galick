# Chromosome workflows

Three exports for splitting a VCF, running Beagle in parallel across
chromosomes, and merging the results back together.

## `split_vcf_by_chrom()`

```r
split_vcf_by_chrom(
  vcf,
  out_dir,
  chrom_filter = NULL,
  bcftools     = NULL,
  overwrite    = FALSE
) -> list(output_dir, chrom_paths)
```

Reads `vcf` (`.vcf` or `.vcf.gz`) and writes one `.vcf` per unique
CHROM into `out_dir`. Pure R; no `bcftools` dependency.

```r
res <- split_vcf_by_chrom("data.vcf.gz", "split_out")
res$chrom_paths
#>      chr1                     chr2                     chr3
#> "split_out/chr1.vcf"     "split_out/chr2.vcf"     "split_out/chr3.vcf"
```

Use `chrom_filter` to keep only listed chromosomes; the rest are
dropped from the output.

## `run_beagle_by_chrom()`

```r
run_beagle_by_chrom(
  gt,
  out_dir,
  ref          = NULL,
  map          = NULL,
  chrom_filter = NULL,
  parallel     = c("none", "future", "mclapply"),
  workers      = 1L,
  ...
) -> list of beagler_run, named by chromosome
```

Splits `gt` under `out_dir/_split/`, then calls `run_beagle()` once
per chromosome. `ref` and `map` are passed to each invocation as-is —
Beagle filters them by chromosome internally.

### Resume logic

If a previous run left a `<out_dir>/<chrom>/run_record.json` with
`status == 0` AND its recorded input SHA-256 matches the current
per-chromosome `gt`, that chromosome is skipped and the cached
`beagler_run` is returned. The hash check guards against silent data
churn (someone edits the input but the run record is still around).

```r
# First call: runs all chromosomes
results <- run_beagle_by_chrom("panel.vcf", "imputed",
                               nthreads = 4L)

# Modify one chromosome's input, then re-run:
results <- run_beagle_by_chrom("panel.vcf", "imputed",
                               nthreads = 4L)
# Only the changed chromosomes re-execute.
```

### Parallel execution

```r
# Multicore via base::parallel
run_beagle_by_chrom(gt, out_dir,
                    parallel = "mclapply", workers = 4L)

# future.apply (declared in Suggests)
future::plan(future::multisession, workers = 4)
run_beagle_by_chrom(gt, out_dir, parallel = "future")
```

`parallel = "future"` requires the `future.apply` package; install it
with `install.packages("future.apply")`.

## `merge_beagle_outputs()`

```r
merge_beagle_outputs(
  run_results,
  output_vcf,
  bcftools  = NULL,
  overwrite = FALSE
) -> path (invisible)
```

Concatenates the `output_paths` from a list of `beagler_run` objects
into a single `.vcf` or `.vcf.gz`. The first non-empty input
contributes the header; subsequent inputs contribute body rows only.

```r
results <- run_beagle_by_chrom("panel.vcf", "imputed", nthreads = 4L)

merged <- "imputed_all.vcf.gz"
merge_beagle_outputs(results, merged)

# Verify
validate_vcf(merged, target = "beagle")$ok
```

Pure R (uses `gzfile()` for `.vcf.gz` IO). A `bcftools concat` fast
path is forward-compatible (the `bcftools` argument is currently
ignored).
