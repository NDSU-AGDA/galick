---
title: Run chromosome-parallel Beagle imputation in R
description: Configure per-chromosome Beagle jobs, explicit parallel resources, resumable run records and validated output merging in Galick.
---
# Chromosome-parallel imputation

Split a multi-chromosome VCF, run Beagle on each chromosome in
parallel, and merge the outputs.

!!! note "Beagle required"

    Like the [Run Beagle](run-beagle.md) example, this needs Beagle
    on `BEAGLE_PATH` or in the user cache.

## Prepare a two-chromosome VCF

The bundled fixture has five samples and twenty markers across two chromosomes.
It demonstrates the workflow, not biological imputation performance.

```r
library(galick)

vcf <- tempfile(fileext = ".vcf")
conversion <- csv_to_vcf(
  system.file("extdata", "tiny.csv", package = "galick"),
  system.file("extdata", "tiny_map.tsv", package = "galick"),
  vcf, encoding = "AHB"
)
stopifnot(conversion$ok)
```

## Just split

If you only want per-chromosome files (e.g. to hand off to a
different tool), use `split_vcf_by_chrom()` directly:

```r
split_dir <- tempfile("galick_split_")
res <- split_vcf_by_chrom(vcf, split_dir)
res$chrom_paths
```

`chrom_filter` keeps a subset:

```r
split_vcf_by_chrom(vcf, split_dir,
                   chrom_filter = "chr1",
                   overwrite = TRUE)
```

## Run Beagle per chromosome

`run_beagle_by_chrom()` splits internally, then calls `run_beagle()`
once per chromosome:

Choose one of the following modes. Each tab creates its own output directory.
Use `mclapply` only on Unix-like systems. The future option requires
`future.apply`; its multisession plan also works on Windows.

=== "Sequential"

    ```r
    out_dir <- tempfile("galick_imp_")
    dir.create(out_dir)
    
    results <- run_beagle_by_chrom(
      gt = vcf,
      out_dir = out_dir,
      mode = "phase",
      nthreads = 1L,
      memory_gb = 1L
    )
    ```

=== "mclapply (UNIX)"

    ```r
    out_dir <- tempfile("galick_multicore_")
    results <- run_beagle_by_chrom(
      gt = vcf,
      out_dir = out_dir,
      mode = "phase",
      nthreads = 1L,
      memory_gb = 1L,
      parallel = "mclapply",
      workers = 2L
    )
    ```

=== "future"

    ```r
    out_dir <- tempfile("galick_future_")
    previous_plan <- future::plan(future::multisession, workers = 2)
    
    results <- run_beagle_by_chrom(
      gt = vcf,
      out_dir = out_dir,
      mode = "phase",
      nthreads = 1L,
      memory_gb = 1L,
      parallel = "future"
    )
    future::plan(previous_plan)
    ```

The result is a named list of `galick_run` objects:

```r
names(results)
stopifnot(all(vapply(results, function(run) identical(run$status, 0L), logical(1))))
results$chr1$status
results$chr2$output_paths
```

## Resume after a failure

If the run was interrupted (or one chromosome failed and you fixed
the input), call `run_beagle_by_chrom()` again with the same
`out_dir`. Successful chromosomes are skipped:

```r
# With unchanged successful inputs, both chromosomes reuse their results.
results <- run_beagle_by_chrom(vcf, out_dir, mode = "phase",
                               nthreads = 1L, memory_gb = 1L)
# If only chr2's input or settings changed, only chr2 would run again.
```

Resume requires success, a matching complete job fingerprint (inputs, reference,
map, settings and runtime), an output checksum match and valid VCF content.
Changing any computational input invalidates the cached chromosome.

## Merge

```r
merged <- tempfile(fileext = ".vcf")
merge_beagle_outputs(results, merged)

validate_vcf(merged, target = "beagle")$ok

vcf_to_dosage(merged) |> dim()
```

`merge_beagle_outputs()` accepts `.vcf` or `.vcf.gz` for the output
extension. Compressed output requires `bgzip`; no `bcftools` dependency.
All expected jobs must have succeeded, with matching sample sets and compatible
headers. Sample order is aligned and DS/GP fields are retained.
