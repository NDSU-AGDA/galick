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

## Build a multi-chrom VCF

```r
library(galick)

vcf <- tempfile(fileext = ".vcf")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=chr1>",
  "##contig=<ID=chr2>",
  "##contig=<ID=chr3>",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
  "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
  "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
  "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
  "chr2\t50\trs4\tT\tA\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1",
  "chr2\t100\trs5\tA\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1",
  "chr3\t75\trs6\tG\tT\t.\tPASS\t.\tGT\t0/1\t0/1\t1/1"
), vcf)
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
                   chrom_filter = c("chr1", "chr3"),
                   overwrite = TRUE)
```

## Run Beagle per chromosome

`run_beagle_by_chrom()` splits internally, then calls `run_beagle()`
once per chromosome:

=== "Sequential"

    ```r
    out_dir <- tempfile("galick_imp_")
    dir.create(out_dir)
    
    results <- run_beagle_by_chrom(
      gt = vcf,
      out_dir = out_dir,
      mode = "phase",
      nthreads = 2L
    )
    ```

=== "mclapply (UNIX)"

    ```r
    results <- run_beagle_by_chrom(
      gt = vcf,
      out_dir = out_dir,
      mode = "phase",
      nthreads = 2L,
      parallel = "mclapply",
      workers = 4L
    )
    ```

=== "future"

    ```r
    library(future)
    plan(multisession, workers = 4)
    
    results <- run_beagle_by_chrom(
      gt = vcf,
      out_dir = out_dir,
      mode = "phase",
      nthreads = 2L,
      parallel = "future"
    )
    ```

The result is a named list of `galick_run` objects:

```r
names(results)
results$chr1$status
results$chr2$output_paths
```

## Resume after a failure

If the run was interrupted (or one chromosome failed and you fixed
the input), call `run_beagle_by_chrom()` again with the same
`out_dir`. Successful chromosomes are skipped:

```r
# Suppose chr2 failed; you fixed it and rerun:
results <- run_beagle_by_chrom(vcf, out_dir, mode = "phase",
                               nthreads = 2L)
# chr1 and chr3 reuse their previous run records (status == 0,
# complete job fingerprint and output checksum match). Only chr2 executes.
```

Resume requires success, a matching complete job fingerprint (inputs, reference,
map, settings and runtime), an output checksum match and valid VCF content.
Changing any computational input invalidates the cached chromosome.

## Merge

```r
merged <- tempfile(fileext = ".vcf.gz")
merge_beagle_outputs(results, merged)

validate_vcf(merged, target = "beagle")$ok

vcf_to_dosage(merged) |> dim()
```

`merge_beagle_outputs()` accepts `.vcf` or `.vcf.gz` for the output
extension. Compressed output requires `bgzip`; no `bcftools` dependency.
All expected jobs must have succeeded, with matching sample sets and compatible
headers. Sample order is aligned and DS/GP fields are retained.
