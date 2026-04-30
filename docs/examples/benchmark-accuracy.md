# Benchmark imputation accuracy

Mask a fraction of known genotypes, run imputation, and measure how
well Beagle recovers them.

!!! note "Beagle required"

    The orchestrator runs Beagle internally. `mask_genotypes()` and
    `compare_imputed_vs_true()` work without Beagle if you bring your
    own imputed VCF.

## Quick benchmark

```r
library(beagler)

# Build a starter VCF
csv <- system.file("extdata", "tiny.csv", package = "beagler")
map <- system.file("extdata", "tiny_map.tsv", package = "beagler")
truth_vcf <- tempfile(fileext = ".vcf")
csv_to_vcf(csv, map, truth_vcf)

# Run the full pipeline
b_dir <- tempfile("bench_")
dir.create(b_dir)
b <- benchmark_imputation(
  vcf = truth_vcf,
  output_dir = b_dir,
  fraction = 0.2,
  seed = 42L,
  mode = "phase",
  nthreads = 2L
)

b$mask$summary$n_masked
#> [1] 20

b$run$status
#> [1] 0

b$comparison$summary$concordance_rate
#> [1] 0.9
b$comparison$summary$dosage_r2
#> [1] 0.93
b$comparison$summary$rmse
#> [1] 0.31
```

## Step by step

If you want more control over each stage:

### 1. Mask

```r
masked_vcf <- file.path(b_dir, "masked.vcf")
mask <- mask_genotypes(
  vcf = truth_vcf,
  output_vcf = masked_vcf,
  fraction = 0.2,
  seed = 42L,
  overwrite = TRUE
)

mask$summary$n_eligible
#> [1] 100      # non-missing cells eligible for masking
mask$summary$n_masked
#> [1] 20       # round(fraction * n_eligible)
```

`seed` makes the masking reproducible — same seed, same masked file.

You can restrict eligibility to specific samples or markers:

```r
mask_genotypes(truth_vcf, masked_vcf, fraction = 0.5, seed = 1L,
               samples = c("s1", "s2"),     # only mask these
               markers = c("rs1", "rs5", "rs10"))
```

### 2. Run Beagle

```r
imputed_prefix <- file.path(b_dir, "imputed")
run <- run_beagle(masked_vcf, imputed_prefix, mode = "phase",
                  nthreads = 2L, overwrite = TRUE)
imputed_vcf <- run$output_paths[1]   # imputed.vcf.gz
```

### 3. Compare

```r
result <- compare_imputed_vs_true(
  imputed_vcf = imputed_vcf,
  true_vcf    = truth_vcf,
  masked_vcf  = masked_vcf      # restrict to cells we hid
)

result$summary
```

```text
$n_compared           20
$n_concordant         18
$concordance_rate      0.9
$dosage_r2             0.93
$rmse                  0.31
$n_samples_compared    5
$n_markers_compared   20
$used_mask_filter   TRUE
```

When `masked_vcf` is provided, the comparison restricts to cells that
were `./.` in masked AND non-missing in truth — i.e., the cells
imputation actually had to fill in. Without it, all overlap cells
contribute (which makes the numbers look better than they should).

## Reproducibility

The whole pipeline is deterministic given a `seed`:

```r
b1 <- benchmark_imputation(truth_vcf, tempfile("b1_"),
                           fraction = 0.1, seed = 7L, mode = "phase")
b2 <- benchmark_imputation(truth_vcf, tempfile("b2_"),
                           fraction = 0.1, seed = 7L, mode = "phase")

identical(b1$mask$summary, b2$mask$summary)
#> [1] TRUE
identical(b1$comparison$summary$concordance_rate,
          b2$comparison$summary$concordance_rate)
#> [1] TRUE
```

(Beagle itself uses its own RNG with its own `seed=` argument; for
fully deterministic Beagle output you'd also pass `seed = ` to
`benchmark_imputation()` via `...`. Without it, Beagle picks a
random seed and the imputed VCF can vary slightly between runs.)
