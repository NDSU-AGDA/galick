# Benchmarking

Three functions for measuring imputation accuracy on your own data:
mask known genotypes, run imputation, compare to the truth.

## `mask_genotypes()`

```r
mask_genotypes(
  vcf,
  output_vcf,
  fraction  = 0.1,
  seed      = NULL,
  samples   = NULL,
  markers   = NULL,
  overwrite = FALSE
) -> beagler_result
```

Reads `vcf`, randomly sets `fraction` of the **non-missing** GT cells
to `./.`, and writes the result to `output_vcf` (a minimal VCF v4.2).

```r
result <- mask_genotypes("data.vcf", "data.masked.vcf",
                         fraction = 0.1, seed = 42L)
result$summary$n_eligible    # number of non-missing cells
result$summary$n_masked      # number actually masked
```

Reproducible: same `seed` → same masked VCF.

`samples` and `markers` restrict which cells are eligible to be masked.

## `compare_imputed_vs_true()`

```r
compare_imputed_vs_true(
  imputed_vcf,
  true_vcf,
  masked_vcf = NULL
) -> beagler_result
```

Computes per-cell accuracy on shared `(sample, marker)` cells:

| metric | meaning |
|---|---|
| `n_compared` | cells with non-missing in both files |
| `n_concordant` | cells where rounded imputed dosage == true dosage |
| `concordance_rate` | `n_concordant / n_compared` |
| `dosage_r2` | squared Pearson correlation of dosages |
| `rmse` | √mean((imputed − true)²) on dosages |

If `masked_vcf` is provided, the comparison is restricted to cells
that were `./.` in masked AND non-missing in truth — i.e., the cells
imputation actually had to fill in.

```r
result <- compare_imputed_vs_true(
  imputed_vcf = "imputed.vcf.gz",
  true_vcf    = "data.vcf",
  masked_vcf  = "data.masked.vcf"
)
result$summary$concordance_rate
#> [1] 0.92
result$summary$dosage_r2
#> [1] 0.95
```

## `benchmark_imputation()`

```r
benchmark_imputation(
  vcf,
  output_dir,
  fraction = 0.1,
  seed     = NULL,
  mode     = c("impute", "phase"),
  ...
) -> list(mask, run, comparison)
```

Convenience orchestrator:

1. Mask `vcf` via `mask_genotypes()` →
   `<output_dir>/masked.vcf`.
2. Call `run_beagle()` on the masked file →
   `<output_dir>/imputed.vcf.gz`. `...` is forwarded (e.g. `nthreads`,
   `memory_gb`, `ref`, `map`).
3. Compare the imputed output to the original (truth) via
   `compare_imputed_vs_true()`, restricted to masked cells.

```r
b <- benchmark_imputation(
  vcf = "data.vcf",
  output_dir = "bench/",
  fraction = 0.2,
  seed = 42L,
  mode = "phase",
  nthreads = 4L
)

b$mask$summary$n_masked
b$run$status
b$comparison$summary$concordance_rate
```

Returns a list with three slots:

- `$mask` — `beagler_result` from `mask_genotypes()`.
- `$run` — `beagler_run` from `run_beagle()`.
- `$comparison` — `beagler_result` from `compare_imputed_vs_true()`,
  with the masked-cells filter applied.

Beagle has to actually run, so this requires a working `BEAGLE_PATH`.
