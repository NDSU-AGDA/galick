# Benchmarking

Functions for measuring imputation accuracy on your own data:
mask known genotypes, run imputation, compare to the truth.

## `evaluate_imputation()`

Compare two to four named configurations using separate tune/test sample
partitions. Candidate lists accept `ref`, `map` and supported `extra_args`.
The baseline must be named explicitly. Use sporadic cell masking or an actual
array marker panel; array evaluation physically removes untyped input records
and requires a reference for every configuration.

```r
study <- evaluate_imputation(
  truth_vcf = "truth.vcf.gz", sample_data = samples,
  candidates = list(default = list(), longer = list(extra_args = "iterations=20")),
  baseline = "default", output_dir = "evaluation", repeats = 3L, seed = 42L
)
study$selection
study$metrics
study$by_group
```

`sample_data` requires unique `sample_id` and `partition` (tune/test), with
optional reporting `group` and holdout `block`. Blocks cannot cross partitions.
RMSE selection requires complete predictions, the declared minimum improvement
over baseline and runner-up, and group regression safeguards. The selected
candidate is frozen before test scoring. Test outcomes never rerank candidates.
Missing DS never falls back to GT. Correlation is descriptive, not the selection
objective. Practical margins are not significance tests.

The `beagler_evaluation` result inherits from `beagler_result`. Decisions are
`confirmed`, `no_clear_winner`, `not_confirmed`, or `incomplete`. Only confirmed
results have a recommended candidate. Each output directory includes isolated
jobs, hashed masks and inputs, metrics and versioned `evaluation.json`.
Nonempty output directories are refused. See the installed
`vignette("evaluate-imputation", package = "beagler")` for a runnable example.

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
  masked_vcf = NULL,
  field = c("GT", "DS")
) -> beagler_result
```

Variants must match CHROM, POS, REF and ALT exactly; allele flips are never
guessed. Evaluation is diploid and biallelic. `field = "DS"` selects expected
ALT dosage explicitly and never falls back to GT. Truth is always GT.
Missing predictions and missing output samples/variants remain in the
eligible-cell denominator:

| metric | meaning |
|---|---|
| `n_eligible` | nonmissing truth cells selected for evaluation |
| `n_compared` | eligible cells with the requested prediction field |
| `n_missing_predictions` | eligible cells without a prediction |
| `completion_rate` | `n_compared / n_eligible` |
| `n_concordant` | eligible cells with matching GT counts; DS is never rounded |
| `concordance_rate` | `n_concordant / n_eligible` |
| `conditional_concordance_rate` | GT concordance among predicted GT cells |
| `dosage_r2` | mean defined per-marker squared correlation |
| `pooled_dosage_r2` | pooled squared correlation, not a substitute for per-marker accuracy |
| `n_undefined_r2` | markers with too few predictions or zero variance |
| `rmse` | √mean((imputed − true)²) on dosages |

If `masked_vcf` is provided, the comparison is restricted to cells
that were `./.` in masked AND non-missing in truth — i.e., the cells
imputation actually had to fill in. The mask must contain every truth variant
and sample. Undefined correlations remain `NA`; inspect `per_marker` and
completion alongside accuracy. Random-cell masking is not a deployment
validation design for population, family or array-transfer performance.

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
  field    = c("GT", "DS"),
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

The seed is passed to both masking and Beagle. The original file is treated
as truth, not independently verified truth. Keep evaluation targets out of
reference panels to avoid leakage.

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
