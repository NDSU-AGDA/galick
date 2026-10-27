# Benchmarking

Functions for measuring imputation accuracy on your own data:
mask known genotypes, run imputation, compare to the truth.

## `design_marker_panel()`

Select an exact number of diploid SNPs with mandatory/excluded sites and minimum
spacing. Supply disjoint candidate and scoring identity tables (`chrom`, `pos`,
`ref`, `alt`), trusted genotypes, a fully called phased reference with no truth
sample IDs, and sample partitions `design`, `tune`, `test`.

Design-only local LD proposes panels; bounded Beagle runs measure recovery.
Scoring sites never enter any panel, and all comparisons use the same original
called-genotype masks. A challenger is frozen before tuning; the tuning choice
is frozen before test confirmation. No engine parameters are optimized.

Every group requires 20 design and 10 tuning/test samples. Declared blocks cannot
cross partitions. Candidate eligibility uses design-only call rate >=95 percent
and pooled MAF >=1 percent. Reference relatedness and renamed duplicates are not
inferred. Without a genetic map Beagle assumes 1 cM/Mb; sparse observed sites can
fail Beagle's per-window requirements.

`(max_proposals + 6) * repeats` jobs must fit `max_runs` before launching. Repeats
change paired engine seeds, not biological samples. Cell and local-pair ceilings
are workload bounds, not total RAM limits. No global optimum or assay feasibility
is claimed. `recommended_panel` is NULL unless test confirmation succeeds.
Inspect metrics, group scores, search history and hashed `panel_design.json`.
See `vignette("design-marker-panel", package = "galick")` for a runnable example.

The development fixture `tools/profile-panel-design.R` records performance, not
biological accuracy. On an arm64 development machine with R 4.5.3 and Beagle 5.4,
120 truth samples, 120 reference samples, 160 candidate sites, 80 scoring sites,
20 selected markers, two proposals and one repeat completed in 16.6 seconds
(seven jobs); the R process reached approximately 206 MiB resident memory.
JVM heap was capped at 1 GiB; the R resident-memory figure is not aggregate
R-plus-Java peak memory. Genome-scale studies require separate profiling.

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

The `galick_evaluation` result inherits from `galick_result`. Decisions are
`confirmed`, `no_clear_winner`, `not_confirmed`, or `incomplete`. Only confirmed
results have a recommended candidate. Each output directory includes isolated
jobs, hashed masks and inputs, metrics and versioned `evaluation.json`.
Nonempty output directories are refused. See the installed
`vignette("evaluate-imputation", package = "galick")` for a runnable example.

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
)
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
)
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
result$summary$dosage_r2
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
)
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

- `$mask` — `galick_result` from `mask_genotypes()`.
- `$run` — `galick_run` from `run_beagle()`.
- `$comparison` — `galick_result` from `compare_imputed_vs_true()`,
  with the masked-cells filter applied.

Beagle has to actually run, so this requires a working `BEAGLE_PATH`.
