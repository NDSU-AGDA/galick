---
title: Compare Beagle configurations on held-out genotypes
description: Set up paired genotype masking, separate tuning and test samples, and confirm or reject an imputation configuration without test-set reranking.
---
# Held-out imputation evaluation

This worked example is also available as an installed vignette. Small constructed
inputs demonstrate the API; they are not evidence of biological performance.


```r
library(galick)
```

`evaluate_imputation()` compares two to four Beagle configurations on paired
held-out genotypes. Tuning samples select a configuration; separate test samples
confirm it against the baseline. A tie or failed confirmation produces no
recommendation. This is called-genotype evaluation, not read-aware inference.

## Define the experiment

The bundled data are small examples, not evidence of biological performance.
Use independently verified calls and deployment-appropriate sample partitions
for a real study. The function does not infer relatedness or find renamed duplicates.

```r
truth <- tempfile(fileext = ".vcf")
csv_to_vcf(
  system.file("extdata", "tiny.csv", package = "galick"),
  system.file("extdata", "tiny_map.tsv", package = "galick"),
  truth
)
samples <- data.frame(
  sample_id = paste0("s", 1:5),
  partition = c("tune", "tune", "test", "test", "test"),
  group = "example"
)
candidates <- list(
  baseline = list(extra_args = "iterations=4"),
  longer = list(extra_args = "iterations=6")
)
```

An optional `block` column keeps declared families or lines within a partition.
The `group` column defines reporting strata, not automatically independent
holdout blocks. Each partition/group needs at least one held-out called genotype.

## Run a paired comparison

This section requires Java and a configured Beagle JAR. It launches real Beagle
jobs and does not download software. Run the blocks in order.

```r
study <- evaluate_imputation(
  truth, samples, candidates, baseline = "baseline",
  output_dir = tempfile("evaluation_"),
  fraction = 0.2, repeats = 1L, seed = 42L,
  field = "GT", nthreads = 1L, memory_gb = 1L,
  min_improvement = 0.01, max_group_regression = 0.01
)
study
study$selection
study$by_group
```

The worst-case budget is `(number of candidates + 2) * repeats` jobs. Jobs run
sequentially. All candidates receive the same partition-specific input and seed
within a repeat. Test samples never enter tuning jobs. Caller RNG state is restored.

## Selection and failed confirmation

Selection requires complete paired predictions, the declared absolute RMSE gain
over the baseline and runner-up, and no group RMSE regression beyond the limit.
RMSE pools squared errors and prediction counts, not group averages. Missing
predictions remain in the denominator for completion and GT concordance, and
prevent selection. DS never falls back to GT. Undefined per-marker R-squared
values remain explicit rather than becoming zero.

This deterministic arithmetic example illustrates failed confirmation; it is
not output from a biological benchmark. Assume all group checks pass:

```r
example <- data.frame(
  partition = c("tune", "test"),
  baseline_rmse = c(0.40, 0.40),
  candidate_rmse = c(0.20, 0.50)
)
example$improvement <- example$baseline_rmse - example$candidate_rmse
frozen_candidate <- "longer" # Selected using tuning only.
decision <- if (example$improvement[example$partition == "test"] >= 0.01) {
  "confirmed"
} else {
  "not_confirmed"
}
data.frame(frozen_candidate, decision)
```

Do not choose a different configuration after seeing test outcomes. A tuning
tie runs only the baseline on test samples and reports `no_clear_winner`.
The margins are practical thresholds, not confidence intervals or significance tests.

## Evaluate an actual marker panel

For array-to-dense evaluation, every candidate needs a fully called, phased
reference with no evaluation sample IDs. Supply exact marker identities:

The following block is a template for your own data, not a continuation using
the bundled fixture. Replace `truth`, `samples`, the panel coordinates and both
reference paths with compatible inputs before running it.

```r
panel <- data.frame(chrom = "chr1", pos = c(100, 2100), ref = c("A", "G"), alt = c("T", "A"))
study <- evaluate_imputation(
  truth, samples,
  candidates = list(
    panel_a = list(ref = "reference_a.vcf.gz"),
    panel_b = list(ref = "reference_b.vcf.gz")
  ),
  baseline = "panel_a", output_dir = "array_evaluation",
  design = "array", observed_markers = panel, field = "DS",
  repeats = 3L, seed = 42L
)
```

Include observed markers for every evaluated chromosome. Untyped records are
physically absent from engine inputs, not just blanked. Missing reference sites
remain unpredicted truth in the evaluation; they are never silently excluded.
Array repeats use the same panel with different paired engine seeds.

## Audit and limits

The output directory contains hashed held-out manifests and GT-only inputs,
isolated job logs/run records, per-marker tables, `metrics.tsv`, `by_group.tsv`,
the frozen `selection.json`, and versioned `evaluation.json`. Source files are
not modified. Existing nonempty directories are refused; there is no automatic resume.

The returned object inherits from `galick_result`; `$ok` reports execution/scoring
errors, not proof of superiority. Inspect `$selection$decision`, completion and
group results. A result with no winner can still have `$ok = TRUE`.

Input must be diploid biallelic SNP calls. No pedigree inference, automatic
parameter search or statistical independence claim is made. Use `field = "GT"`
when the engine does not produce expected DS for the evaluated records. Do not
interpret model-estimated DR2 as measured accuracy on held-out truth.
