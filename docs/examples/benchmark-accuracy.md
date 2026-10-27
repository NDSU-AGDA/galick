# Benchmark imputation accuracy

This instructional example masks known GT calls, runs Beagle and compares the
recovered calls with the originals. It needs Java and Beagle. Its five samples
and 20 markers are far too small to establish biological accuracy.

## Prepare and benchmark

```r
library(galick)

csv <- system.file("extdata", "tiny.csv", package = "galick")
map <- system.file("extdata", "tiny_map.tsv", package = "galick")
truth <- tempfile(fileext = ".vcf")
conversion <- csv_to_vcf(csv, map, truth, encoding = "AHB")
stopifnot(conversion$ok)

benchmark <- benchmark_imputation(
  vcf = truth,
  output_dir = tempfile("benchmark_"),
  fraction = 0.2,
  seed = 42L,
  mode = "phase",
  field = "GT",
  nthreads = 1L
)
stopifnot(identical(benchmark$run$status, 0L))
benchmark$comparison$summary
```

The comparison is restricted to calls masked from nonmissing truth. Inspect
`n_eligible`, `n_missing_predictions`, `completion_rate`, concordance and RMSE.
Undefined per-marker correlations remain undefined; they are not evidence of
perfect recovery. The package does not guarantee particular scores.

## Run each step explicitly

```r
masked <- tempfile(fileext = ".vcf")
mask <- mask_genotypes(
  vcf = truth,
  output_vcf = masked,
  fraction = 0.2,
  seed = 42L
)
stopifnot(mask$ok)

run <- run_beagle(
  gt = masked,
  out = tempfile("benchmark_imputed_"),
  mode = "phase",
  seed = 42L,
  nthreads = 1L
)
stopifnot(identical(run$status, 0L))

comparison <- compare_imputed_vs_true(
  imputed_vcf = run$output_paths[[1]],
  true_vcf = truth,
  masked_vcf = masked,
  field = "GT"
)
comparison$summary
comparison$summary$per_marker
```

`field = "DS"` requests expected ALT dosage explicitly and never falls back to
GT. Use it only when the imputed output contains DS for the records being scored.

## Design a real evaluation

Random-cell masking measures recovery under that missingness experiment. It
does not establish transfer to a different family, generation, population or
low-density assay. Do not use evaluation samples in a reference panel, and do
not treat previously imputed calls as independently verified truth.

Use [Held-out evaluation](evaluate-imputation.md) to select a configuration on
tuning samples and confirm it on separate test samples. For array deployment,
physically thin the input to observed markers and supply an independent phased
reference. Repeat engine seeds are not independent biological replications.
