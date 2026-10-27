# API signatures

Generated from the package source by `tools/build-reference.R`.
Use `?function_name` in R for complete parameter documentation.
Signatures below are declarations; required inputs are not example data.

## `beagle_version()`

```r
beagle_version(beagle_jar = NULL)
```

## `benchmark_imputation()`

```r
benchmark_imputation(
  vcf,
  output_dir,
  fraction = 0.1,
  seed = NULL,
  mode = c("impute", "phase"),
  field = c("GT", "DS"),
  ...
)
```

## `compare_imputed_vs_true()`

```r
compare_imputed_vs_true(
  imputed_vcf,
  true_vcf,
  masked_vcf = NULL,
  field = c("GT", "DS")
)
```

## `csv_to_vcf()`

```r
csv_to_vcf(
  csv,
  marker_map,
  output_vcf,
  sample_col = NULL,
  encoding = c("auto", "AHB", "012", "centered", "AABB", "VCF", "REFHET", "IUPAC"),
  ploidy = 2L,
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
  ahb_ref = c("A", "B"),
  aabb_ref_homo = c("AA", "BB"),
  missing_codes = c(NA, ".", "--", "", "N", "NN", "./."),
  overwrite = FALSE
)
```

## `design_marker_panel()`

```r
design_marker_panel(
  truth_vcf,
  reference_vcf,
  sample_data,
  candidate_markers,
  scoring_markers,
  n_markers,
  output_dir,
  mandatory_markers = NULL,
  excluded_markers = NULL,
  min_spacing_bp = 0L,
  ld_window_bp = 1e+06,
  map = NULL,
  extra_args = character(),
  field = c("DS", "GT"),
  repeats = 3L,
  max_proposals = 12L,
  max_runs = 60L,
  min_improvement = 0.01,
  max_group_regression = 0.01,
  max_genotype_cells = 2e+07,
  max_ld_pairs = 2e+06,
  seed = 2026L,
  nthreads = 1L,
  memory_gb = 4L,
  beagle_jar = NULL,
  java = NULL
)
```

## `detect_genotype_ploidy()`

```r
detect_genotype_ploidy(vcf)
```

## `detect_subgenome_labels()`

```r
detect_subgenome_labels(vcf)
```

## `evaluate_imputation()`

```r
evaluate_imputation(
  truth_vcf,
  sample_data,
  candidates,
  baseline,
  output_dir,
  design = c("sporadic", "array"),
  observed_markers = NULL,
  fraction = 0.1,
  field = c("GT", "DS"),
  repeats = 3L,
  seed = 2026L,
  min_improvement = 0.01,
  max_group_regression = 0.01,
  nthreads = 1L,
  memory_gb = 4L,
  max_runs = 20L,
  beagle_jar = NULL,
  java = NULL
)
```

## `export_genotypes()`

```r
export_genotypes(
  vcf,
  output,
  format = c("vcf", "plink", "hapmap", "csv", "tsv", "dosage"),
  downcode = NULL,
  overwrite = FALSE,
  digits = 17L,
  field = c("GT", "DS")
)
```

## `hapmap_to_vcf()`

```r
hapmap_to_vcf(
  hmp,
  output_vcf,
  build = NA_character_,
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
  missing_codes = c(NA, ".", "--", "", "N", "NN"),
  overwrite = FALSE
)
```

## `imputation_report()`

```r
imputation_report(run, file)
```

## `install_beagle()`

```r
install_beagle(version = "5.4", dest = NULL)
```

## `mask_genotypes()`

```r
mask_genotypes(
  vcf,
  output_vcf,
  fraction = 0.1,
  seed = NULL,
  samples = NULL,
  markers = NULL,
  overwrite = FALSE
)
```

## `merge_beagle_outputs()`

```r
merge_beagle_outputs(
  run_results,
  output_vcf,
  bcftools = NULL,
  overwrite = FALSE
)
```

## `merge_vcf()`

```r
merge_vcf(
  vcfs,
  output_vcf,
  mode = c("by_sample", "by_marker"),
  sample_collision = c("error", "suffix", "coalesce"),
  overwrite = FALSE
)
```

## `plink_to_vcf()`

```r
plink_to_vcf(
  plink_prefix,
  vcf,
  plink = NULL,
  overwrite = FALSE
)
```

## `repair_vcf()`

```r
repair_vcf(
  vcf,
  output_vcf,
  steps = c("drop_invalid_records", "filter_biallelic", "fix_marker_ids", "remove_duplicates",      "sort"),
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
  overwrite = FALSE
)
```

## `run_beagle()`

```r
run_beagle(
  gt,
  out,
  ref = NULL,
  map = NULL,
  mode = c("impute", "phase"),
  nthreads = 1L,
  memory_gb = 4L,
  seed = NULL,
  beagle_jar = NULL,
  java = NULL,
  extra_args = character(),
  run_dir = NULL,
  overwrite = FALSE
)
```

## `run_beagle_by_chrom()`

```r
run_beagle_by_chrom(
  gt,
  out_dir,
  ref = NULL,
  map = NULL,
  chrom_filter = NULL,
  parallel = c("none", "future", "mclapply"),
  workers = 1L,
  ...
)
```

## `simulate_vcf()`

```r
simulate_vcf(
  founder_vcf,
  map,
  output_dir,
  scenario = c("outbred", "f2", "selfed"),
  inheritance = c("diploid", "autotetraploid_bivalent"),
  founder_data = NULL,
  parents = NULL,
  n_target = 100L,
  n_reference = 100L,
  population_size = 500L,
  generations = NULL,
  observed_markers = NULL,
  missingness = list(rate = 0.05, sample_sd = 0, marker_sd = 0),
  error_rate = 0.001,
  seed = 2026L,
  max_allele_cells = 2e+07
)
```

## `split_vcf_by_chrom()`

```r
split_vcf_by_chrom(
  vcf,
  out_dir,
  chrom_filter = NULL,
  bcftools = NULL,
  overwrite = FALSE
)
```

## `validate_genotype_table()`

```r
validate_genotype_table(
  data,
  sample_col = NULL,
  encoding = c("auto", "AHB", "012", "centered", "AABB", "VCF", "REFHET", "IUPAC"),
  ploidy = 2L,
  missing_codes = c(NA, ".", "--", "", "N", "NN", "./.")
)
```

## `validate_marker_map()`

```r
validate_marker_map(
  map,
  genotype_table = NULL,
  marker_ids = NULL,
  require_alleles = FALSE,
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is")
)
```

## `validate_reference_match()`

```r
validate_reference_match(
  target_vcf,
  reference_vcf,
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
  require_phased = FALSE
)
```

## `validate_vcf()`

```r
validate_vcf(
  vcf,
  target = c("general", "beagle", "plink"),
  max_records_to_scan = Inf,
  progress = interactive()
)
```

## `vcf_to_dosage()`

```r
vcf_to_dosage(
  vcf,
  samples = NULL,
  markers = NULL
)
```

## `vcf_to_hardcalls()`

```r
vcf_to_hardcalls(
  vcf,
  samples = NULL,
  markers = NULL,
  encoding = c("012", "AHB", "AABB")
)
```

## `vcf_to_plink()`

```r
vcf_to_plink(
  vcf,
  plink_prefix,
  plink = NULL,
  overwrite = FALSE
)
```
