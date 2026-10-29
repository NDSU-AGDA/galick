---
title: Inspect Galick validation reports and experiment results
description: Read errors, warnings, VCF diagnostics, Beagle run logs and simulation, evaluation or marker-panel results returned by Galick functions.
---
# Result classes

Galick uses S3 objects and ordinary lists, matrices and paths. Not every function
returns a result class: matrix extraction returns a matrix, for example.

## `galick_result`

Validation, conversion, repair, export and comparison return `galick_result`.

| Field | Meaning |
|---|---|
| `ok` | No collected errors; warnings alone do not make this false |
| `errors`, `warnings` | Tibbles with `code`, `message` and structured `where` |
| `summary` | Function-specific counts, paths and diagnostics |
| `metadata` | Captured call, package version and time |

```r
library(galick)

path <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
result <- validate_vcf(path, target = "beagle")
result$ok
result$errors
as.data.frame(result)
```

For ordinary results, `print()` shows up to three issues of each severity and
`summary()` also shows remaining issues. VCF validation has a dedicated report
for issue overlap, affected records, unsupported-allele dependencies and optional
quality-field comparisons. `as.data.frame()` flattens scalar context when possible and
otherwise retains JSON context; it is not a genotype matrix.

## `galick_run`

Beagle and the PLINK bridge return `galick_run`, with process `status`,
`runtime_sec`, `command`, `version`, `output_paths` and log/record paths.
Beagle runs retain JSON records; the current PLINK bridge does not.

After launching Beagle, inspect `run$status` before using its outputs.
`print(run)` shows command and paths. **`summary(run)`**, not `print(run)`,
adds the stderr tail after a failure. Preflight input or tool-discovery errors
raise before a run object is returned.

Beagle records use schema version **2.0.0**, with required `versions$galick`,
`versions$r` and `versions$os`. Input files have SHA-256 hashes. The package
does not certify regulatory compliance or independently verify the biological
quality of inputs.

## Experimental subclasses

| Class | Function | Inspect |
|---|---|---|
| `galick_simulation` | `simulate_vcf()` | `paths`, `samples`, `observation`, `diagnostics`, `model` |
| `galick_evaluation` | `evaluate_imputation()` | `selection`, `metrics`, `by_group`, `runs` |
| `galick_panel_design` | `design_marker_panel()` | `selection`, `recommended_panel`, `metrics`, `search_history` |

These inherit from `galick_result`. For evaluation and panel design, `$ok` reports
execution/scoring errors, not improvement. A completed experiment may have
`no_clear_winner` or `not_confirmed`; only `confirmed` supplies a recommendation.

## Raised conditions and collected issues

Raised conditions use `galick_error_<domain>_<reason>`. Collected issue codes do
not have that prefix. Handle the specific raised class you intend to recover:

```r
tryCatch(
  validate_genotype_table(NULL),
  galick_error_invalid_input = function(error) {
    message(conditionMessage(error))
    NULL
  }
)
```
