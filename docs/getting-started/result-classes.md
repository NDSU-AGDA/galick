# The result classes

beagler is S3-only. Two classes carry every function's output:

- `beagler_result` for validators, repair, conversion, and comparison.
- `beagler_run` for execution (Beagle, PLINK).

Both classes have `print`, `summary`, `format`, and `as.data.frame`
S3 methods.

## `beagler_result`

```r
structure(
  list(
    ok        = logical(1),                  # TRUE iff zero errors
    errors    = tibble(code, message, where),
    warnings  = tibble(code, message, where),
    summary   = list(),                      # function-specific
    metadata  = list(call, package_version, time)
  ),
  class = c("beagler_result", "list")
)
```

`ok` is `TRUE` iff `nrow(errors) == 0`. Warnings do **not** flip it
to `FALSE`.

`errors` and `warnings` are tibbles with the same three columns:

| column | type | meaning |
|---|---|---|
| `code` | character | stable identifier (e.g. `vcf_unsorted_records`) |
| `message` | character | human-readable, generated via `cli` |
| `where` | list-column | structured context per issue |

```r
r <- validate_vcf("data.vcf", target = "beagle")
r$ok                              # logical(1)
r$errors                          # tibble of errors
r$summary                         # function-specific named list
print(r)                          # short on-screen summary
as.data.frame(r)                  # flattens errors + warnings
```

## `beagler_run`

```r
structure(
  list(
    status        = integer(1),              # process exit code
    runtime_sec   = numeric(1),
    command       = character(),             # full argv vector
    version       = list(beagle, java, ...),
    output_paths  = character(),
    stdout_log    = character(1),
    stderr_log    = character(1),
    run_record    = character(1),            # path to JSON
    metadata      = list(call, package_version, time)
  ),
  class = c("beagler_run", "list")
)
```

`status == 0` means the external tool exited cleanly. Beagle's own
complaints (after Java launches) come back as a `beagler_run` with
`status != 0` — you can inspect stderr via `print()`:

```r
run <- run_beagle("data.vcf", "out", extra_args = c("nonsense_flag=true"))
run$status                 # !=0
print(run)                 # shows command + tail of stderr
```

The run record is a v1.0.0 JSON document with required
`schema_version`, `command`, `versions`, `runtime_sec`, `status`,
`started_at`, `ended_at`, plus optional `args`, `inputs`, and
`outputs` arrays. Inputs and outputs include SHA-256 hashes — share
them alongside the imputed VCF for reproducibility audits.

```r
jsonlite::fromJSON(run$run_record)$versions
#> $beagler  "0.1.0"
#> $r        "4.5.3"
#> $os       "darwin..."
#> $beagle   "5.4"
#> $java     "17.0.5"
```

## Error and warning conventions

User-facing errors raise via `cli::cli_abort()` with a `class`
attribute, so callers and tests can dispatch on them:

```r
tryCatch(
  validate_genotype_table(NULL),
  beagler_error_invalid_input = function(e) {
    message("Bad input. Code: ", class(e)[1])
  }
)
```

The full error-class taxonomy follows two patterns:

- `beagler_error_<domain>_<reason>` for raised errors (programmer
  problems: `NULL` inputs, wrong types, missing files).
- `<domain>_<reason>` for collected `$errors` codes (data-quality
  problems: encoding mismatch, duplicate IDs, multiallelic records).

Issue codes returned in `$errors` and `$warnings` are stable across
versions; if we change semantics we add a new code, not redefine an
existing one.
