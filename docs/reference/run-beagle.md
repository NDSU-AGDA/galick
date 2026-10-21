# Run Beagle

Three exports: the runner plus two convenience helpers.

## `run_beagle()`

```r
run_beagle(
  gt,
  out,
  ref         = NULL,
  map         = NULL,
  mode        = c("impute", "phase"),
  nthreads    = 1L,
  memory_gb   = 4L,
  seed        = NULL,
  beagle_jar  = NULL,
  java        = NULL,
  extra_args  = character(),
  run_dir     = NULL,
  overwrite   = FALSE
) -> beagler_run
```

Wraps Beagle 5.x via `processx`. Returns a `beagler_run` even when
Beagle exits non-zero — call `print(run)` to see the captured stderr.

### Inputs

- `gt` — plain `.vcf` or `.vcf.gz` target genotype file.
- `out` — output **prefix** (no `.vcf` / `.vcf.gz` suffix).
  Beagle writes `<out>.vcf.gz`.
- `ref` *(optional)* — `.vcf` or `.vcf.gz` reference panel; no index is
  required. Calls must be complete, phased and diploid. Sample IDs must
  not overlap the target, and shared variant alleles must match.
- `map` *(optional)* — PLINK-format genetic map (cM positions).

### Pre-flight

`run_beagle()` runs `validate_vcf(gt, target = "beagle")` before
launching Beagle. If validation finds errors, the function aborts
with `beagler_error_run_beagle_invalid_input` and points you at
`repair_vcf()`.

### Resource controls

```r
run_beagle(
  gt = "panel.vcf.gz",
  out = "imputed",
  memory_gb = 16L,        # JVM -Xmx16g
  nthreads = 8L,
  seed = 42L              # for reproducibility
)
```

### `extra_args`

Unique `key=value` strings for additional Beagle options. Typed arguments
(`gt`, `out`, `ref`, `map`, `nthreads`, `seed`, `impute`) cannot be overridden:

```r
run_beagle(gt, out, extra_args = c("burnin=10", "iterations=20"))
```

Beagle validates option names and values beyond those typed arguments.

### What gets captured

Inside `run_dir` (auto-created if `NULL`):

- `stdout.log`, `stderr.log` — captured streams.
- `run_record.json` — v1.0.0 audit JSON with command, versions,
  args, runtime, status, and SHA-256 hashes of inputs and outputs.

The `beagler_run$output_paths` slot lists the validated `<out>.vcf.gz`.
Failed runs have no successful output path. On a failed overwrite the previous
VCF and indexes are restored; a stale file cannot satisfy post-flight checks.

### Indexing

After a successful run, the function looks for `tabix -p vcf` first,
then `bcftools index -t`. If neither tool is on `PATH`, the run
record marks indexing as `skipped` (not an error).

## `install_beagle()`

```r
install_beagle(version = "5.4", dest = NULL) -> path
```

Downloads a pinned tested Beagle JAR to
`tools::R_user_dir("beagler", "cache")` (or `dest`).

```r
install_beagle()
beagle_version()
#> [1] "5.4"
```

`version` must equal `"5.4"` in v0.1.0; other values raise
`beagler_error_beagle_unsupported_version`. The pinned URL is
`beagle.29Oct24.c8e.jar` — see [Installation](../getting-started/installation.md#beagle)
for caveats.

## `beagle_version()`

```r
beagle_version(beagle_jar = NULL) -> character(1)
```

Resolves the JAR the same way `run_beagle()` does and parses its
version banner. Useful as a sanity check after `install_beagle()` or
when verifying `BEAGLE_PATH`:

```r
beagle_version()                            # auto-discovered
beagle_version("/opt/beagle/beagle.jar")    # explicit path
```
