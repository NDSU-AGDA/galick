# Run Beagle

Phase a tiny VCF, inspect the `beagler_run`, and read the JSON run
record.

!!! note "Beagle required"

    This example assumes Java + Beagle 5.x are reachable. See
    [Installation](../getting-started/installation.md#beagle) — the
    fastest path is `install_beagle()`.

## Set up

```r
library(beagler)

# If you don't have Beagle yet:
# install_beagle()
beagle_version()
#> [1] "5.4"

# Build a tiny gt VCF from the bundled CSV fixture
csv <- system.file("extdata", "tiny.csv", package = "beagler")
map <- system.file("extdata", "tiny_map.tsv", package = "beagler")
gt <- tempfile(fileext = ".vcf")
csv_to_vcf(csv, map, gt)
```

## Phase

```r
out_dir <- tempfile("beagler_demo_")
dir.create(out_dir)
out_prefix <- file.path(out_dir, "phased")

run <- run_beagle(gt, out_prefix, mode = "phase", nthreads = 2L)
run$status
#> [1] 0
run$runtime_sec
#> [1] 1.42
run$output_paths
#> [1] ".../phased.vcf.gz"
```

## What's in `beagler_run`

```r
run
```

```text
── beagler_run: OK ─────────────────────────────────────────────────────
Runtime: 1.42 sec
Command: /usr/bin/java -Xmx4g -jar /opt/beagle/beagle.jar gt=/tmp/.../gt.vcf out=/tmp/.../phased nthreads=2 impute=false

── Outputs
• /tmp/beagler_demo_xxxxxx/phased.vcf.gz

── Logs
• stdout: /tmp/beagler_demo_xxxxxx/stdout.log
• stderr: /tmp/beagler_demo_xxxxxx/stderr.log
```

## Run record

```r
record <- jsonlite::fromJSON(run$run_record)
record$schema_version
#> [1] "1.0.0"

record$versions
#> $beagler  "0.1.0"
#> $r        "4.5.3"
#> $os       "darwin22.6.0"
#> $beagle   "5.4"
#> $java     "17.0.5"

record$inputs
#>     path  sha256
#> 1   ...   4f29a18b3c127cc8...

record$outputs
#>     path  sha256
#> 1   ...   7d2e9a1b48ff21de...
```

## Generate a Markdown report

```r
imputation_report(run, "report.md")
cat(readLines("report.md")[1:5], sep = "\n")
#> # beagler run report
#>
#> **Status:** OK (0)  |  **Runtime:** 1.42 sec
```

## Reference panels

To impute against a reference panel, pass `ref = `:

```r
run_beagle(
  gt = "panel.vcf.gz",
  out = "imputed",
  ref = "ref.vcf.gz",          # must have ref.vcf.gz.tbi sibling
  map = "genetic_map.txt",     # PLINK genetic map (cM)
  mode = "impute",
  nthreads = 8L,
  memory_gb = 16L
)
```

`ref` **must** be `.vcf.gz` (BGZF) with a tabix index. Plain `.vcf`
references are a hard error — Beagle requires random access. See
[Run Beagle](../reference/run-beagle.md#run_beagle) for the full list
of options.

Use `validate_reference_match()` first if you're unsure whether your
target and reference align:

```r
validate_reference_match("panel.vcf.gz", "ref.vcf.gz")
```

## Failure modes

Pre-flight problems with the input raise immediately and point you
at `repair_vcf()`:

```r
broken_vcf <- system.file("extdata", "tiny_dirty.vcf", package = "beagler")
tryCatch(
  run_beagle(broken_vcf, "out"),
  error = function(e) class(e)[1]
)
#> [1] "beagler_error_run_beagle_invalid_input"
```

Beagle's own complaints come back as a `beagler_run` with
`status != 0`:

```r
run <- run_beagle(gt, "out", extra_args = c("nonsense_flag=true"))
run$status
#> [1] 13
print(run)         # shows the tail of stderr
```
