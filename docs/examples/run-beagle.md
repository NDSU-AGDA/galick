# Run Beagle

This example requires Java 8+ and a configured Beagle 5 JAR. It demonstrates
phasing of a tiny fixture, not a scientifically representative imputation study.
See [Installation](../getting-started/installation.md#beagle).

## Prepare the input

```r
library(galick)

csv <- system.file("extdata", "tiny.csv", package = "galick")
map <- system.file("extdata", "tiny_map.tsv", package = "galick")
vcf <- tempfile(fileext = ".vcf")
conversion <- csv_to_vcf(csv, map, vcf, encoding = "AHB")
stopifnot(conversion$ok)
```

## Phase and inspect

```r
run <- run_beagle(
  gt = vcf,
  out = tempfile("phased_"),
  mode = "phase",
  nthreads = 1L,
  seed = 42L
)
stopifnot(identical(run$status, 0L))

summary(run)
calls <- vcf_to_hardcalls(run$output_paths[[1]], encoding = "012")
record <- jsonlite::fromJSON(run$run_record)
record$versions
record$args$index_action
```

Use the returned path instead of constructing an assumed output filename.
Failed runs have no successful output path. `summary(run)` includes the stderr
tail when available. Runtime and accuracy are not fixed properties of this example.

## Reference-based imputation

For your own data, pass a reference and a genetic map explicitly:

```r
compatibility <- validate_reference_match(
  target_vcf = "target.vcf.gz",
  reference_vcf = "reference.vcf.gz",
  require_phased = TRUE
)
stopifnot(compatibility$ok)

run <- run_beagle(
  gt = "target.vcf.gz",
  out = "imputed",
  ref = "reference.vcf.gz",
  map = "genetic.map",
  nthreads = 4L,
  memory_gb = 8L,
  seed = 42L
)
stopifnot(identical(run$status, 0L))
```

These paths are placeholders. Target and reference must share compatible variant
identities; reference calls must be complete, phased and diploid, with disjoint
sample IDs. Plain VCF references are accepted; a tabix index is not a runner
precondition. No FASTA-based allele normalization or relationship inference is
performed. Supply a representative reference and verify assembly identity.

## Keep the run record

```r
report <- tempfile(fileext = ".md")
imputation_report(run, report)
```

Keep the full JSON record and logs with the imputed VCF. They document the job,
not its empirical accuracy. Follow [Held-out evaluation](evaluate-imputation.md)
to compare configurations without selecting on test outcomes.
