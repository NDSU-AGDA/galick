---
title: Get started with Galick genotype workflows in R
description: Convert bundled CSV genotypes to VCF, validate inputs, run Beagle and extract allele counts in a concise Galick quick-start workflow.
---
# Quick tour

Run this example after installing Galick. The bundled five-sample, 20-marker data
demonstrate conversion and validation; they are not an accuracy benchmark.

## Convert and inspect

```r
library(galick)

csv <- system.file("extdata", "tiny.csv", package = "galick")
map <- system.file("extdata", "tiny_map.tsv", package = "galick")
vcf <- tempfile(fileext = ".vcf")

conversion <- csv_to_vcf(csv, map, vcf, encoding = "AHB")
stopifnot(conversion$ok)

validation <- validate_vcf(vcf, target = "beagle")
stopifnot(validation$ok)
validation$warnings

counts <- vcf_to_dosage(vcf)
hardcalls <- vcf_to_hardcalls(vcf, encoding = "AHB")
stopifnot(identical(dim(counts), c(5L, 20L)))
counts[1:3, 1:4]
```

Rows are samples and columns are markers. `counts` comes from GT allele counts,
not expected DS dosage. Warnings do not make `$ok` false; inspect them separately.

## Repair a separate dirty fixture

```r
dirty <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
clean <- tempfile(fileext = ".vcf")

before <- validate_vcf(dirty, target = "beagle")
before$errors

repair <- repair_vcf(dirty, clean)
stopifnot(repair$ok)
stopifnot(validate_vcf(clean, target = "beagle")$ok)
repair$summary$step_summary
```

Repair can remove records and writes a minimal GT-only VCF. It is not a general
lossless normalization step. Inspect the step summary and warnings before using
the result; failed post-repair validation prevents publication of the output.

## Run Beagle when configured

This section requires Java and a Beagle JAR. It uses the converted VCF directly;
the preceding dirty-fixture example is independent.

```r
run <- run_beagle(
  gt = vcf,
  out = tempfile("tour_phased_"),
  mode = "phase",
  nthreads = 1L,
  seed = 42L
)
stopifnot(identical(run$status, 0L))

calls <- vcf_to_hardcalls(run$output_paths[[1]], encoding = "012")
record <- jsonlite::fromJSON(run$run_record)
record$versions
```

Keep the JSON run record with the output. It records the command, versions,
settings, input/output hashes and indexing outcome. `summary(run)` includes the
stderr tail for a failed run.

## Next steps

- [Result classes](result-classes.md): errors, warnings, statuses and decisions.
- [Linked simulation](../examples/simulate-vcf.md): founder-conditioned data.
- [Held-out evaluation](../examples/evaluate-imputation.md): tune and test.
- [Panel design](../examples/design-marker-panel.md): select measured SNPs.
- [Reference](../reference/index.md): all 28 public functions.
