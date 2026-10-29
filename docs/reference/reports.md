---
title: Generate Beagle imputation run reports with Galick
description: Write concise Markdown reports of Beagle status, commands, versions, runtime, output files and hashed run-record provenance.
---
# Reports

`imputation_report(run, file)` writes a plain Markdown summary and invisibly
returns the output path. It accepts a `galick_run`, not an evaluation or panel
design result. Render the Markdown separately using your existing tools.

The report includes process status, runtime, command, output paths, versions and
log locations. Input paths and abbreviated hashes are read from the run record
when present. The complete output hashes remain in the JSON run record; the
Markdown report currently lists output paths without hashes.

After a Beagle run:

```r
stopifnot(identical(run$status, 0L))
report <- tempfile(fileext = ".md")
imputation_report(run, report)
cat(readLines(report), sep = "\n")
```

Neither the report nor the JSON record establishes the biological accuracy of
imputation. Use [held-out evaluation](benchmarking.md) for measured recovery.
