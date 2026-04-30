# Reports

## `imputation_report()`

```r
imputation_report(run, file) -> path (invisible)
```

Writes a one-page Markdown summary of a `beagler_run` to `file`.
Captures status, runtime, the full command, input/output paths with
truncated SHA-256 hashes, the version block, and log file paths.

When `run$run_record` exists, the function pulls input/output hashes
straight from the v1.0.0 JSON; otherwise it falls back to whatever
the in-memory `beagler_run` provides.

```r
run <- run_beagle("data.vcf", "imputed", mode = "phase")
imputation_report(run, "report.md")

cat(readLines("report.md"), sep = "\n")
```

```markdown
# beagler run report

**Status:** OK (0)  |  **Runtime:** 12.30 sec

_Generated: 2026-04-28T15:34:12-0500_

## Command

```
/usr/bin/java -Xmx4g -jar /opt/beagle/beagle.jar gt=data.vcf out=imputed nthreads=1 impute=false
```

## Inputs

- `data.vcf` (sha256: `4f29a18b3c12`)

## Outputs

- `imputed.vcf.gz` (sha256: `7d2e9a1b48ff`)

## Versions

- **beagler**: 0.1.0
- **r**: 4.5.3
- **os**: darwin22.6.0
- **beagle**: 5.4
- **java**: 17.0.5

## Logs

- stdout: `/tmp/.../stdout.log`
- stderr: `/tmp/.../stderr.log`
- run record: `/tmp/.../run_record.json`
```

The output is plain Markdown; render it to HTML/PDF with whatever
tooling you already use (pandoc, Quarto, GitHub).
