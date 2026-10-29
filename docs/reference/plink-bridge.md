---
title: Convert between VCF and genetics PLINK files in R
description: Use Galick to convert VCF and PLINK BED, BIM and FAM files with executable checks, staged outputs and explicit conversion status.
---
# PLINK bridge

These functions require a genetics PLINK executable supporting the PLINK 1.9
commands shown below. They return `galick_run` objects with status, command
and output paths. Do not assume every PLINK 2 build supports these legacy commands.

Galick checks `plink1.9`, then `plink` on `PATH`, or accepts an explicit
`plink =` path. It verifies the genetics version banner and rejects PuTTY's
unrelated `plink.exe` SSH client.

## `vcf_to_plink()`

```r
vcf_to_plink(
  vcf,
  plink_prefix,
  plink     = NULL,
  overwrite = FALSE
)
```

Runs the equivalent of `plink --vcf <vcf> --make-bed --out <plink_prefix>`.
Successful conversion publishes `<prefix>.bed`, `<prefix>.bim` and `<prefix>.fam`.
Files are staged before publication; PLINK's temporary `.log` file and a separate
JSON run record are not retained by this bridge.

```r
result <- vcf_to_plink("data.vcf.gz", "data")
result$status
result$output_paths
```

If any of `<prefix>.{bed,bim,fam}` already exist and `overwrite =
FALSE`, the function aborts with `galick_error_vcf_to_plink_output_exists`.

## `plink_to_vcf()`

```r
plink_to_vcf(
  plink_prefix,
  vcf,
  plink     = NULL,
  overwrite = FALSE
)
```

Runs the equivalent of `plink --bfile <plink_prefix> --recode vcf`, then publishes
the validated VCF at the requested `vcf` path. The actual subprocess command
uses a temporary output prefix rather than writing directly over the destination.

```r
result <- plink_to_vcf("data", "data.vcf")
result$status
file.exists("data.vcf")
```

## Scope

These functions convert files. They do not perform PLINK quality-control
analyses or expose arbitrary PLINK flags; use genetics PLINK directly for those tasks.
