# PLINK bridge

Two thin wrappers around PLINK 1.9 / 2.0. Both shell out via
`processx::run()` and return `galick_run` objects.

PLINK must be on `PATH` or passed explicitly via the `plink = ` argument.

## `vcf_to_plink()`

```r
vcf_to_plink(
  vcf,
  plink_prefix,
  plink     = NULL,
  overwrite = FALSE
)
```

Runs `plink --vcf <vcf> --make-bed --out <plink_prefix>`. PLINK writes
`<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam`, plus a `.log` file.

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

Runs `plink --bfile <plink_prefix> --recode vcf --out <derived_prefix>`
where `derived_prefix` is `vcf` with the `.vcf` extension stripped.

```r
result <- plink_to_vcf("data", "data.vcf")
result$status
file.exists("data.vcf")
```

## What's not in v0.1.0

A generic `run_plink()` is **explicitly out of scope**. PLINK has
hundreds of flags; wrapping them all would be a separate package.
For PLINK QC steps (`--maf`, `--geno`, `--mind`, `--hwe`, etc.),
shell out via `processx::run()` directly:

```r
processx::run(
  "plink",
  args = c("--bfile", "data",
           "--maf", "0.05",
           "--geno", "0.1",
           "--make-bed",
           "--out", "data_qc")
)
```
