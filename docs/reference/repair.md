---
title: Repair VCF files for Beagle in R
description: Sort VCF records, fill missing marker IDs and apply explicit duplicate or invalid-record filters, with post-repair validation and loss reporting.
---
# Repair VCF files

`repair_vcf()` applies five selectable technical steps and publishes a minimal
VCF only if Galick's Beagle-target validation passes. This is not a guarantee
of biological correctness: filtering can remove genuine variation, and the
meaning of unsupported alleles must be established before recoding them.

## `repair_vcf()`

```r
repair_vcf(
  vcf,
  output_vcf,
  steps = c("drop_invalid_records", "filter_biallelic",
            "fix_marker_ids", "remove_duplicates", "sort"),
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
  overwrite   = FALSE
)
```

The steps run in order. The user controls the set; nothing is run
secretly.

| Step | What it does |
|---|---|
| `drop_invalid_records` | Drops rows with bad column count, invalid POS, empty CHROM, invalid REF/ALT, REF == ALT, an absent GT field, malformed/non-diploid GT, or inconsistent FORMAT/sample subfields. |
| `filter_biallelic` | Drops rows whose ALT contains a comma. |
| `fix_marker_ids` | Replaces missing `ID` values (`.`, empty, whitespace) with `chrom_pos_ref_alt`. |
| `remove_duplicates` | Drops duplicate `marker_id` rows and duplicate `(chrom, pos)` rows, keeping the first occurrence after prior steps. |
| `sort` | Orders by `(chrom, pos)`. |

Missing diploid calls such as `./.` are retained; this is not complete-case
filtering. Dropping an invalid record removes that marker for all samples.

## What's **not** repairable

These show up as `$errors` and the function writes nothing:

- `repair_vcf_missing_chrom_header` — no `#CHROM` line.
- `repair_vcf_no_samples` — Beagle needs at least one sample.
- `repair_vcf_duplicated_samples` — no obvious resolution.

## Output

The output is a **minimal** VCF: only `CHROM`, `POS`, `ID`, `REF`,
`ALT`, and `FORMAT=GT`. Other FORMAT subfields (e.g. `DP`, `GQ`) are
dropped on output — `repair_vcf()` extracts the GT subfield from each
sample cell.

```r
# Input has FORMAT=GT:DP, sample cells like "0/1:30"
# Output has FORMAT=GT, sample cells like "0/1"
```

## Post-repair validation

After applying steps, the function writes the output to a tempfile,
runs `validate_vcf(target = "beagle")` on it, and only promotes the
tempfile to `output_vcf` if **zero errors** remain.

If something went wrong (e.g., the user excluded a step that would
have addressed a condition), the result reports `repair_vcf_remaining_*`
errors that wrap the failed validator codes:

```r
# Multiallelic record present; user excluded filter_biallelic
result <- repair_vcf("dirty.vcf", "out.vcf",
                     steps = "drop_invalid_records")
result$ok
result$errors$code
file.exists("out.vcf")
```

## Chromosome-style normalization

`chrom_style != "auto"` rewrites CHROM labels in the output **and**
uses the normalized form for duplicate detection and sorting. This
catches `chr1` + `1` collisions:

```r
# Input has both chr1:100 and 1:100 — collide under "numeric"
repair_vcf("messy.vcf", "out.vcf", chrom_style = "numeric")
```

## Example: full repair

```r
dirty <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
clean <- tempfile(fileext = ".vcf")
result <- repair_vcf(dirty, clean)

result$summary$step_summary
```

`$summary$step_summary` records what each step did, and step actions
also surface as `$warnings` so they show up in `print(result)`.
