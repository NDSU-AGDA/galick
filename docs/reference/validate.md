# Validate

Four validators. Each returns a `galick_result` with `$errors`,
`$warnings`, and a function-specific `$summary`.

## `validate_genotype_table()`

Inspect a wide genotype table before conversion.

```r
validate_genotype_table(
  data,
  sample_col      = NULL,
  encoding        = c("auto", "AHB", "012", "centered",
                       "AABB", "VCF", "REFHET", "IUPAC"),
  ploidy          = 2L,
  missing_codes   = c(NA, ".", "--", "", "N", "NN", "./.")
)
```

**What it checks**

- Shape: data frame, ≥ 1 row, ≥ 2 columns.
- Sample column: auto-detected as the first character/factor column
  with no NA, 100 % unique, and tokens that don't look like a genotype
  encoding.
- Encoding: detects from data or validates against `encoding =`.
  Mixed/unknown returns errors; small alphabet wins ties (`{A}` →
  `AHB`, not `IUPAC`).
- Duplicates: sample IDs and marker names — both are **errors** (they
  break downstream conversion).
- Monomorphic markers and ≥ 50 % missing — **warnings**.

**Example**

```r
data <- data.frame(
  sample_id = c("plant_1", "plant_2", "plant_3"),
  rs1 = c("A", "H", "B"),
  rs2 = c("H", "B", "A")
)

result <- validate_genotype_table(data)
result$ok
result$summary$inferred_encoding   # "AHB"
result$summary$missing_rate_overall
```

## `validate_marker_map()`

```r
validate_marker_map(
  map,
  genotype_table   = NULL,
  marker_ids       = NULL,
  require_alleles  = FALSE,
  chrom_style      = c("auto", "numeric", "chr_prefixed", "as_is")
)
```

**What it checks**

- Required columns: `marker_id`, `chrom`, `pos`. Optional: `ref`,
  `alt`, `cm`.
- `pos` is finite, positive, integer-valued.
- `chrom` is non-empty.
- `marker_id` is non-empty (NA / "" / whitespace = error).
- Duplicates on `marker_id` and `(chrom, pos)` — both errors.
- Sortedness within chromosome — warning.
- Optional `genotype_table` (M2 wide format) or `marker_ids` (bare
  list) — cross-checks coverage.

`chrom_style != "auto"` normalizes labels first, so duplicate detection
catches `chr1` vs `1` collisions.

## `validate_vcf()`

```r
validate_vcf(
  vcf,
  target              = c("general", "beagle", "plink"),
  max_records_to_scan = Inf
)
```

Accepts plain `.vcf` or `.vcf.gz` (read transparently via `gzfile()`).

**Targets**

- `"general"` — VCF v4.2 essentials. Multiallelic / no-samples are warnings.
- `"beagle"` — adds: diploid throughout, REF ≠ ALT, multiallelic is an
  error, `vcf_beagle_not_bgzf` warning when input is plain `.vcf`.
- `"plink"` — adds: `vcf_plink_unsupported_chrom`,
  `vcf_plink_missing_ids` (PLINK requires unique non-empty IDs).

**Performance**

Full validation walks every record. For very large files, set
`max_records_to_scan` to bound the scan. Sortedness and duplicate
`(chrom, pos)` checks are skipped in partial-scan mode (they would
false-positive at chunk boundaries) and `summary$scanned_partial`
is `TRUE`.

## `validate_reference_match()`

```r
validate_reference_match(
  target_vcf,
  reference_vcf,
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
  require_phased = FALSE
)
```

Compare a target VCF to a reference panel for compatibility before
calling `run_beagle(..., ref = reference_vcf)`.

**What it checks**

- Both files parseable, with `#CHROM` headers.
- Chromosome label style: warns under `"auto"` if styles differ;
  normalizes under `"numeric"` / `"chr_prefixed"`.
- Shared chromosomes: error if zero overlap.
- REF/ALT consistency at shared `(chrom, pos)`: errors on mismatch.

`$summary` reports counts of shared / target-only / reference-only
chromosomes and markers, plus the inferred chromosome styles on each
side.

Use `require_phased = TRUE` for the runner's stricter preflight: complete phased
diploid reference calls and disjoint sample IDs. Assembly metadata is checked
for conflicts where present, but no reference FASTA or biological relationships
are checked.
