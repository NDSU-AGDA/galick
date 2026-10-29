---
title: Convert CSV and HapMap genotype data to VCF in R
description: Convert breeder genotype tables and TASSEL HapMap files to validated VCF using explicit encodings, sample IDs and reference-allele maps.
---
# Convert to VCF

Two functions that take breeder-style genotype tables and emit VCF v4.2.

Both run the relevant validators internally and **refuse partial
output**: if validation collects any errors, no VCF is written and
the result carries the aggregated codes.

## `csv_to_vcf()`

```r
csv_to_vcf(
  csv,
  marker_map,
  output_vcf,
  sample_col    = NULL,
  encoding      = "auto",
  ploidy        = 2L,
  chrom_style   = "auto",
  ahb_ref       = c("A", "B"),
  aabb_ref_homo = c("AA", "BB"),
  missing_codes = c(NA, ".", "--", "", "N", "NN", "./."),
  overwrite     = FALSE
)
```

`csv` may be a path to a CSV file or an in-memory data frame.
`marker_map` may be a path (TSV/CSV) or a data frame. Both must align
on marker IDs.

**Pipeline**

1. `validate_genotype_table()` (M2)
2. `validate_marker_map(map, genotype_table = csv, require_alleles = TRUE)` (M3)
3. M4-specific checks (additive on top of M2/M3):
    - `csv_to_vcf_unmappable_markers` — markers in CSV but not in map.
    - `csv_to_vcf_invalid_alleles` — REF/ALT outside `{A, C, G, T}`.
    - `csv_to_vcf_unsupported_extension` — `output_vcf` not ending in `.vcf`.
    - `csv_to_vcf_output_exists` — refuses to clobber unless `overwrite = TRUE`.
4. Recode each cell via the marker's REF/ALT.
5. Write VCF v4.2 with `##fileformat`, `##fileDate`, `##source`,
   `##contig=<ID=...>` per chromosome, `##FORMAT=<ID=GT,...>`.

**REF/ALT mapping conventions**

For non-DNA encodings, `ahb_ref` and `aabb_ref_homo` swap which token
codes homo-REF:

```r
csv_to_vcf(csv, map, "out.vcf", encoding = "AHB",
           ahb_ref = "B")     # data 'A' now codes homo-ALT
csv_to_vcf(csv, map, "out.vcf", encoding = "AABB",
           aabb_ref_homo = "BB")
```

For `IUPAC`, the REF/ALT pair comes from `marker_map` per marker, so
no convention argument is needed.

## `hapmap_to_vcf()`

```r
hapmap_to_vcf(
  hmp,
  output_vcf,
  build         = NA_character_,
  chrom_style   = c("auto", "numeric", "chr_prefixed", "as_is"),
  missing_codes = c(NA, ".", "--", "", "N", "NN"),
  overwrite     = FALSE
)
```

Reads TASSEL HapMap (`.hmp.txt`): standard 11 metadata columns
(`rs#`, `alleles`, `chrom`, `pos`, ...) followed by per-sample columns
of single-character IUPAC codes.

The `alleles` column ("A/T") gives REF/ALT directly; the function
decodes IUPAC heterozygous codes per marker (e.g. `W` is the het code
for `A/T`).

Alleles are used as declared. The converter does not use the HapMap `strand`
column to reverse-complement calls or check them against a reference FASTA.
Normalize strand and reference orientation upstream before biological analysis.

**Optional `build`**

```r
hapmap_to_vcf("input.hmp.txt", "out.vcf", build = "GRCh38")
```

Writes `##reference=GRCh38` to the VCF header.

**v0.1.0 limitations**

- Biallelic SNPs only — REF/ALT must be from `{A, C, G, T}` and
  distinct.
- Allele field uses `/` as the separator; comma-delimited variants
  (`A,T`) are rejected.
- Multiallelic markers (`A/T/G`) error out.
