# Roadmap

The full milestone log lives in
[`beagler_roadmap.md`](https://github.com/DidierMurilloF/beagleR/blob/main/beagler_roadmap.md)
in the repository.

## v0.1.0 (this release)

- [x] **M0** — Pivot from a Rust `extendr` backend to pure R.
- [x] **M1** — Foundation utilities (S3 result classes, process
  runner, path utils, chrom utils, run-record JSON schema).
- [x] **M2** — `validate_genotype_table()` with encoding detection.
- [x] **M3** — `validate_marker_map()` with chrom-style normalization.
- [x] **M4** — `csv_to_vcf()` and `hapmap_to_vcf()` with strict-fail
  validation.
- [x] **M5** — `validate_vcf()` with the tolerant `parse_vcf_raw`
  parser and target-aware checks.
- [x] **M6** — `repair_vcf()` with post-repair re-validation.
- [x] **M7** — `run_beagle()`, `install_beagle()`, `beagle_version()`.
- [x] **M10** — Vignettes, fixtures, NEWS, MkDocs site.
- [x] **M11** — BGZF / `.vcf.gz` validation path,
  `imputation_report()`, `validate_reference_match()`.
- [x] **M12 + M13** — `vcf_to_dosage()`, `vcf_to_hardcalls()`,
  `vcf_to_plink()`, `plink_to_vcf()`.
- [x] **M14** — Chromosome workflows: `split_vcf_by_chrom()`,
  `run_beagle_by_chrom()`, `merge_beagle_outputs()`.
- [x] **M15** — Masked-genotype benchmarking: `mask_genotypes()`,
  `compare_imputed_vs_true()`, `benchmark_imputation()`.

## Permanently out of scope

These are downstream analytical tasks. Other R packages serve them
well (`AGHmatrix`, `rrBLUP`, `BGLR`, `sommer`, `GAPIT`, `BWGS`).
beagler will not add them:

- GRM / kinship-matrix construction
- GWAS, PCA, mixed models
- Genomic prediction, breeding-value estimation, selection indices
- Heritability estimation
- Polyploid imputation (Beagle is diploid-only)
- Beagle 4.x support
- A generic `run_plink()` wrapper

## Possibly later

- Excel input — currently deferred. Save as CSV first.
- A faster `bcftools concat` path inside `merge_beagle_outputs()`.
- Streaming / sampled VCF validation for files with millions of
  records.
- Marker-map merging across builds.
- Templated Rmd / Quarto reports.
