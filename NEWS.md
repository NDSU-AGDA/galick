# beagler 0.1.0

First public release. 21 user-facing exports across validation,
conversion, repair, Beagle execution, dosage / hardcall extraction,
PLINK bridge, chromosome workflows, and masked-genotype benchmarking.

## Validate

- `validate_genotype_table()` — wide genotype table: encoding detection,
  missingness, monomorphic markers, duplicate sample/marker checks.
- `validate_marker_map()` — required columns, position validity, chromosome
  style, sortedness, optional cross-check against a genotype table or a
  marker-id list.
- `validate_vcf()` — `.vcf` or `.vcf.gz` under `target = "general"`,
  `"beagle"`, or `"plink"`. Tolerant raw parser; per-target severity for
  multiallelic / no-samples; `max_records_to_scan` for partial scans.
- `validate_reference_match()` — compare a target VCF and a reference
  panel for chromosome compatibility, shared marker positions, and
  REF/ALT consistency at shared sites. `chrom_style` normalization
  catches `chr1` vs `1` collisions.

## Convert (input → VCF)

- `csv_to_vcf()` — CSV to VCF v4.2; runs M2 + M3 validators plus
  M4-specific stricter checks (biallelic A/C/G/T, marker coverage).
  Supports `AHB` / `012` / `centered` / `AABB` / `VCF` / `REFHET` /
  `IUPAC` encodings. Strict-fail: no partial output ever.
- `hapmap_to_vcf()` — TASSEL HapMap (`.hmp.txt`) to VCF v4.2; per-cell
  IUPAC consistency checks; optional `build` writes `##reference=`.

## Convert (extract from VCF)

- `vcf_to_dosage()` — numeric matrix of ALT-allele counts
  (samples × markers). Reads `.vcf` or `.vcf.gz`.
- `vcf_to_hardcalls()` — character matrix of hard calls under `"012"`,
  `"AHB"`, or `"AABB"`.

## Repair

- `repair_vcf()` — five user-controllable steps (`drop_invalid_records`,
  `filter_biallelic`, `fix_marker_ids`, `remove_duplicates`, `sort`)
  that rewrite a minimal VCF. Post-repair re-validates against
  `target = "beagle"` and refuses to write anything that wouldn't pass.

## Run Beagle

- `run_beagle()` — Beagle 5.x via `processx`. Java/JAR discovery, JVM
  heap, pre-flight VCF validation, BGZF reference panel + `.tbi`
  requirement, structured argv (no shell), v1.0.0 JSON run record with
  SHA-256 hashes, best-effort `tabix` or `bcftools` indexing.
- `install_beagle()` — downloads pinned tested 5.4 build
  (`beagle.29Oct24.c8e.jar`).
- `beagle_version()` — version string from the resolved JAR.

## Chromosome workflows

- `split_vcf_by_chrom()` — pure-R splitter; one `.vcf` per CHROM.
- `run_beagle_by_chrom()` — splits, runs Beagle per chromosome, with
  resume logic (skips chromosomes whose `run_record.json` reports
  `status == 0` and matching input SHA-256). `parallel` supports
  `"none"`, `"future"` (via `future.apply`), `"mclapply"`.
- `merge_beagle_outputs()` — concatenates per-chromosome outputs into a
  single VCF (`.vcf` or `.vcf.gz`). Pure R; bcftools fast path is
  forward-compatible.

## PLINK bridge

- `vcf_to_plink()` — shells out to PLINK with `--vcf <vcf> --make-bed`.
- `plink_to_vcf()` — `--bfile <prefix> --recode vcf`.

## Reports

- `imputation_report()` — Markdown summary of a `beagler_run` (status,
  runtime, command, inputs/outputs with truncated SHA-256s, version
  block, log paths).

## Benchmarking

- `mask_genotypes()` — randomly mask a fraction of non-missing GT
  cells; reproducible with `seed`.
- `compare_imputed_vs_true()` — per-cell concordance, dosage R²,
  RMSE on shared cells. Optional `masked_vcf` filters to the cells
  imputation had to fill.
- `benchmark_imputation()` — orchestrator: mask → Beagle → compare.

## Object model

S3 only (no R6, no S4). Two classes:

- `beagler_result` for validators, repair, conversion, comparison.
  Slots: `ok`, `errors`, `warnings`, `summary`, `metadata`. Methods:
  `print`, `summary`, `format`, `as.data.frame`.
- `beagler_run` for execution. Slots: `status`, `runtime_sec`,
  `command`, `version`, `output_paths`, `stdout_log`, `stderr_log`,
  `run_record`, `metadata`. Same S3 methods.

## Notable design decisions

- Tolerant parser (`parse_vcf_raw`) backs every VCF reader; the
  validator translates parse diagnostics into `$errors` / `$warnings`
  codes rather than aborting on shape problems.
- Errors raise via `cli::cli_abort()` only for programmer/shape problems
  (NULL inputs, wrong types, unsupported file extensions). Data-quality
  issues are *collected* into `$errors` so callers can inspect.
- Run records are JSON Schema 2020-12 v1.0.0 with required
  `schema_version`; validated before write.
- `centered` (-1/0/1) is supported alongside the other six encodings.
- chromosome-style normalization (`"chr1"` ↔ `"1"`) catches duplicate
  collisions across naming styles.
- BGZF input handled transparently via `gzfile()` (sequential reads
  work without BGZF-specific tooling; tabix is only needed for random
  access).

## Pre-history (0.0.0.9000)

The package was originally a thin R wrapper over a Rust `extendr`
backend. M0 removed the Rust scaffolding entirely and pivoted to pure R
with a test-first workflow. See `beagler_roadmap.md` for the M0–M15
milestone log and `AGENTS.md` for project conventions.
