# beagler 0.1.0

First public release. 27 user-facing exports across validation,
conversion, repair, Beagle execution, dosage / hardcall extraction,
PLINK bridge, chromosome workflows, and masked-genotype benchmarking.

## Licensing and authors

- The package license is now MIT. External tools and dependencies retain their
  own licenses.
- Muyideen Yusuf is included in the package authors.

## Founder-conditioned simulation

- `simulate_vcf()` adds dependency-free Haldane meiosis and forward breeding
  simulation from phased founders and a supplied cM map. Diploid scenarios
  include within-population random mating, biparental F2 and selfed lineages.
- Random-bivalent autotetraploid inheritance is limited to outbred populations;
  quadrivalents, double reduction, preferential pairing and allopolyploids are
  unsupported. Polyploid Beagle/evaluation restrictions remain unchanged.
- Bundles include complete truth/reference VCFs, unphased observed targets,
  pedigrees, replay-checked recombination events, calibrated heterogeneous
  missingness, controlled calling errors, LD diagnostics and hashed provenance.
- Independent AlphaSimR 2.1.0 comparisons are archived with a reproducible
  development script. AlphaSimR is not a runtime dependency.

## Configuration evaluation

- `evaluate_imputation()` compares two to four Beagle configurations using
  paired sporadic masks or physical array thinning, separate tuning/test sample
  partitions, a fixed run budget and shared resource settings. It selects by
  tuning RMSE with complete predictions and group safeguards, then confirms
  against a baseline without test-set reranking.
- `beagler_evaluation` retains group and per-marker metrics, failure diagnostics,
  held-out manifests, frozen selection and versioned, hashed experiment records.
  A runnable vignette covers paired comparison and failed confirmation.

## Release validation fixes

- Benchmark alignment now requires compatible coordinates and REF/ALT,
  distinguishes GT from DS, retains missing predictions in coverage counts,
  and reports per-marker correlations with undefined cases explicitly.
- Chromosome resume verifies a complete input/settings/runtime fingerprint,
  output checksum and VCF validity. Legacy cache records are rerun.
- Chromosome concatenation rejects failed or incomplete jobs, aligns sample
  columns, reconciles headers, retains DS/GP and writes BGZF when compressed.
- GT allele indices are checked against ALT. Reference preflight checks
  duplicate coordinates, declared assemblies, phasing and sample overlap.
  Plain reference VCFs are supported without a mandatory tabix index.
- SHA-256 uses portable `digest`; the minimum supported R version is 4.1.

## Genotype export and VCF merge

- `export_genotypes()` adds VCF, PLINK PED/MAP, HapMap, CSV/TSV,
  and numeric export. VCF retains FORMAT fields, phase sets and reference
  dictionaries; compressed output requires bgzip. HapMap uses diploid IUPAC
  SNP calls. PED and HapMap report phase loss and fully mask partial calls.
  Polyploid calls remain intact in VCF and GT counts; arbitrary truncation
  is refused. `field = "DS"` explicitly selects expected dosage.
- `merge_vcf()` adds explicit `by_sample` and `by_marker` callset merges.
  It aligns compatible variant identities before duplicate occurrences,
  reindexes GT, recomputes sample-merge AC/AN, and reports omitted annotations.
  Reference dictionaries must agree and independent phase blocks stay distinct.
  Valid overlapping SNPs and indels survive marker merges. Sample collisions
  support error, suffix, or unphased coalesce policies.
- `detect_subgenome_labels()` reports chromosome naming labels without
  inferring species ploidy. `detect_genotype_ploidy()` reports the maximum
  observed GT allele-slot count.

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
milestone log and `docs/about/conventions.md` for project conventions.
