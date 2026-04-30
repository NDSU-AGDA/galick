# beagler Development Roadmap

## Vision

**beagler** is an R toolkit for converting, validating, repairing, imputing, and exporting plant- and animal-breeding genotype data using **Beagle** and **PLINK** as external command-line tools.

The package focuses on **genotype preprocessing and imputation workflows**. It is not a downstream-analysis package.

## Audience

Plant and animal breeders working with messy real-world inputs:

- CSV genotype tables exported from arrays or pipelines
- HapMap (`.hmp.txt`) files from TASSEL / GBS workflows
- PLINK binary sets (`.bed` / `.bim` / `.fam`) from collaborators
- VCF files of varying compliance from multiple vendors

## Permanently out of scope

These are downstream analytical tasks. Other R packages serve them well (`AGHmatrix`, `rrBLUP`, `BGLR`, `sommer`, `GAPIT`, etc.). They will not be added to beagler:

- GRM / kinship-matrix construction
- GWAS, PCA, mixed models
- Genomic prediction, breeding-value estimation, selection indices
- Heritability estimation
- Population-genetics statistics
- Polyploid imputation (Beagle is diploid-only; we detect and refuse politely)
- Beagle 4.x support (we target Beagle 5.x only)
- A generic `run_plink()` wrapper (PLINK has hundreds of flags; we are not wrapping them)

---

## Strategic decisions (settled before any code is written)

### S1. Package name

The repo is `beagleR`; `DESCRIPTION` says `beagler`. **Decision:** the R package name is **`beagler`** (lowercase, idiomatic). The repository may stay `beagleR` for branding; only the package name and `library()` call must match.

### S2. Dependency policy

- **CRAN-only R dependencies for v0.1.0**: `cli`, `data.table`, `jsonlite`, `processx`, `R.utils`, `tibble`, `vroom`, `withr`. `jsonlite` lands in M1 (run-record JSON); the rest are added per milestone as they're consumed.
- For VCF / BGZF I/O, **shell out to `bcftools` / `tabix` / `bgzip`** at runtime when present; fall back to a small in-house reader for plain `.vcf` and `R.utils::gunzip` for `.vcf.gz`.
- **Do not** depend on `Rsamtools` / `VariantAnnotation` in v0.1.0 (Bioconductor pulls in a heavy tree). Re-evaluate for v0.2.0 if VCF I/O becomes the bottleneck.
- Java and the Beagle JAR are **runtime requirements**, never bundled. Detected at run time.

`SystemRequirements` in `DESCRIPTION` will read:

```
Java (>= 8); Beagle (>= 5.x); optional: bcftools, tabix, bgzip, plink (>= 1.9)
```

### S3. Distribution

- v0.1.0: GitHub-only (`pak::pak("DidierMurilloF/beagleR")` / `remotes::install_github`).
- v0.2.0: target CRAN once the dependency surface and SystemRequirements are stable.
- The Beagle JAR cannot ship in a CRAN tarball. An optional helper `install_beagle(version = "5.4")` will download the **pinned default version (Beagle 5.4)** to `tools::R_user_dir("beagler", "cache")` on first use. The version string is a literal default, **not** a "latest version" lookup; auto-detection of the newest upstream release is deferred until a stable release feed (or a maintained mirror) is identified. Users who want a different Beagle 5.x release pass `version = "5.x"` explicitly.

### S4. Return-object schemas (frozen at M1, used everywhere)

Every validation, repair, conversion, and execution function returns a typed S3 object. Schemas are defined once in M1 and never re-invented per function.

```r
structure(
  list(
    ok = logical(1),
    errors = data.frame(code, message, where),
    warnings = data.frame(code, message, where),
    summary = list(),
    metadata = list(call, package_version, time)
  ),
  class = c("beagler_result", "list")
)
```

`ok` is `TRUE` iff `nrow(errors) == 0` (warnings do not flip it). `summary` is a function-specific named list (encoding details, missingness rates, per-step counts, etc.). `metadata` records the call expression, package version, and timestamp.

```r
structure(
  list(
    status = integer(1),
    runtime_sec = numeric(1),
    command = character(),
    version = list(beagle, java),
    output_paths = character(),
    stdout_log = character(1),
    stderr_log = character(1),
    run_record = character(1),
    metadata = list(call, package_version, time)
  ),
  class = c("beagler_run", "list")
)
```

`status` is the process exit code (0 = success). `command` is the full argv as a vector — never a shell-composed string. `output_paths` are the files actually written. `stdout_log`, `stderr_log`, and `run_record` are absolute paths.

S3 methods on both classes (defined in M1): `print`, `summary`, `format`, `as.data.frame`.

### S5. Error and warning conventions

- All user-facing errors raised via `cli::cli_abort()` with a `class` attribute (e.g. `beagler_error_invalid_encoding`) so users and tests can dispatch on them.
- Validation issues are *collected* into the `$warnings` slot of `beagler_result`, **not** raised as R warnings, so they don't pollute logs in batch pipelines.
- Catastrophic problems (file not found, Java missing, JAR missing, Beagle 4.x detected) raise immediately.

---

## Milestone 0 — Pivot cleanup ("decontamination")

**Goal:** make the repo a clean, pure-R package skeleton with CI green and a real test suite scaffold.

This must finish before M1 starts. Nothing in M1+ assumes a Rust toolchain.

### 0.1 Files and directories to remove

```
src/                                # entire tree
  beagler.so
  entrypoint.c
  entrypoint.o
  Makevars
  Makevars.ucrt
  Makevars.win
  beagler-win.def
  .gitignore
  rust/                             # entire subtree
    Cargo.toml
    Cargo.lock
    src/lib.rs
    target/                         # all build artefacts
R/extendr-wrappers.R
R/convert_csv_to_vcf.R              # see 0.4 for rationale
man/run_beagle_internal.Rd
```

### 0.2 `DESCRIPTION` changes

- Remove `SystemRequirements: Cargo (Rust's package manager), rustc`
- Remove `Config/rextendr/version`
- Add `SystemRequirements: Java (>= 8); Beagle (>= 5.x); optional bcftools, tabix, plink`
- Replace `Imports:` with `data.table, vroom, processx, withr, R.utils, tibble, cli`
- Add `Suggests: testthat (>= 3.0.0), covr, knitr, rmarkdown, lintr`
- Add `Config/testthat/edition: 3`
- Set `Version: 0.0.0.9000`
- Keep `License: GPL (>= 3)` and the existing `Authors@R`.

### 0.3 `NAMESPACE` changes

- Remove `useDynLib(beagler, .registration = TRUE)`
- Regenerate from roxygen2; let it stay almost empty until M1 introduces the first real exports.

### 0.4 Existing R code

- `R/beagle.R::run_beagle()` — replace the body with a placeholder that calls `cli::cli_abort()` saying "use M7 implementation". Move the BEAGLE_PATH check into the future `find_beagle_jar()`. The actual reimplementation lands in M7.
- `R/convert_csv_to_vcf.R::convert_csv_to_vcf()` — **delete**. The function name is misleading: it requires the input CSV to already contain `CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT`, so it is a row-formatter, not a converter. Worse, its `recode_snp()` helper attempts to handle both `0/1/2` and `-1/0/1` codings simultaneously by combining their R-value matches with `||`, which is logically impossible: in the `-1/0/1` encoding, `0` should mean *heterozygous* and `1` should mean *homozygous-alt*, but the first matching branch always claims them for the `0/1/2` interpretation, so the `-1/0/1` path is silently mis-coded. It will be reimplemented correctly in M4 with an explicit `encoding` argument and full test coverage.

### 0.5 `.Rbuildignore`

```
^src/$
^.*\.Rproj$
^\.Rproj\.user$
^LICENSE\.md$
^\.github$
^beagler_roadmap\.md$
^vignettes/precomputed/
^tests/fixtures/large/
^README\.Rmd$
```

### 0.6 `.gitignore` correction

The current `.gitignore` excludes `*.csv`, `*.vcf`, `*.gz` globally — which will silently drop any test fixture you add. Replace with:

```
.Rproj.user
.Rhistory
.RData
.Ruserdata
src/*.o
src/*.so
src/*.dll
*.log
/tmp/
/scratch/
# Fixtures under inst/extdata/ and tests/testthat/fixtures/ are NOT ignored.
```

### 0.7 Test scaffold

```
tests/testthat.R
tests/testthat/_snaps/
tests/testthat/helper-skip.R
tests/testthat/test-genotype-recoding.R   # see 0.9
```

### 0.8 CI

`.github/workflows/R-CMD-check.yaml` — `r-lib/actions/setup-r-dependencies@v2`, matrix over `{ubuntu-latest, macos-latest, windows-latest} × {release, devel}`.

`.github/workflows/lint.yaml` — `lintr` on every push.

A separate `.github/workflows/integration.yaml` is added in M7 (it caches the Beagle JAR and PLINK).

### 0.9 Day-one test (written before any feature code)

`tests/testthat/test-genotype-recoding.R` is written first and asserts the *correct* mapping for both encodings:

| input value | encoding `0/1/2` | encoding `-1/0/1` |
|---|---|---|
| `-1` | (invalid) | `0/0` |
| `0` | `0/0` | `0/1` |
| `1` | `0/1` | `1/1` |
| `2` | `1/1` | (invalid) |

The test fails until M4 lands the new recoder, and then turns green. It is the canary that demonstrates the package is actually testing behaviour, not stubs.

### 0.10 Documentation hygiene

- Rewrite `README.md` to reflect the new direction: install via `pak::pak("DidierMurilloF/beagleR")`, one short example, link to vignettes (placeholder until M10).
- Create `NEWS.md` with a single "0.0.0.9000 — pivot to pure R" entry.

### Exit criteria for M0

- `R CMD check --as-cran` passes with 0 errors, 0 warnings on all three OSes (one NOTE about no exports is acceptable here).
- `git grep -i 'extendr\|rust\|cargo\|rextendr'` returns nothing.
- The day-one recoding test exists and (correctly) **fails** with a clear "not yet implemented" message.
- CI green on lint + R-CMD-check matrix.

---

## Milestone 1 — Foundation utilities

**Goal:** internal infrastructure that every later milestone reuses.

### Deliverables

- Result-object constructors and S3 methods (per the schemas in §S4 above):
  - `new_beagler_result()`, `new_beagler_run()`
  - `print`, `summary`, `format`, `as.data.frame` for both classes
- Process runner: `run_cli(cmd, args, stdout_log, stderr_log, timeout, env)` wrapping `processx::run()` with capture, timeout, and structured return.
- Path utilities: `normalize_path_safe()`, `assert_readable()`, `assert_writable()`, `with_run_dir()`.
- Logging: `log_step(run_dir, step, ...)` writes one JSON-line entry per step.
- Chromosome utilities: `normalize_chrom_label()`, `chrom_style()` (detects `chr1` vs `1` vs `Chr01`).
- Run record: `write_run_record_json(run_dir, ...)`; schema in `inst/schemas/run_record.schema.json`.
- Test helpers in `tests/testthat/helper-skip.R`:
  - `skip_if_no_java()`, `skip_if_no_beagle()`, `skip_if_no_plink()`,
  - `skip_if_no_bcftools()`, `skip_if_no_tabix()`.

### Exports

None. Everything in M1 is internal — declared via roxygen `@keywords internal` and absent from `NAMESPACE`. Visibility is enforced by `NAMESPACE`, never by a leading dot in the function name.

### Exit criteria

- ≥90% line coverage on M1 internals.
- Snapshot tests for `print` / `summary` of both result classes.
- `run_cli()` tested against `Rscript` invocations cross-platform (success, non-zero exit, timeout, stderr capture).

---

## Milestone 2 — Genotype-table validation

**Goal:** validate breeder-style wide genotype tables before any conversion.

### Function

```r
validate_genotype_table(
  data,
  sample_col = NULL,
  encoding = c("auto", "AHB", "012", "centered", "AABB", "VCF", "REFHET", "IUPAC"),
  ploidy = 2L,
  missing_codes = c(NA, ".", "--", "", "N", "NN", "./.")
) -> beagler_result
```

`sample_col = NULL` triggers auto-detection.

`$summary` includes:

- `n_samples`, `n_markers`, `inferred_encoding`, `inferred_ploidy`
- `missing_rate_overall`, `missing_rate_by_sample`, `missing_rate_by_marker`
- `monomorphic_markers`, `duplicated_samples`, `duplicated_markers`
- Mixed-encoding flag (e.g. some cells look like `AHB`, others like `012`)

### Exports

`validate_genotype_table()`.

### Exit criteria

- Tests for every supported encoding × every shape of malformed input.
- Snapshot tests on `print()` output.

---

## Milestone 3 — Marker-map validation

**Goal:** validate the marker metadata required for a real CSV → VCF conversion.

### Function

```r
validate_marker_map(
  map,
  genotype_table = NULL,
  require_alleles = FALSE,
  chrom_style = c("auto", "numeric", "chr_prefixed", "as_is")
) -> beagler_result
```

`genotype_table` (when provided) cross-checks marker coverage. `require_alleles = TRUE` demands `REF` and `ALT` columns in the map.

`$summary` includes:

- `n_markers_in_map`
- `markers_in_table_not_in_map`, `markers_in_map_not_in_table`
- `duplicated_marker_ids`, `duplicated_chrom_pos`
- `unsorted_within_chrom`
- `chromosome_summary` (counts and label style per chromosome)

### Exports

`validate_marker_map()`.

---

## Milestone 4 — CSV and HapMap → VCF

**Goal:** convert breeder genotype tables to valid VCF.

Excel input is **deferred to v0.2.0**. Document a "save as CSV" workflow.

### Functions

```r
csv_to_vcf(
  csv,
  marker_map,
  output_vcf,
  sample_col = NULL,
  encoding = "auto",
  ploidy = 2L,
  chrom_style = "auto",
  ahb_ref = "A",
  aabb_ref_homo = "AA",
  missing_codes = c(NA, ".", "--", "", "N", "NN", "./."),
  overwrite = FALSE
) -> beagler_result

hapmap_to_vcf(
  hmp,
  output_vcf,
  build = NA_character_,
  chrom_style = "auto",
  missing_codes = c(NA, ".", "--", "", "N", "NN"),
  overwrite = FALSE
) -> beagler_result
```

`csv_to_vcf()` requires `marker_map` because REF and ALT cannot be inferred from `A/H/B` or `0/1/2` without it.

`hapmap_to_vcf()` decodes IUPAC codes (`R, Y, S, W, K, M, N`) and uses HapMap's embedded `chrom` / `pos` / `alleles` columns. The optional `build` argument is recorded only in the VCF header.

### Exports

`csv_to_vcf()`, `hapmap_to_vcf()`.

### Exit criteria

- Roundtrip: CSV → VCF → re-read → genotypes match cell-for-cell.
- The day-one recoding test from M0 turns green here.
- HapMap fixture from a public TASSEL example.
- Output VCF passes `validate_vcf(target = "general")` from M5.

---

## Milestone 5 — VCF validation

**Goal:** one validator with multiple targets, replacing the planned `check_vcf` / `check_vcf_for_beagle` / `check_vcf_samples` / `check_vcf_markers` quartet.

### Function

```r
validate_vcf(
  vcf,
  target = c("general", "beagle", "plink"),
  max_records_to_scan = Inf
) -> beagler_result
```

Always-on checks:

- Required header lines and column count
- Sample columns present, no duplicates
- Marker IDs present, no duplicate IDs, no duplicated `(chrom, pos)`
- REF and ALT non-empty; biallelic where required
- Genotype field structure consistent with declared FORMAT
- Records sorted by `(chrom, pos)`
- Multiallelic sites flagged
- Non-DNA characters flagged

`target = "beagle"` adds:

- Diploid throughout (Beagle is diploid-only)
- File is BGZF (warn if not; error only when also passed as `ref=`)
- REF ≠ ALT

`target = "plink"` adds:

- Chromosome labels supported by PLINK
- Variant IDs unique

### Performance note

Full validation walks every record of the input VCF. On very large files (millions of variants on hundreds of samples) this can dominate runtime even when no errors are present. The `max_records_to_scan` argument exists today as a cheap escape hatch (header in full + a fixed number of records). A richer **sampled-mode validation** — full header + a deterministic sample of records, with the sample size and seed recorded in the result object — is on the v0.3.0 list. v0.1.0 and v0.2.0 ship full-scan validation only; users with very large files should rely on `max_records_to_scan` or stage validation against a per-chromosome split (M9, v0.2.0).

### Exports

`validate_vcf()`.

---

## Milestone 6 — VCF repair

**Goal:** fix common VCF issues that block Beagle, exposed as **one** function.

The recommended workflow is `validate_vcf()` → `repair_vcf()` → `validate_vcf()` → `run_beagle()`. M5's `validate_vcf(target = "beagle")` is the contract that defines what M6 must repair: every condition the validator reports as an error has a corresponding repair step here.

### Function

```r
repair_vcf(
  vcf,
  output_vcf,
  steps = c(
    "sort",
    "fix_marker_ids",
    "remove_duplicates",
    "filter_biallelic",
    "drop_invalid_records"
  ),
  chrom_style = "auto",
  overwrite = FALSE
) -> beagler_result
```

Each step is implemented by an internal helper:

```r
sort_vcf()
fix_marker_ids() # uses chrom_pos_ref_alt naming for missing IDs
remove_duplicate_markers()
filter_biallelic()
drop_invalid_records()
```

- Always writes a new file; never modifies in place.
- `$summary` reports per-step counts (`n_records_in`, `n_records_out`, `n_changed`).

### Exports

`repair_vcf()` only. The five step helpers stay internal.

### Exit criteria

- A repaired version of `inst/extdata/tiny_dirty.vcf` passes `validate_vcf(target = "beagle")`.
- Each step has a "does the right thing" test and an "idempotent on a clean file" test.
- For every error code the M5 validator can emit under `target = "beagle"`, either a repair step exists or the roadmap explicitly notes the case is out of scope (e.g. multiallelic split with custom REF/ALT logic — left to bcftools).

---

## Milestone 7 — Beagle execution (single VCF)

**Goal:** run Beagle from R safely, with full discovery, capture, and audit.

By the time a VCF reaches `run_beagle()` it has already passed `validate_vcf(target = "beagle")` (post-`repair_vcf()` if needed). M7 owns Java/JAR discovery, JVM heap configuration, the run record, and output indexing — **not** VCF cleanliness.

### Functions

```r
run_beagle(
  gt,
  out,
  ref = NULL,
  map = NULL,
  mode = c("impute", "phase"),
  nthreads = 1L,
  memory_gb = 4L,
  seed = NULL,
  beagle_jar = NULL,
  java = NULL,
  extra_args = character(),
  run_dir = NULL,
  overwrite = FALSE
) -> beagler_run

install_beagle(version = "5.4", dest = NULL)
beagle_version(beagle_jar = NULL)
```

`gt` is the target VCF path; `out` is the output prefix. `ref` is a pre-built reference panel VCF (see "Reference panels" below). `map` is a PLINK-format genetic map. `beagle_jar` and `java` are auto-discovered when `NULL`. `install_beagle()` downloads the pinned default Beagle 5.4 JAR.

### Reference panels

Reference panels (`ref=`) are supplied as **pre-built** VCF files; beagler does **not** construct reference panels in v0.1.0. A dedicated `validate_reference_match(target_vcf, reference_vcf)` — checking that chromosome labels, allele orientation, and build coordinates agree between target and panel — is on the **v0.2.0** list. Until it lands, target/reference compatibility is the user's responsibility.

### Behaviour

- `beagle_jar` resolution order: `beagle_jar` arg → `BEAGLE_PATH` env → `tools::R_user_dir("beagler", "cache")` → common locations.
- `java` resolution order: `java` arg → `JAVA_HOME/bin/java` → `Sys.which("java")`. Refuse Java < 8.
- Refuse Beagle < 5 (detected from the JAR's first stdout line).
- Build argv as a *vector*; never compose shell strings (correctness on Windows depends on this).
- Write three files into `run_dir`:
  - `stdout.log`, `stderr.log`
  - `run_record.json` — command, versions, args, input SHA-256, runtime, exit code
- After success, BGZIP-compress and `tabix`-index the output if `bcftools` / `tabix` are available; otherwise leave the Beagle-default `.vcf.gz` (BGZF already) and warn that indexing was skipped.
- Always return a `beagler_run`, success or failure. Non-zero status surfaces stderr text in `print()`.

### Exports

`run_beagle()`, `install_beagle()`, `beagle_version()`.

### Exit criteria

- Integration test: tiny VCF (≤ 100 markers, ≤ 10 samples) phased successfully on Linux and macOS CI.
- Negative-path test: invalid VCF → non-zero `$status`, stderr captured, useful `print()` output.
- `find_java()` and `find_beagle_jar()` unit-tested with `withr::with_envvar()`.
- Run record JSON validates against `inst/schemas/run_record.schema.json`.
- Polyploidy guard: VCF with non-diploid GT fields refused with an actionable message before Java is even launched.

---

## Milestone 8 — VCF ↔ PLINK and output extraction (deferred to v0.2.0)

**Goal:** practical interoperability and getting imputed results back into R.

### Functions

```r
vcf_to_plink(vcf, plink_prefix, plink = NULL, overwrite = FALSE) -> beagler_run
plink_to_vcf(plink_prefix, vcf, plink = NULL, overwrite = FALSE) -> beagler_run

vcf_to_dosage(vcf, samples = NULL, markers = NULL) -> matrix
vcf_to_hardcalls(
  vcf,
  samples = NULL,
  markers = NULL,
  encoding = c("012", "AHB", "AABB")
) -> matrix
```

### Out of scope here

`run_plink()`, `plink_filter()`, `plink_qc()` are **not** added. PLINK has hundreds of flags; we are not wrapping them.

### Exit criteria

- Roundtrip: VCF → PLINK → VCF preserves genotypes.
- Dosage extraction matches a hand-computed matrix on a 5 × 5 fixture.
- Tests skip cleanly when `plink` is not on PATH.

---

## Milestone 9 — Chromosome workflows (deferred to v0.2.0)

**Goal:** scale to large VCFs by splitting, running, and merging per chromosome.

### Functions

```r
split_vcf_by_chrom(
  vcf,
  out_dir,
  chrom_filter = NULL,
  bcftools = NULL
) -> list

run_beagle_by_chrom(
  gt,
  out_dir,
  ref = NULL,
  map = NULL,
  parallel = c("none", "future", "mclapply"),
  workers = 1L,
  ...
) -> list

merge_beagle_outputs(run_results, output_vcf, bcftools = NULL) -> beagler_run
```

`run_beagle_by_chrom()` returns a list with one `beagler_run` per chromosome.

### Behaviour

- Requires `bcftools` for split and merge (a slow R-only fallback exists for split, with a loud warning; not for merge).
- Per-chromosome runs each get their own `run_dir` and log file.
- **Resume logic:** if `run_dir/<chrom>/run_record.json` reports `status == 0` and the input SHA-256 still matches, skip that chromosome.

### Exit criteria

- 21-chromosome fixture run merges to a single VCF that passes `validate_vcf(target = "beagle")`.
- A failure injected on chromosome 7 leaves 1–6 reusable on rerun.

---

## Milestone 10 — Reports, vignettes, v0.1.0 cut

### Deliverables

- `print.beagler_run()` — one-screen run summary.
- `summary.beagler_result()` — one-screen QC summary.
- `imputation_report(run, file)` — one-page Markdown summary (no Rmd templating; `cat()`-style for v0.1.0).
- Vignettes:
  - `01_csv_to_vcf.Rmd`
  - `02_hapmap_to_vcf.Rmd`
  - `03_validate_for_beagle.Rmd`
  - `04_run_beagle.Rmd`
  - `05_extract_imputed_genotypes.Rmd` (uses M8; ships in v0.2.0)
- `pkgdown` site, GitHub Pages.
- `NEWS.md` updated with v0.1.0 entry.
- Tag and release **v0.1.0**.

### Exit criteria

- All vignettes build without `eval = FALSE` chunks (small bundled fixtures + skip predicates in chunk options).
- `R CMD check --as-cran` clean on all three OSes.
- Coverage ≥ 85 % on non-integration code.

---

## v0.1.0 target scope (the only thing that ships)

**Exported user-facing functions (7):**

```r
validate_genotype_table()
validate_marker_map()
validate_vcf()
csv_to_vcf()
hapmap_to_vcf()
repair_vcf()
run_beagle()
```

**Helpers exported but secondary:**

```r
install_beagle()
beagle_version()
```

**S3 methods:** `print`, `summary`, `format`, `as.data.frame` for `beagler_result` and `beagler_run`.

That is the entire public API of v0.1.0. Everything else is `@keywords internal`. `vcf_to_dosage()`, `vcf_to_hardcalls()`, the PLINK bridge, chromosome workflows, and `validate_reference_match()` all land in v0.2.0 with M8, M9, and the reference-compatibility checker.

---

## Internal helpers (not exported)

```r
# Process and CLI
run_cli()
build_beagle_cmd()
with_run_dir()

# Discovery
find_beagle_jar()
find_java()
detect_java_version()
detect_beagle_version()
find_bcftools()
find_tabix()
find_plink() # M8

# Path / IO
normalize_path_safe()
assert_readable()
assert_writable()
is_bgzf()
bgzip_in_place()
tabix_index()

# Genotype encoding
detect_genotype_encoding()
recode_to_vcf_gt() # replaces deleted recode_snp()
decode_iupac()
is_biallelic_row()

# Chromosome
normalize_chrom_label()
chrom_style()

# VCF
read_vcf_header()
parse_vcf_records()
write_vcf()
summarize_missingness()
detect_ploidy()

# Validation primitives
assert_unique()
assert_in_set()
collect_issues() # builder for beagler_result

# Repair primitives
sort_vcf()
fix_marker_ids()
remove_duplicate_markers()
filter_biallelic()
drop_invalid_records()

# Logging
log_step()
write_run_record_json()
```

All marked `@keywords internal` and absent from `NAMESPACE`. Internal status is enforced by `NAMESPACE` and roxygen, never by a leading dot.

---

## Deferred to later versions

| When | What |
|---|---|
| **v0.2.0** | M8 (PLINK bridge + dosage / hardcalls), M9 (chromosome workflows), `validate_reference_match()` (target ↔ reference VCF compatibility check), Excel input, CRAN submission |
| **v0.3.0** | Marker-map merging across builds, masked-genotype benchmarking (`mask_genotypes`, `benchmark_imputation`, `compare_imputed_vs_true`), templated Rmd reports, sampled-mode VCF validation |
| **Under review** | PLINK QC wrappers; reference-panel **construction** helpers (no v0.x commitment — `bcftools` is the right tool today, and beagler's role is to *validate* compatibility with `validate_reference_match()`, not to build panels) |

---

## Testing strategy

### Layers

1. **Unit tests** — pure R, no external tools. ≥ 85 % coverage. Run on every push.
2. **Roundtrip tests** — table → VCF → re-read; VCF → PLINK → VCF; dosage matches expectation.
3. **Integration tests** — Beagle, PLINK, bcftools. Skip cleanly via `skip_if_no_*` when tools are missing. Run in a dedicated CI job that installs and caches them.
4. **Snapshot tests** — `print` / `summary` output for both result classes; user-facing error messages (so message changes are deliberate, not accidental).

### Fixtures

`inst/extdata/` (shipped):

```
tiny.csv                     50 markers × 10 samples
tiny_map.tsv                 marker map for tiny.csv
tiny.hmp.txt                 HapMap equivalent
tiny.vcf                     clean, validates as 'beagle' target
tiny_dirty.vcf               deliberately broken; used by repair tests
tiny_imputed_expected.vcf    small reference output for regression
```

`tests/testthat/fixtures/` (not shipped):

```
medium_chrom.vcf.gz          3 chromosomes, ~5k markers (M9)
```

### CI matrix

- **Lint** — `lintr` on every push.
- **R-CMD-check** — `{ubuntu, macos, windows} × {release, devel}`. CRAN deps only.
- **Integration** — Linux only, runs after R-CMD-check passes. Caches Beagle JAR, PLINK, bcftools. Allowed to fail on `devel`.
- **Coverage** — `covr::codecov()` on Linux release.

### Day-one canary (in M0)

`tests/testthat/test-genotype-recoding.R` is written *before* any new feature code and asserts the correct mapping table for both the `0/1/2` and `-1/0/1` encodings (the bug in the deleted `recode_snp()`). It fails until M4 lands and turns green there. It is the package's first proof of life.

---

## Risks and assumptions log

| # | Assumption | Risk if wrong | Mitigation |
|---|---|---|---|
| 1 | Pure R + shelled-out `bcftools` / `tabix` is fast enough | Slow on million-marker VCFs | Re-evaluate `Rsamtools` / `VariantAnnotation` for v0.2.0 |
| 2 | `processx::run()` quoting is correct on Windows | Beagle launch fails on Windows | Argv-vector calls only; Windows tested in CI |
| 3 | Users have Java 8+ | Beagle won't start | Detect early, refuse with clear install message |
| 4 | Users provide a marker map for `csv_to_vcf()` | REF / ALT cannot be inferred | `marker_map` is a required argument; fail loud |
| 5 | Users' data is diploid | Polyploid genotypes break Beagle silently | Detect ploidy in `validate_vcf()` and refuse |
| 6 | `BEAGLE_PATH` points at Beagle 5.x | Beagle 4.x has different args | Detect version, refuse < 5 |
| 7 | Chromosome labels are clean | Mixed `chr1` / `1` styles silently break joins | `normalize_chrom_label()` used in every conversion path |
| 8 | CRAN allows our dependency set | Submission rejected | Stay GitHub-only for v0.1.0; Bioconductor option held in reserve |
| 9 | Beagle 5.4 download URL stays stable | `install_beagle()` fails behind firewalls or after upstream URL churn | Document manual install; cache once downloaded; revisit pinned version each minor release |
| 10 | `bcftools` is widely available | Index step skipped, downstream tools complain | Warn loudly; document install per OS |
| 11 | Full VCF validation is fast enough for v0.1.0 audiences | Validators dominate runtime on huge files | `max_records_to_scan` available; sampled-mode validation deferred to v0.3.0 |

---

## Suggested execution order

The development order mirrors the user-facing workflow `validate_vcf()` → `repair_vcf()` → `validate_vcf()` → `run_beagle()`:

1. **M0** — Pivot cleanup, CI bring-up, day-one canary test
2. **M1** — Result classes, CLI runner, discovery scaffolding
3. **M2** — Genotype-table validator
4. **M3** — Marker-map validator
5. **M4** — `csv_to_vcf()` and `hapmap_to_vcf()`
6. **M5** — `validate_vcf()` (general / beagle / plink targets)
7. **M6** — `repair_vcf()` driven by M5's `target = "beagle"` checks
8. **M7** — `run_beagle()` with discovery, JVM heap, run record
9. **M10** — Reports, vignettes, **release v0.1.0**
10. **M8** — PLINK bridge + dosage / hardcalls → v0.2.0
11. **M9** — Chromosome workflows → v0.2.0
12. **`validate_reference_match()`** — bundled with v0.2.0

---

## Appendix A — Files removed during M0

```
src/beagler.so
src/entrypoint.c
src/entrypoint.o
src/Makevars
src/Makevars.ucrt
src/Makevars.win
src/beagler-win.def
src/.gitignore
src/rust/Cargo.toml
src/rust/Cargo.lock
src/rust/src/lib.rs
src/rust/target/                           (entire tree)
R/extendr-wrappers.R
R/convert_csv_to_vcf.R
man/run_beagle_internal.Rd
```

`R/beagle.R` is rewritten, not deleted; `man/run_beagle.Rd` is regenerated by roxygen2.

## Appendix B — Files added during M0

```
NEWS.md
.github/workflows/R-CMD-check.yaml
.github/workflows/lint.yaml
tests/testthat.R
tests/testthat/helper-skip.R
tests/testthat/test-genotype-recoding.R
inst/schemas/run_record.schema.json        (placeholder; populated in M1)
```

---

*Last updated: 2026-04-28*
