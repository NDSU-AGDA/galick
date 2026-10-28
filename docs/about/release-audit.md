# Release audit

## Scope

The Galick rename was checked against the source namespace, S3 dispatch,
condition classes, installed resources, JSON records, tests, vignettes and
MkDocs website. The 28 public function names are unchanged.

## Findings addressed

| Finding | Resolution |
|---|---|
| Package identity was embedded in classes, conditions, caches, options and fixtures. | Renamed consistently to `galick`, with an explicit migration guide. |
| Website omitted recent exports and misdescribed result types. | Updated navigation, complete API signatures, experimental subclasses and examples. |
| Homepage could use a nonexistent repaired file. | Replaced with an executable bundled-data workflow and explicit status checks. |
| Static example output differed from fixtures or implied accuracy evidence. | Removed hand-written outputs; installed vignettes generate their own results. |
| GT counts and expected DS dosage could be confused. | Clarified extraction and explicit DS export, with no fallback claim. |
| Reference examples incorrectly required an indexed compressed panel. | Matched the actual runner's plain/compressed reference contract. |
| Renamed version key changes the run-record contract. | Bumped run-record schema to 2.0.0; other experiment schemas remain 1.0.0. |
| Record serialization could emit a scalar command or an array for empty args. | Preserved JSON types and added regression tests for malformed record fields. |
| Repair output promotion did not check the final file move. | Reused transactional file publication with rollback. |
| PLINK could report stale or partial output after a failed overwrite. | Stage fresh outputs, validate the set, then publish with rollback. |
| Help pages had no executable examples. | Added examples for all 28 exports, a small phased-founder fixture and documentation regression tests. |
| Repository and installation references still used the previous owner or feature branch. | Updated to `NDSU-AGDA/galick` and the default branch. |
| Package metadata advertised unavailable public documentation and issue URLs. | Omitted inaccessible URLs while the repository remains private; installed vignettes remain available. |

## Reproducible checks

- `Rscript tools/build-reference.R` regenerates public API declarations.
- `Rscript tools/check-documentation.R` parses website/README R fences and
  checks named arguments against the current public functions.
- Add `--run` to execute the curated bundled-data examples. Beagle sections
  run only with a configured JAR and Java; installation snippets are never run.
- `mkdocs build --strict` checks the website build and internal file links.
- `devtools::test()` runs regressions and available external-tool integration.
- `R CMD check` validates the built source package, documentation and vignettes.

The audit does not establish a global algorithmic optimum, biological prediction
accuracy, regulatory compliance or genome-scale performance.

## Remaining boundaries

The installer checks the downloaded JAR's version, not a published checksum.
PLINK bridge objects do not currently retain separate JSON run records.
Most VCF operations are in-memory. General validation is not full conformance
testing for every symbolic structural-variant representation.
HapMap conversion uses the declared alleles and does not reverse-complement
calls from strand metadata or validate them against a reference FASTA.

Real-data experiments with appropriate family/generation holdouts remain needed
before making biological superiority or breeding-value claims. The repository
is now `NDSU-AGDA/galick`; public hosting and the release tag remain separate steps.

## CRAN preparation

Checks are recorded in the repository's `cran-comments.md`, which is excluded
from the source package. All 28 help examples have been exercised, including
the optional multi-job examples with `GALICK_RUN_EXAMPLES=true` and configured
Java/Beagle. The installer example is deliberately not executed; its separate
network regression test runs only when explicitly enabled.

Windows R-release and R-devel builds have been requested from win-builder.
Review the reports emailed to the maintainer before CRAN submission. Uploading
to win-builder does not submit the package to CRAN.
