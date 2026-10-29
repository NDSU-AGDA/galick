## Submission status

Preparation for the first CRAN submission of galick 0.1.0. Not submitted.
Windows reports must be reviewed before submission.

## Test environments

- macOS 27.0.1, arm64, R 4.5.3.
- Windows R-release 4.6.1: the 2026-10-28 build exposed a PuTTY/genetics PLINK
  executable-name collision. Discovery now verifies the genetics version banner;
  tests skip an unavailable genetics executable rather than invoking PuTTY.
  The prior Windows report does not validate the corrected source.
- Updated Windows R-release and R-devel checks were uploaded on 2026-10-29; review the
  new reports at the maintainer email address before submission.

## R CMD check --as-cran

0 errors | 0 warnings | 1 note

The only note is:

    New submission

All suggested packages were installed and _R_CHECK_FORCE_SUGGESTS_ was true.
PDF and HTML manuals, examples, tests and vignette rebuilding passed.
The external system-clock probe was disabled, matching the CI check environment.

## R CMD check

0 errors | 0 warnings | 0 notes

The ordinary check used the same source tarball, all suggested packages and
available integration tools, including Java/Beagle. PDF manual generation and
vignette rebuilding passed.

## Additional checks

- Full test suite with Java/Beagle, PLINK, bcftools, bgzip and tabix available:
  1,647 assertions passed; 0 failures, 0 warnings, 0 skips. The opt-in installer
  download test was enabled and passed using a temporary destination.
- All 28 exported functions have help examples. The multi-job evaluation and
  marker-panel examples were explicitly enabled and exercised separately.
  The installer help example is marked dontrun; checks do not download software.
- lintr: no lints.
- urlchecker: all remaining package URLs valid.
- Website: mkdocs build --strict passed; 153 R fences across 36 documentation
  files parsed and checked against public function arguments. A separate static
  check validates canonical URLs, unique metadata, sitemap and search coverage
  for all 30 pages and 28 exported functions.

## External software and access

Java and Beagle are required only for Beagle workflows. External software is
not bundled or downloaded during package installation, loading or ordinary
checks. Integration tests and examples are conditional on available tools;
longer examples and the network installer test require explicit opt-in.

The repository and issue tracker at NDSU-AGDA/galick are now public and listed
in DESCRIPTION. README links use public absolute URLs rather than file URIs
pointing to files excluded from the source package. The previous Windows
spelling note identified valid domain vocabulary: HapMap, PLINK, VCF/VCFs and
autotetraploid. These terms describe file formats, external software and an
explicit inheritance model, not typographical errors.

## Checked source archive

The corrected source tarball was uploaded to both win-builder targets.
SHA-256:

    674098e205490966cc1ce97b810d6c13f801164daa77639ba349a0a569e51932

Only instructional fixtures are bundled; local research VCF files are excluded.
