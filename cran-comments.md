## Submission status

Resubmission preparation for galick 0.1.0 following CRAN feedback.
The revised archive has not yet been resubmitted.

## Changes requested by CRAN

- Enclosed 'HapMap', 'Beagle' and 'PLINK' in single quotes in Description.
  Title contains no software names and is unchanged.
- Added Browning, Zhou and Browning (2018)
  <doi:10.1016/j.ajhg.2018.07.015> as the reference for the external
  'Beagle' imputation method, not for the package's simulation or panel design.

The check results and archive checksums below describe previous checks,
before these metadata revisions.

## Test environments

- macOS 27.0.1, arm64, R 4.5.3.
- Windows Server 2022, R-release 4.6.2: 0 errors, 0 warnings, 1 NOTE.
  Report: https://win-builder.r-project.org/OHvHX53J40HA/00check.log
  The NOTE contains "New submission" and valid domain terms flagged by spelling
  checks (HapMap, PLINK, VCF/VCFs and autotetraploid). Tests, vignettes and both
  manual formats passed. The earlier PuTTY/genetics PLINK collision is resolved.
- R-devel report: pending. The Windows report above covers the archived source
  checksum below; subsequent documentation and diagnostic-wording changes are
  being checked locally.

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

- Latest full test suite with Java/Beagle, PLINK, bcftools, bgzip and tabix:
  1,646 assertions passed; 0 failures, 0 warnings, 1 expected skip for the
  opt-in installer download test. That unchanged network test passed in the
  earlier preflight run using a temporary destination.
- All 28 exported functions have help examples. The multi-job evaluation and
  marker-panel examples were explicitly enabled and exercised separately.
  The installer help example is marked dontrun; checks do not download software.
- lintr: no lints.
- urlchecker: all remaining package URLs valid.
- Website: mkdocs build --strict passed; 154 R fences across 35 documentation
  files parsed and checked against public function arguments. A separate static
  check validates canonical URLs, unique metadata, sitemap and search coverage
  for all 29 pages and 28 exported functions.
- Executed 68 runnable tutorial blocks, including simulation, evaluation,
  panel design and sequential/multicore/future chromosome jobs. Installation
  commands and templates requiring user data are explicitly excluded.

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

Latest local R CMD check archive after the documentation audit (0 errors,
0 warnings, 0 notes), SHA-256:

    1877f3671a2c08d708e5bcf8c1eb140dfa1da7d3d55a703acfc04ee1f53409c8
