## Submission status

Preparation for the first CRAN submission of galick 0.1.0. Not submitted.
Windows reports must be reviewed before submission.

## Test environments

- macOS 27.0.1, arm64, R 4.5.3.
- Windows R-release and R-devel: uploaded to win-builder on 2026-10-28;
  reports pending at the maintainer email address in DESCRIPTION.

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
  1,638 assertions passed; 0 failures, 0 warnings, 0 skips. The opt-in installer
  download test was enabled and passed using a temporary destination.
- All 28 exported functions have help examples. The multi-job evaluation and
  marker-panel examples were explicitly enabled and exercised separately.
  The installer help example is marked dontrun; checks do not download software.
- lintr: no lints.
- urlchecker: all remaining package URLs valid.
- Website: mkdocs build --strict passed; 153 R fences across 34 documentation
  files parsed and checked against public function arguments.

## External software and access

Java and Beagle are required only for Beagle workflows. External software is
not bundled or downloaded during package installation, loading or ordinary
checks. Integration tests and examples are conditional on available tools;
longer examples and the network installer test require explicit opt-in.

The repository is NDSU-AGDA/galick and currently requires access. Unavailable
public website and issue-tracker URLs have been omitted from DESCRIPTION.
Installed vignettes provide the tutorials without requiring repository access.

## Checked source archive

The same source tarball was uploaded to both win-builder targets.
SHA-256:

    1d7d3692ed3e5b43573124ed8f8cb9d9e3df004232336197e759dd04370203ad

Only instructional fixtures are bundled; local research VCF files are excluded.
