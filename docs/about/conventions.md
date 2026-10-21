# Project conventions

These conventions describe the package's current development practices.

## Object model

**S3 only.** No R6, no S4. All beagler objects are S3 classes built
with `structure(list(...), class = c("beagler_xxx", "list"))`.

## Tests

- Every feature, refactor, and bug fix ships with tests.
  **testthat 3rd edition.**
- Tests live in `tests/testthat/test-<topic>.R`; helpers in
  `tests/testthat/helper-*.R`.
- Bug fixes: write a failing test first, then the fix.
- Integration tests that need external tools (Beagle, PLINK,
  bcftools, tabix, java) call the matching `skip_if_no_*()` helper
  at the top.

## Naming

- **No leading dot** in function or variable names. Use `find_java()`,
  not `.find_java()`. Internal status is conferred by roxygen
  `@keywords internal` and absence from `NAMESPACE`, never by a
  leading dot.
- The only legitimate `.` in a function name is **S3 method
  dispatch**: `print.beagler_run`, `summary.beagler_result`,
  `as.data.frame.beagler_run`.
- `snake_case` for everything: functions, arguments, locals, classes,
  file names.

## Style (R source code)

Follow the [tidyverse style guide](https://style.tidyverse.org/) with
these specifics:

- **Native pipe `|>`.** Don't use `%>%`.
- **Vertical pipe layout** — one operation per line.
- **No alignment whitespace.** Single space around `<-` and `=`.
  Same rule applies inside argument lists — never pad `=` signs into
  a column.
- **Comments only where the logic isn't self-evident.** No "what"
  comments. Reserve comments for the *why*: non-obvious invariants,
  workarounds for upstream bugs, cross-references to Beagle quirks.
- Double quotes for strings.
- The lint gate uses a 120-character line limit, including VCF fixture strings.

## Errors and warnings

- User-facing errors via `cli::cli_abort()` with a `class` attribute
  (e.g. `beagler_error_invalid_encoding`) so users and tests can
  dispatch on them.
- User-facing warnings via `cli::cli_warn()`, also classed.
- Validation issues are *collected* into the `$warnings` slot of
  `beagler_result`, never raised as live R warnings.

## Subprocess execution

- External tools (`java`, `beagle.jar`, `plink`, `bcftools`, `tabix`,
  `bgzip`) are invoked through `processx::run()` (preferred) or
  `system2()` with the argument list as a **vector**.
- Never compose shell command strings via `paste()` / `sprintf()`
  and hand them to `system()`. Windows quoting and paths with spaces
  will silently break.

## Documentation

- roxygen2 for every function in `R/`.
- `@export` only the user-facing API. The v0.1.0 public surface is
  currently contains 27 exported functions.
- Runnable examples by default; use `\dontrun{}` only when an
  example genuinely requires an external tool.

## Dependencies

- v0.1.0 R dependencies are CRAN-only: `cli`, `data.table`,
  `digest`, `jsonlite`, `processx`, `tibble`, `vroom`, `withr`.
- Adding a dependency requires a `DESCRIPTION` update and a roadmap
  note.
- For VCF / BGZF I/O, prefer shelling out to `bcftools` / `tabix` /
  `bgzip` over reimplementing parsers in R.
