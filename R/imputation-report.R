#' Write a Markdown summary report of a `galick_run`
#'
#' Produces a one-page Markdown file capturing status, runtime, command,
#' input/output paths, version block, and log/run-record locations.
#' Pulls SHA-256 hashes from the run record JSON when present.
#'
#' @param run A `galick_run` returned by [run_beagle()].
#' @param file Output path (typically ending in `.md`).
#' @return The output path, invisibly.
#' @examplesIf file.exists(Sys.getenv("BEAGLE_PATH")) && nzchar(Sys.which("java"))
#' vcf <- system.file("extdata", "tiny_founders.vcf", package = "galick")
#' directory <- tempfile("report_")
#' dir.create(directory)
#' run <- run_beagle(vcf, file.path(directory, "phased"), mode = "phase",
#'                   nthreads = 1L, memory_gb = 1L,
#'                   run_dir = file.path(directory, "logs"))
#' report <- file.path(directory, "report.md")
#' imputation_report(run, report)
#' readLines(report, n = 5L)
#' unlink(directory, recursive = TRUE)
#' @export
imputation_report <- function(run, file) {
  if (!inherits(run, "galick_run")) {
    cli::cli_abort(
      "{.arg run} must be a {.cls galick_run} object.",
      class = "galick_error_imputation_report_invalid_input"
    )
  }
  if (length(file) != 1L || !is.character(file) || is.na(file)) {
    cli::cli_abort(
      "{.arg file} must be a single non-NA character string.",
      class = "galick_error_imputation_report_invalid_input"
    )
  }

  rec <- if (!is.na(run$run_record) && file.exists(run$run_record)) {
    jsonlite::fromJSON(run$run_record, simplifyVector = FALSE)
  } else {
    list()
  }

  status_str <- if (is.na(run$status)) {
    "status NA"
  } else if (run$status == 0L) {
    "OK (0)"
  } else {
    sprintf("FAILED (%d)", run$status)
  }
  runtime_str <- if (is.na(run$runtime_sec)) {
    "NA"
  } else {
    sprintf("%.2f sec", run$runtime_sec)
  }

  inputs_block <- if (!is.null(rec$inputs) && length(rec$inputs) > 0L) {
    vapply(rec$inputs, function(inp) {
      sprintf("- `%s` (sha256: `%s`)", inp$path, substr(inp$sha256, 1, 12))
    }, character(1))
  } else {
    "_(none recorded)_"
  }
  outputs_block <- if (length(run$output_paths) > 0L && !all(is.na(run$output_paths))) {
    vapply(run$output_paths, function(p) sprintf("- `%s`", p), character(1))
  } else {
    "_(none written)_"
  }
  versions_block <- if (length(run$version) > 0L) {
    vapply(names(run$version), function(k) {
      val <- run$version[[k]]
      sprintf("- **%s**: %s", k, if (is.null(val)) "n/a" else as.character(val))
    }, character(1))
  } else {
    "_(none recorded)_"
  }

  lines <- c(
    "# galick run report",
    "",
    sprintf("**Status:** %s  |  **Runtime:** %s", status_str, runtime_str),
    "",
    sprintf("_Generated: %s_", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    "",
    "## Command",
    "",
    "```",
    paste(run$command, collapse = " "),
    "```",
    "",
    "## Inputs",
    "",
    inputs_block,
    "",
    "## Outputs",
    "",
    outputs_block,
    "",
    "## Versions",
    "",
    versions_block,
    "",
    "## Logs",
    "",
    sprintf(
      "- stdout: `%s`",
      if (length(run$stdout_log) > 0L) run$stdout_log[1] else "n/a"
    ),
    sprintf(
      "- stderr: `%s`",
      if (length(run$stderr_log) > 0L) run$stderr_log[1] else "n/a"
    ),
    sprintf(
      "- run record: `%s`",
      if (length(run$run_record) > 0L) run$run_record[1] else "n/a"
    )
  )

  writeLines(lines, file)
  invisible(file)
}
