#' Logging and run record (internal)
#' @name utils-logging
#' @keywords internal
NULL

#' Append a JSON Lines entry to `<run_dir>/run.log`
#'
#' Each call writes one line containing `time`, `step`, and any extra
#' key/value pairs passed via `...`.
#'
#' @keywords internal
#' @noRd
log_step <- function(run_dir, step, ...) {
  if (length(run_dir) != 1L || !is.character(run_dir) || is.na(run_dir)) {
    cli::cli_abort(
      "{.field run_dir} must be a single non-NA character string.",
      class = "galick_error_invalid_argument"
    )
  }
  if (!dir.exists(run_dir)) {
    cli::cli_abort(
      "{.field run_dir} does not exist: {.file {run_dir}}",
      class = "galick_error_invalid_argument"
    )
  }
  if (length(step) != 1L || !is.character(step) || is.na(step)) {
    cli::cli_abort(
      "{.field step} must be a single non-NA character string.",
      class = "galick_error_invalid_argument"
    )
  }

  log_path <- file.path(run_dir, "run.log")
  entry <- c(
    list(
      time = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
      step = step
    ),
    list(...)
  )

  json <- jsonlite::toJSON(entry, auto_unbox = TRUE, null = "null", na = "null")
  cat(json, "\n", file = log_path, append = TRUE, sep = "")
  invisible(NULL)
}

#' Write a versioned run record JSON to `<run_dir>/run_record.json`
#'
#' The record is validated against the v2.0.0 contract before writing; an
#' invalid record raises `galick_error_invalid_run_record`.
#'
#' @keywords internal
#' @noRd
write_run_record_json <- function(run_dir,
                                  command,
                                  versions,
                                  args = list(),
                                  inputs = list(),
                                  outputs = list(),
                                  runtime_sec,
                                  status,
                                  started_at,
                                  ended_at) {
  if (!dir.exists(run_dir)) {
    cli::cli_abort(
      "{.field run_dir} does not exist: {.file {run_dir}}",
      class = "galick_error_invalid_argument"
    )
  }

  record <- list(
    schema_version = "2.0.0",
    command = command,
    versions = as.list(versions),
    args = as.list(args),
    inputs = as.list(inputs),
    outputs = as.list(outputs),
    runtime_sec = runtime_sec,
    status = status,
    started_at = as.character(started_at),
    ended_at = as.character(ended_at)
  )

  validate_run_record(record)
  record$command <- I(record$command)
  if (!length(record$args)) record$args <- stats::setNames(list(), character())

  out_path <- file.path(run_dir, "run_record.json")
  json <- jsonlite::toJSON(
    record,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  writeLines(json, out_path)
  invisible(out_path)
}

#' Validator for the Galick run-record schema (v2.0.0)
#'
#' Mirrors the required-fields and type rules from
#' `inst/schemas/run_record.schema.json`. We keep this in lockstep with the
#' JSON Schema; if either drifts, tests catch it.
#'
#' @keywords internal
#' @noRd
validate_run_record <- function(record) {
  required <- c(
    "schema_version", "command", "versions",
    "runtime_sec", "status", "started_at", "ended_at"
  )
  missing_fields <- setdiff(required, names(record))
  if (length(missing_fields) > 0L) {
    cli::cli_abort(
      "Run record missing required fields: {.field {missing_fields}}",
      class = "galick_error_invalid_run_record"
    )
  }

  if (!identical(record$schema_version, "2.0.0")) {
    cli::cli_abort(
      "Run record {.field schema_version} must be {.val 2.0.0}; got {.val {record$schema_version}}.",
      class = "galick_error_invalid_run_record"
    )
  }

  if (!is.character(record$command) || length(record$command) < 1L || anyNA(record$command)) {
    cli::cli_abort(
      "Run record {.field command} must be a non-empty character vector.",
      class = "galick_error_invalid_run_record"
    )
  }

  if (!is.list(record$versions)) {
    cli::cli_abort(
      "Run record {.field versions} must be a list.",
      class = "galick_error_invalid_run_record"
    )
  }

  versions_required <- c("galick", "r", "os")
  versions_missing <- setdiff(versions_required, names(record$versions))
  if (length(versions_missing) > 0L) {
    cli::cli_abort(
      "Run record {.field versions} missing required fields: {.field {versions_missing}}",
      class = "galick_error_invalid_run_record"
    )
  }

  scalar_text <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (any(!vapply(record$versions, scalar_text, logical(1)))) {
    cli::cli_abort("Run-record versions must be nonempty strings.",
      class = "galick_error_invalid_run_record"
    )
  }

  if (!is.numeric(record$runtime_sec) || length(record$runtime_sec) != 1L ||
        !is.finite(record$runtime_sec) || record$runtime_sec < 0) {
    cli::cli_abort(
      "Run record {.field runtime_sec} must be a non-negative number.",
      class = "galick_error_invalid_run_record"
    )
  }

  if (!is.numeric(record$status) || length(record$status) != 1L || !is.finite(record$status) ||
        record$status != floor(record$status) || abs(record$status) > .Machine$integer.max) {
    cli::cli_abort(
      "Run record {.field status} must be an integer.",
      class = "galick_error_invalid_run_record"
    )
  }

  for (field in c("started_at", "ended_at")) {
    value <- record[[field]]
    if (!scalar_text(value) ||
          !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}([.][0-9]+)?Z$", value) ||
          is.na(as.POSIXct(value, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))) {
      cli::cli_abort("Run-record timestamps must be UTC ISO 8601 strings.",
        class = "galick_error_invalid_run_record"
      )
    }
  }
  if (as.POSIXct(record$ended_at, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC") <
        as.POSIXct(record$started_at, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")) {
    cli::cli_abort("Run-record end time precedes start time.", class = "galick_error_invalid_run_record")
  }
  bad_names <- is.null(names(record$args)) || any(!nzchar(names(record$args))) ||
    anyDuplicated(names(record$args))
  if (!is.list(record$args) || (length(record$args) && bad_names)) {
    cli::cli_abort("Run-record args must be a named list.", class = "galick_error_invalid_run_record")
  }
  for (field in c("inputs", "outputs")) {
    entries <- record[[field]]
    if (!is.list(entries)) {
      cli::cli_abort("Run-record file entries must be lists.", class = "galick_error_invalid_run_record")
    }
    for (entry in entries) {
      valid <- is.list(entry) && scalar_text(entry$path)
      hash <- entry$sha256
      if (!valid || (field == "inputs" && is.null(hash)) ||
            (!is.null(hash) && (!scalar_text(hash) || !grepl("^[0-9a-f]{64}$", hash)))) {
        cli::cli_abort("Run-record paths and SHA-256 hashes are invalid.",
          class = "galick_error_invalid_run_record"
        )
      }
    }
  }

  invisible(TRUE)
}
