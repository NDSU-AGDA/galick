#' galick_run S3 class (internal)
#' @name galick-run
#' @keywords internal
NULL

#' Construct a galick_run
#'
#' Captures the result of a single external-tool invocation: status,
#' runtime, the argv that ran, captured log paths, and the path to the
#' run-record JSON.
#'
#' @keywords internal
#' @noRd
new_galick_run <- function(status = NA_integer_,
                           runtime_sec = NA_real_,
                           command = character(),
                           version = list(),
                           output_paths = character(),
                           stdout_log = NA_character_,
                           stderr_log = NA_character_,
                           run_record = NA_character_,
                           call = NULL) {
  metadata <- list(
    call = call,
    package_version = utils::packageVersion("galick"),
    time = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  )

  structure(
    list(
      status = as.integer(status),
      runtime_sec = as.numeric(runtime_sec),
      command = as.character(command),
      version = as.list(version),
      output_paths = as.character(output_paths),
      stdout_log = as.character(stdout_log),
      stderr_log = as.character(stderr_log),
      run_record = as.character(run_record),
      metadata = metadata
    ),
    class = c("galick_run", "list")
  )
}

#' @export
print.galick_run <- function(x, ...) {
  status_str <- if (is.na(x$status)) {
    "status NA"
  } else if (x$status == 0L) {
    "OK"
  } else {
    sprintf("status %d", x$status)
  }
  cli::cli_h1(sprintf("galick_run: %s", status_str))

  rt <- if (is.na(x$runtime_sec)) "NA" else sprintf("%.2f sec", x$runtime_sec)
  cli::cli_text(sprintf("Runtime: %s", rt))

  if (length(x$command) > 0L) {
    cli::cli_text("Command: {paste(x$command, collapse = ' ')}")
  }

  outs <- x$output_paths[!is.na(x$output_paths)]
  if (length(outs) > 0L) {
    cli::cli_h3("Outputs")
    for (p in outs) cli::cli_li("{.file {p}}")
  }

  has_stdout_log <- length(x$stdout_log) > 0L && !is.na(x$stdout_log[1])
  has_stderr_log <- length(x$stderr_log) > 0L && !is.na(x$stderr_log[1])
  if (has_stdout_log || has_stderr_log) {
    cli::cli_h3("Logs")
    if (has_stdout_log) cli::cli_li("stdout: {.file {x$stdout_log}}")
    if (has_stderr_log) cli::cli_li("stderr: {.file {x$stderr_log}}")
  }

  invisible(x)
}

#' @export
summary.galick_run <- function(object, ...) {
  print(object)
  failed <- !is.na(object$status) && object$status != 0L
  has_stderr <- length(object$stderr_log) > 0L &&
    !is.na(object$stderr_log[1]) &&
    file.exists(object$stderr_log[1])
  if (failed && has_stderr) {
    cli::cli_h3("stderr (last 20 lines)")
    err_lines <- utils::tail(readLines(object$stderr_log[1], warn = FALSE), 20L)
    for (line in err_lines) cli::cli_text("{line}")
  }
  invisible(object)
}

#' @export
format.galick_run <- function(x, ...) {
  sprintf(
    "<galick_run: status=%s, runtime=%.2fs>",
    as.character(x$status),
    x$runtime_sec
  )
}

#' @export
as.data.frame.galick_run <- function(x, ...) {
  data.frame(
    status = x$status,
    runtime_sec = x$runtime_sec,
    n_outputs = length(x$output_paths),
    stdout_log = if (length(x$stdout_log) > 0L) x$stdout_log[1] else NA_character_,
    stderr_log = if (length(x$stderr_log) > 0L) x$stderr_log[1] else NA_character_,
    stringsAsFactors = FALSE
  )
}
