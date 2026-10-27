#' Subprocess runner (internal)
#' @name utils-cli
#' @keywords internal
NULL

#' Run an external command and capture its output
#'
#' Capture-then-write semantics: stdout/stderr are captured in memory by
#' [processx::run()] and *then* written to `stdout_log` / `stderr_log` if
#' provided. Real-time tee would require [processx::process()] and is not
#' part of the M1 contract.
#'
#' Argv must be a vector. The function never composes shell strings.
#'
#' @keywords internal
#' @noRd
run_cli <- function(cmd,
                    args = character(),
                    env = NULL,
                    timeout = Inf,
                    stdout_log = NULL,
                    stderr_log = NULL) {
  if (length(cmd) != 1L || !is.character(cmd) || is.na(cmd)) {
    cli::cli_abort(
      "{.field cmd} must be a single non-NA character string.",
      class = "galick_error_invalid_argument"
    )
  }

  resolved <- Sys.which(cmd)
  if (resolved == "" || is.na(resolved)) {
    if (file.exists(cmd)) {
      resolved <- normalize_path_safe(cmd, must_exist = TRUE)
    } else {
      cli::cli_abort(
        "Executable not found: {.field {cmd}}",
        class = "galick_error_executable_not_found"
      )
    }
  }

  start_time <- Sys.time()
  result <- processx::run(
    command = unname(resolved),
    args = as.character(args),
    env = env,
    timeout = as.numeric(timeout),
    error_on_status = FALSE
  )
  end_time <- Sys.time()
  runtime_sec <- as.numeric(difftime(end_time, start_time, units = "secs"))

  status <- if (is.null(result$status)) -1L else as.integer(result$status)
  stdout <- if (is.null(result$stdout)) "" else result$stdout
  stderr <- if (is.null(result$stderr)) "" else result$stderr

  if (isTRUE(result$timeout)) {
    status <- -1L
    stderr <- if (nzchar(stderr)) paste0(stderr, "\ntimeout reached") else "timeout reached"
  }

  if (!is.null(stdout_log)) {
    assert_writable(stdout_log)
    writeLines(stdout, stdout_log, sep = "")
  }
  if (!is.null(stderr_log)) {
    assert_writable(stderr_log)
    writeLines(stderr, stderr_log, sep = "")
  }

  list(
    status = status,
    runtime_sec = runtime_sec,
    stdout = stdout,
    stderr = stderr,
    command = c(unname(resolved), as.character(args))
  )
}
