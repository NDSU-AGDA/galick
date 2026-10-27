#' Path utilities (internal)
#'
#' @name utils-paths
#' @keywords internal
NULL

#' Normalize a path to an absolute form
#' @keywords internal
#' @noRd
normalize_path_safe <- function(path, must_exist = FALSE) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    cli::cli_abort(
      "Path must be a single non-NA character string.",
      class = "galick_error_invalid_path"
    )
  }
  expanded <- path.expand(path)
  normalizePath(expanded, mustWork = if (must_exist) TRUE else FALSE)
}

#' Assert a path is readable
#' @keywords internal
#' @noRd
assert_readable <- function(path) {
  abs <- normalize_path_safe(path, must_exist = TRUE)
  if (file.access(abs, mode = 4L) != 0L) {
    cli::cli_abort(
      "Path is not readable: {.file {abs}}",
      class = "galick_error_path_not_readable"
    )
  }
  invisible(abs)
}

#' Assert a path (or its parent directory) is writable
#' @keywords internal
#' @noRd
assert_writable <- function(path) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    cli::cli_abort(
      "Path must be a single non-NA character string.",
      class = "galick_error_invalid_path"
    )
  }
  expanded <- path.expand(path)
  if (file.exists(expanded)) {
    if (file.access(expanded, mode = 2L) != 0L) {
      cli::cli_abort(
        "Path is not writable: {.file {expanded}}",
        class = "galick_error_path_not_writable"
      )
    }
  } else {
    parent <- dirname(expanded)
    if (!dir.exists(parent) || file.access(parent, mode = 2L) != 0L) {
      cli::cli_abort(
        "Parent directory is not writable: {.file {parent}}",
        class = "galick_error_path_not_writable"
      )
    }
  }
  invisible(expanded)
}
