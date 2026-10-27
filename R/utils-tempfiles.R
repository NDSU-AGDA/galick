#' Run-directory lifecycle (internal)
#' @name utils-tempfiles
#' @keywords internal
NULL

#' Evaluate an expression with `getOption("galick.run_dir")` set
#'
#' If `dir = NULL`, a fresh temp directory is created. The directory is *not*
#' deleted on exit; the caller decides its fate. The previous value of the
#' `galick.run_dir` option is restored even if `expr` errors.
#'
#' @keywords internal
#' @noRd
with_run_dir <- function(dir = NULL, expr) {
  if (is.null(dir)) {
    dir <- tempfile("galick_run_")
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  abs_dir <- normalize_path_safe(dir, must_exist = TRUE)

  withr::local_options(list(galick.run_dir = abs_dir))
  force(expr)
}
