#' Portable SHA-256 file checksums
#' @keywords internal
#' @noRd
file_sha256 <- function(paths) {
  vapply(paths, function(path) {
    assert_readable(path)
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
  }, character(1), USE.NAMES = FALSE)
}
