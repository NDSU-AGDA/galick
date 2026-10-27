#' Pinned tested Beagle 5.4 download URL
#'
#' Updated when needed. We deliberately pin a *tested* build rather than
#' tracking the latest upstream release, because Beagle's JAR filenames
#' include date+hash and change unpredictably. v0.1.0 ships against
#' `beagle.29Oct24.c8e.jar`.
#'
#' @keywords internal
#' @noRd
beagle_pinned_url <- function() {
  "https://faculty.washington.edu/browning/beagle/beagle.29Oct24.c8e.jar"
}

#' Download a pinned tested Beagle JAR to the user cache
#'
#' Beagle is distributed separately under its own license. This helper
#' downloads a pinned tested release into the user's galick cache (or a
#' custom destination) and verifies it by parsing the version line.
#'
#' v0.1.0 supports `version = "5.4"` only; other values raise. The pinned
#' build is `beagle.29Oct24.c8e.jar` from Brian Browning's official Beagle
#' page. We call this a "pinned tested 5.4 jar" rather than "the canonical
#' 5.4 release"; upstream may list a different build by the time you read
#' this. If your environment requires a different build, set `BEAGLE_PATH`
#' or pass `beagle_jar = ` to `run_beagle()`.
#'
#' @param version Currently must equal `"5.4"`.
#' @param dest Destination path. Defaults to
#'   `tools::R_user_dir("galick", "cache")/beagle.jar`.
#' @return The destination path, invisibly.
#' @export
install_beagle <- function(version = "5.4", dest = NULL) {
  if (!identical(version, "5.4")) {
    cli::cli_abort(
      c(
        "Only Beagle {.val 5.4} is supported in v0.1.0; got {.val {version}}.",
        i = "Set {.envvar BEAGLE_PATH} or pass {.arg beagle_jar} for other builds."
      ),
      class = "galick_error_beagle_unsupported_version"
    )
  }

  if (is.null(dest)) {
    dest_dir <- tools::R_user_dir("galick", "cache")
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    }
    dest <- file.path(dest_dir, "beagle.jar")
  } else {
    if (length(dest) != 1L || !is.character(dest) || is.na(dest)) {
      cli::cli_abort(
        "{.arg dest} must be a single non-NA character string.",
        class = "galick_error_invalid_argument"
      )
    }
    dest_dir <- dirname(dest)
    if (!dir.exists(dest_dir)) {
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    }
  }

  url <- beagle_pinned_url()
  status <- tryCatch(
    utils::download.file(url, dest, mode = "wb", quiet = TRUE),
    error = function(e) {
      cli::cli_abort(
        c(
          "Could not download Beagle from {.url {url}}.",
          i = conditionMessage(e),
          i = "Check network or download the JAR manually."
        ),
        class = "galick_error_beagle_download_failed"
      )
    }
  )
  if (!is.numeric(status) || status != 0L) {
    cli::cli_abort(
      "Beagle download returned non-zero status {.val {status}}.",
      class = "galick_error_beagle_download_failed"
    )
  }

  java <- find_java()
  detect_beagle_version(dest, java)

  invisible(dest)
}

#' Get the version of the resolved Beagle JAR
#'
#' Resolves the JAR the same way [run_beagle()] does (arg, then
#' `BEAGLE_PATH`, then `tools::R_user_dir("galick", "cache")`) and parses
#' the version banner.
#'
#' @param beagle_jar Optional explicit path.
#' @return A character scalar like `"5.4"`.
#' @export
beagle_version <- function(beagle_jar = NULL) {
  jar <- find_beagle_jar(beagle_jar)
  java <- find_java()
  detect_beagle_version(jar, java)$raw
}
