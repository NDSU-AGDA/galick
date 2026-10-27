#' Tool discovery (internal)
#' @name find-tools
#' @keywords internal
NULL

#' Locate the Java executable
#'
#' Resolution order: explicit `java` arg → `JAVA_HOME/bin/java` →
#' `Sys.which("java")`. Aborts with `galick_error_java_not_found` if
#' nothing usable is found.
#'
#' @keywords internal
#' @noRd
find_java <- function(java = NULL) {
  if (!is.null(java)) {
    if (length(java) != 1L || !is.character(java) || is.na(java)) {
      cli::cli_abort(
        "{.arg java} must be a single non-NA character string.",
        class = "galick_error_java_not_found"
      )
    }
    if (!file.exists(java)) {
      cli::cli_abort(
        "Java path does not exist: {.file {java}}",
        class = "galick_error_java_not_found"
      )
    }
    return(normalize_path_safe(java, must_exist = TRUE))
  }

  java_home <- Sys.getenv("JAVA_HOME", unset = "")
  if (nzchar(java_home)) {
    candidate <- file.path(java_home, "bin", "java")
    if (.Platform$OS.type == "windows") {
      candidate <- paste0(candidate, ".exe")
    }
    if (file.exists(candidate)) {
      return(normalize_path_safe(candidate, must_exist = TRUE))
    }
  }

  via_path <- Sys.which("java")
  if (nzchar(via_path)) {
    return(unname(via_path))
  }

  cli::cli_abort(
    c(
      "Could not find {.code java}.",
      i = "Install Java 8+ or set {.envvar JAVA_HOME}.",
      i = "macOS: {.code brew install openjdk}; Ubuntu: {.code apt install default-jre}."
    ),
    class = "galick_error_java_not_found"
  )
}

#' Detect Java version (and refuse < 8)
#' @keywords internal
#' @noRd
detect_java_version <- function(java) {
  res <- run_cli(java, c("-version"))
  parsed <- parse_java_version_output(res$stderr)
  if (is.null(parsed)) {
    cli::cli_abort(
      "Could not parse Java version from {.code {java} -version}.",
      class = "galick_error_java_version_unknown"
    )
  }
  if (parsed$major < 8L) {
    cli::cli_abort(
      "Java 8 or later required; found {.val {parsed$raw}}.",
      class = "galick_error_java_version_too_old"
    )
  }
  parsed
}

#' Parse `java -version` stderr output
#'
#' Java prints the version line to stderr in one of two forms:
#' - `openjdk version "11.0.x" 2021-...` (Java 9+)
#' - `java version "1.8.0_xxx"` (Java 8 and earlier)
#'
#' Returns `list(major, minor, raw)` or `NULL` on parse failure.
#'
#' @keywords internal
#' @noRd
parse_java_version_output <- function(text) {
  joined <- paste(text, collapse = "\n")
  m <- regmatches(joined, regexpr('"[0-9._]+"', joined))
  if (length(m) == 0L) {
    return(NULL)
  }
  raw <- gsub('"', "", m[1])
  parts <- strsplit(raw, "[._]")[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    return(NULL)
  }
  if (parts[1] == "1" && length(parts) >= 2L) {
    major <- suppressWarnings(as.integer(parts[2]))
  } else {
    major <- suppressWarnings(as.integer(parts[1]))
  }
  if (is.na(major)) {
    return(NULL)
  }
  minor <- if (length(parts) >= 2L) suppressWarnings(as.integer(parts[2])) else 0L
  if (is.na(minor)) minor <- 0L
  list(major = major, minor = minor, raw = raw)
}

#' Locate the Beagle JAR
#'
#' Resolution order: explicit `beagle_jar` arg → `BEAGLE_PATH` env →
#' `tools::R_user_dir("galick", "cache")/beagle.jar`. Aborts with
#' `galick_error_beagle_not_found` and points at `install_beagle()`
#' if nothing is found.
#'
#' @keywords internal
#' @noRd
find_beagle_jar <- function(beagle_jar = NULL) {
  if (!is.null(beagle_jar)) {
    if (length(beagle_jar) != 1L || !is.character(beagle_jar) || is.na(beagle_jar)) {
      cli::cli_abort(
        "{.arg beagle_jar} must be a single non-NA character string.",
        class = "galick_error_beagle_not_found"
      )
    }
    if (!file.exists(beagle_jar)) {
      cli::cli_abort(
        "Beagle JAR does not exist: {.file {beagle_jar}}",
        class = "galick_error_beagle_not_found"
      )
    }
    return(normalize_path_safe(beagle_jar, must_exist = TRUE))
  }

  env_path <- Sys.getenv("BEAGLE_PATH", unset = "")
  if (nzchar(env_path) && file.exists(env_path)) {
    return(normalize_path_safe(env_path, must_exist = TRUE))
  }

  cache_path <- file.path(tools::R_user_dir("galick", "cache"), "beagle.jar")
  if (file.exists(cache_path)) {
    return(normalize_path_safe(cache_path, must_exist = TRUE))
  }

  cli::cli_abort(
    c(
      "Could not find a Beagle JAR.",
      i = "Pass {.arg beagle_jar} explicitly,",
      i = "set {.envvar BEAGLE_PATH} to its location,",
      i = "or run {.code galick::install_beagle()} to download it."
    ),
    class = "galick_error_beagle_not_found"
  )
}

#' Detect Beagle version (and refuse < 5)
#' @keywords internal
#' @noRd
detect_beagle_version <- function(beagle_jar, java) {
  res <- run_cli(java, c("-jar", beagle_jar))
  parsed <- parse_beagle_version_output(c(res$stdout, res$stderr))
  if (is.null(parsed)) {
    cli::cli_abort(
      "Could not parse Beagle version from {.file {beagle_jar}}.",
      class = "galick_error_beagle_version_unknown"
    )
  }
  if (parsed$major < 5L) {
    cli::cli_abort(
      "Beagle 5.x or later required; found {.val {parsed$raw}}.",
      class = "galick_error_beagle_version_too_old"
    )
  }
  parsed
}

#' Parse Beagle's no-args banner for the `(version X.Y)` line
#'
#' Beagle 5.x prints something like:
#'   `beagle.<date>.<hash>.jar  (version 5.4)`
#'
#' Returns `list(major, minor, raw)` or `NULL` on parse failure.
#'
#' @keywords internal
#' @noRd
parse_beagle_version_output <- function(text) {
  joined <- paste(text, collapse = "\n")
  m <- regmatches(joined, regexpr("\\(version[ \t]+[0-9.]+\\)", joined))
  if (length(m) == 0L) {
    return(NULL)
  }
  raw <- gsub("[()]", "", m[1])
  raw <- sub("^version[ \t]+", "", raw)
  parts <- strsplit(raw, ".", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0L) {
    return(NULL)
  }
  major <- suppressWarnings(as.integer(parts[1]))
  if (is.na(major)) {
    return(NULL)
  }
  minor <- if (length(parts) >= 2L) suppressWarnings(as.integer(parts[2])) else 0L
  if (is.na(minor)) minor <- 0L
  list(major = major, minor = minor, raw = raw)
}

#' Locate the PLINK executable
#'
#' Resolution order: explicit `plink` arg -> `Sys.which("plink")`.
#' Aborts with `galick_error_plink_not_found` if nothing usable.
#'
#' @keywords internal
#' @noRd
find_plink <- function(plink = NULL) {
  if (!is.null(plink)) {
    if (length(plink) != 1L || !is.character(plink) || is.na(plink)) {
      cli::cli_abort(
        "{.arg plink} must be a single non-NA character string.",
        class = "galick_error_plink_not_found"
      )
    }
    if (!file.exists(plink)) {
      cli::cli_abort(
        "PLINK path does not exist: {.file {plink}}",
        class = "galick_error_plink_not_found"
      )
    }
    return(normalize_path_safe(plink, must_exist = TRUE))
  }
  via <- Sys.which("plink")
  if (nzchar(via)) {
    return(unname(via))
  }
  cli::cli_abort(
    c(
      "Could not find {.code plink}.",
      i = "Install PLINK 1.9 or 2.0 and ensure it is on PATH.",
      i = "macOS: {.code brew install plink}; Ubuntu: {.code apt install plink1.9}."
    ),
    class = "galick_error_plink_not_found"
  )
}

#' @keywords internal
#' @noRd
find_bcftools <- function() {
  via <- Sys.which("bcftools")
  if (nzchar(via)) unname(via) else NULL
}

#' @keywords internal
#' @noRd
find_tabix <- function() {
  via <- Sys.which("tabix")
  if (nzchar(via)) unname(via) else NULL
}
