#' Resolve and fingerprint every input that determines a Beagle computation
#' @keywords internal
#' @noRd
prepare_beagle_job <- function(gt, out, ref = NULL, map = NULL,
                               mode = c("impute", "phase"), nthreads = 1L,
                               memory_gb = 4L, seed = NULL, beagle_jar = NULL,
                               java = NULL, extra_args = character()) {
  mode <- match.arg(mode)
  for (name in c("nthreads", "memory_gb")) {
    value <- get(name)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
          value < 1 || value != floor(value) || value > .Machine$integer.max) {
      cli::cli_abort("{.arg {name}} must be a positive integer.",
        class = "beagler_error_run_beagle_invalid_input"
      )
    }
  }
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L ||
                           !is.finite(seed) || seed != floor(seed) || abs(seed) > .Machine$integer.max)) {
    cli::cli_abort("{.arg seed} must be NULL or a finite integer.",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  if (!is.character(extra_args)) {
    cli::cli_abort("extra_args must be a character vector of key=value settings.",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  keys <- sub("=.*$", "", extra_args)
  reserved <- c("gt", "out", "ref", "map", "nthreads", "seed", "impute")
  if (!is.character(extra_args) || anyNA(extra_args) ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_-]*=[^\r\n]+$", extra_args)) ||
        any(keys %in% reserved) || anyDuplicated(keys)) {
    cli::cli_abort("extra_args must contain unique key=value settings, without overriding typed arguments.",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  java_path <- find_java(java)
  java_version <- detect_java_version(java_path)
  jar_path <- find_beagle_jar(beagle_jar)
  beagle_version <- detect_beagle_version(jar_path, java_path)
  files <- list(gt = gt, ref = ref, map = map, jar = jar_path, java = java_path)
  files <- files[!vapply(files, is.null, logical(1))]
  files <- lapply(files, function(path) {
    list(
      path = normalize_path_safe(path, must_exist = TRUE), sha256 = file_sha256(path)
    )
  })
  argv <- build_beagle_argv(jar_path,
    gt = files$gt$path, out = out,
    ref = if (!is.null(ref)) files$ref$path else NULL,
    map = if (!is.null(map)) files$map$path else NULL,
    mode = mode, nthreads = nthreads, memory_gb = memory_gb,
    seed = seed, extra_args = extra_args
  )
  versions <- list(
    beagler = as.character(utils::packageVersion("beagler")),
    r = as.character(getRversion()), os = R.version$os,
    beagle = beagle_version$raw, java = java_version$raw
  )
  fingerprint <- digest::digest(jsonlite::toJSON(
    list(contract = 1L, files = files, argv = argv, versions = versions),
    auto_unbox = TRUE, null = "null"
  ), algo = "sha256", serialize = FALSE)
  list(
    java = java_path, jar = jar_path, argv = argv, versions = versions,
    files = files, fingerprint = fingerprint
  )
}
