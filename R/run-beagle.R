#' Run Beagle for phasing or imputation
#'
#' Wraps Beagle 5.x with safe defaults, structured argv (no shell strings),
#' and a captured run record. Returns a `beagler_run` even when Beagle
#' itself exits non-zero — callers can inspect stderr via `print()`.
#'
#' Pre-flight: `validate_vcf(gt, target = "beagle")` must succeed. Pre-flight
#' failure raises `beagler_error_run_beagle_invalid_input` (repair the input
#' with [repair_vcf()] first). Post-flight failures (Beagle exits non-zero)
#' come back as a `beagler_run` with `status != 0`.
#'
#' v0.1.0 limitations:
#' - `gt` must be plain `.vcf`. v0.1.0's [validate_vcf()] does not decompress;
#'   `gunzip` your file first if needed.
#' - `ref` must be `.vcf.gz` (BGZF) with a sibling `.tbi` index. Beagle
#'   requires tabix-indexed reference panels.
#' - `out` is a Beagle output **prefix**, not a filename. Passing `.vcf` or
#'   `.vcf.gz` is rejected. Beagle writes `<out>.vcf.gz`.
#'
#' @param gt Path to plain `.vcf` target genotype file.
#' @param out Output prefix (no extension).
#' @param ref Optional `.vcf.gz` reference panel (must have `.tbi` index).
#' @param map Optional PLINK-format genetic map (`<chrom>.map` per Beagle).
#' @param mode `"impute"` (default) sets `impute=true` in Beagle args;
#'   `"phase"` sets `impute=false`.
#' @param nthreads Number of threads.
#' @param memory_gb JVM `-Xmx` heap size in GB.
#' @param seed Optional integer random seed.
#' @param beagle_jar Optional path; auto-discovered if `NULL`.
#' @param java Optional path; auto-discovered if `NULL`.
#' @param extra_args Character vector of extra Beagle arguments, passed
#'   verbatim *after* the package's defaults (so the user can override them).
#'   No validation; the caller is responsible.
#' @param run_dir Directory for `stdout.log`, `stderr.log`, `run_record.json`.
#'   Auto-created if `NULL`.
#' @param overwrite Refuse to clobber `<out>.vcf.gz` unless `TRUE`.
#'
#' @return A `beagler_run` S3 object.
#' @export
run_beagle <- function(gt,
                      out,
                      ref = NULL,
                      map = NULL,
                      mode = c("impute", "phase"),
                      nthreads = 1L,
                      memory_gb = 4L,
                      seed = NULL,
                      beagle_jar = NULL,
                      java = NULL,
                      extra_args = character(),
                      run_dir = NULL,
                      overwrite = FALSE) {
  the_call <- match.call()
  mode <- match.arg(mode)

  if (length(gt) != 1L || !is.character(gt) || is.na(gt)) {
    cli::cli_abort(
      "{.arg gt} must be a single non-NA character string.",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  if (!file.exists(gt)) {
    cli::cli_abort(
      "{.arg gt} not found: {.file {gt}}",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", gt, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg gt} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }

  if (length(out) != 1L || !is.character(out) || is.na(out)) {
    cli::cli_abort(
      "{.arg out} must be a single non-NA character string.",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  if (grepl("\\.vcf$|\\.vcf\\.gz$", out, ignore.case = TRUE)) {
    cli::cli_abort(
      c(
        "{.arg out} is a Beagle output prefix; it must not end in {.val .vcf} or {.val .vcf.gz}.",
        i = "Beagle appends {.val .vcf.gz} to whatever you pass."
      ),
      class = "beagler_error_run_beagle_invalid_input"
    )
  }
  expected_output <- paste0(out, ".vcf.gz")
  if (file.exists(expected_output) && !overwrite) {
    cli::cli_abort(
      "Output already exists: {.file {expected_output}} (use {.code overwrite = TRUE}).",
      class = "beagler_error_run_beagle_output_exists"
    )
  }

  if (!is.null(ref)) {
    if (length(ref) != 1L || !is.character(ref) || is.na(ref) || !file.exists(ref)) {
      cli::cli_abort(
        "{.arg ref} must be a path to an existing VCF file.",
        class = "beagler_error_run_beagle_invalid_input"
      )
    }
    if (!grepl("\\.vcf\\.gz$", ref, ignore.case = TRUE)) {
      cli::cli_abort(
        c(
          "{.arg ref} must be a {.val .vcf.gz} (BGZF-compressed) file.",
          i = "Beagle requires BGZF + tabix-indexed reference panels.",
          i = "Run {.code bgzip ref.vcf && tabix -p vcf ref.vcf.gz} first."
        ),
        class = "beagler_error_run_beagle_invalid_input"
      )
    }
    if (!file.exists(paste0(ref, ".tbi"))) {
      cli::cli_abort(
        c(
          "{.arg ref} is not tabix-indexed: {.file {paste0(ref, '.tbi')}} not found.",
          i = "Run {.code tabix -p vcf {ref}} to create it."
        ),
        class = "beagler_error_run_beagle_invalid_input"
      )
    }
  }

  if (!is.null(map) && !file.exists(map)) {
    cli::cli_abort(
      "{.arg map} not found: {.file {map}}",
      class = "beagler_error_run_beagle_invalid_input"
    )
  }

  java_path <- find_java(java)
  java_version <- detect_java_version(java_path)

  beagle_path <- find_beagle_jar(beagle_jar)
  beagle_ver <- detect_beagle_version(beagle_path, java_path)

  gt_validation <- validate_vcf(gt, target = "beagle")
  blocking <- gt_validation$errors$code
  if (length(blocking) > 0L) {
    cli::cli_abort(
      c(
        "{.arg gt} failed pre-flight validation ({length(blocking)} error(s)).",
        i = "Repair with {.code repair_vcf()} before calling {.fn run_beagle}.",
        i = "Codes: {.val {blocking}}"
      ),
      class = "beagler_error_run_beagle_invalid_input"
    )
  }

  if (is.null(run_dir)) {
    run_dir <- tempfile("beagler_run_")
  }
  if (!dir.exists(run_dir)) {
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  }
  run_dir <- normalize_path_safe(run_dir, must_exist = TRUE)

  args <- build_beagle_argv(
    beagle_jar = beagle_path,
    gt = gt,
    out = out,
    ref = ref,
    map = map,
    mode = mode,
    nthreads = nthreads,
    memory_gb = memory_gb,
    seed = seed,
    extra_args = extra_args
  )

  stdout_log <- file.path(run_dir, "stdout.log")
  stderr_log <- file.path(run_dir, "stderr.log")
  started_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")

  result <- run_cli(
    java_path,
    args = args,
    stdout_log = stdout_log,
    stderr_log = stderr_log
  )

  ended_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")

  output_paths <- if (file.exists(expected_output)) expected_output else character()

  index_action <- attempt_index(expected_output, output_paths)

  inputs <- list(list(path = gt, sha256 = unname(tools::sha256sum(gt))))
  if (!is.null(ref)) {
    inputs[[length(inputs) + 1L]] <- list(
      path = ref, sha256 = unname(tools::sha256sum(ref))
    )
  }
  if (!is.null(map)) {
    inputs[[length(inputs) + 1L]] <- list(
      path = map, sha256 = unname(tools::sha256sum(map))
    )
  }
  outputs <- lapply(output_paths, function(p) {
    list(path = p, sha256 = unname(tools::sha256sum(p)))
  })

  versions <- compact_list(list(
    beagler = as.character(utils::packageVersion("beagler")),
    r = as.character(getRversion()),
    os = R.version$os,
    beagle = beagle_ver$raw,
    java = java_version$raw,
    bcftools = if (!is.null(find_bcftools())) "available" else NULL,
    tabix = if (!is.null(find_tabix())) "available" else NULL
  ))

  args_record <- compact_list(list(
    gt = gt,
    out = out,
    ref = ref,
    map = map,
    mode = mode,
    nthreads = as.integer(nthreads),
    memory_gb = as.integer(memory_gb),
    seed = if (!is.null(seed)) as.integer(seed) else NULL,
    extra_args = as.list(extra_args),
    index_action = index_action
  ))

  run_record <- write_run_record_json(
    run_dir = run_dir,
    command = c(java_path, args),
    versions = versions,
    args = args_record,
    inputs = inputs,
    outputs = outputs,
    runtime_sec = result$runtime_sec,
    status = result$status,
    started_at = started_at,
    ended_at = ended_at
  )

  new_beagler_run(
    status = result$status,
    runtime_sec = result$runtime_sec,
    command = c(java_path, args),
    version = versions,
    output_paths = output_paths,
    stdout_log = stdout_log,
    stderr_log = stderr_log,
    run_record = run_record,
    call = the_call
  )
}

# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

#' Construct the argv vector for `java -jar beagle.jar ...`
#'
#' Always emits, in order: JVM heap, `-jar`, jar path, then Beagle's
#' `key=value` arguments. `extra_args` is appended verbatim *after* package
#' defaults so users can override.
#'
#' @keywords internal
#' @noRd
build_beagle_argv <- function(beagle_jar, gt, out,
                              ref = NULL, map = NULL,
                              mode = "impute",
                              nthreads = 1L,
                              memory_gb = 4L,
                              seed = NULL,
                              extra_args = character()) {
  args <- c(
    sprintf("-Xmx%dg", as.integer(memory_gb)),
    "-jar",
    beagle_jar,
    sprintf("gt=%s", gt),
    sprintf("out=%s", out),
    sprintf("nthreads=%d", as.integer(nthreads))
  )
  if (!is.null(ref)) args <- c(args, sprintf("ref=%s", ref))
  if (!is.null(map)) args <- c(args, sprintf("map=%s", map))
  if (!is.null(seed)) args <- c(args, sprintf("seed=%d", as.integer(seed)))
  if (mode == "phase") args <- c(args, "impute=false")
  c(args, as.character(extra_args))
}

#' Best-effort tabix indexing of Beagle output
#'
#' Tries `tabix -p vcf` first, then `bcftools index -t`. Returns the action
#' name (`"tabix"`, `"bcftools"`) or `"skipped"` if neither tool is found
#' or if Beagle did not write the expected output.
#'
#' @keywords internal
#' @noRd
attempt_index <- function(expected_output, output_paths) {
  if (length(output_paths) == 0L || !file.exists(expected_output)) {
    return("skipped")
  }
  tabix <- find_tabix()
  if (!is.null(tabix)) {
    res <- run_cli(tabix, c("-p", "vcf", expected_output))
    if (res$status == 0L) return("tabix")
  }
  bcftools <- find_bcftools()
  if (!is.null(bcftools)) {
    res <- run_cli(bcftools, c("index", "-t", expected_output))
    if (res$status == 0L) return("bcftools")
  }
  "skipped"
}

#' @keywords internal
#' @noRd
compact_list <- function(x) x[!vapply(x, is.null, logical(1))]
