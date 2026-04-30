#' Split a VCF into per-chromosome files
#'
#' Reads `vcf` via the tolerant parser and writes one `.vcf` per
#' unique CHROM into `out_dir`. Pure R; no `bcftools` dependency.
#'
#' @param vcf Path to a `.vcf` or `.vcf.gz` file.
#' @param out_dir Directory to write per-chromosome files into. Created if
#'   it doesn't exist.
#' @param chrom_filter Optional character vector of chromosome labels to
#'   keep; everything else is dropped from the output.
#' @param bcftools Optional path to `bcftools`. Currently unused (R-only
#'   path always taken in v0.1.0); kept in the signature for forward
#'   compatibility with a fast path.
#' @param overwrite Refuse to clobber per-chromosome files unless `TRUE`.
#'
#' @return A list with `output_dir` and a named character vector
#'   `chrom_paths` (names are CHROM labels, values are file paths).
#' @export
split_vcf_by_chrom <- function(vcf, out_dir,
                               chrom_filter = NULL,
                               bcftools = NULL,
                               overwrite = FALSE) {
  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf) || !file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a path to an existing VCF file.",
      class = "beagler_error_split_vcf_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_split_vcf_invalid_input"
    )
  }
  if (length(out_dir) != 1L || !is.character(out_dir) || is.na(out_dir)) {
    cli::cli_abort(
      "{.arg out_dir} must be a single non-NA character string.",
      class = "beagler_error_split_vcf_invalid_input"
    )
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  parsed <- parse_vcf_raw(vcf)
  if (is.na(parsed$chrom_header_line)) {
    cli::cli_abort(
      "{.arg vcf} has no #CHROM header; cannot split.",
      class = "beagler_error_split_vcf_invalid_input"
    )
  }

  header_lines <- c(parsed$meta, parsed$chrom_header_line)

  n_expected <- length(parsed$header_cols)
  body_chroms <- vapply(parsed$body_split, function(x) {
    if (length(x) >= 1L) x[1] else NA_character_
  }, character(1))
  body_col_ok <- parsed$body_col_counts == n_expected

  keepable <- !is.na(body_chroms) & body_col_ok
  unique_chroms <- unique(body_chroms[keepable])
  if (!is.null(chrom_filter)) {
    unique_chroms <- intersect(unique_chroms, chrom_filter)
  }

  paths <- character(length(unique_chroms))
  names(paths) <- unique_chroms
  for (ch in unique_chroms) {
    out_path <- file.path(out_dir, sprintf("%s.vcf", safe_chrom_filename(ch)))
    if (file.exists(out_path) && !overwrite) {
      cli::cli_abort(
        c(
          "Per-chromosome output already exists: {.file {out_path}}",
          i = "Use {.code overwrite = TRUE} to clobber."
        ),
        class = "beagler_error_split_vcf_output_exists"
      )
    }
    keep_idx <- which(keepable & body_chroms == ch)
    writeLines(c(header_lines, parsed$body_lines[keep_idx]), out_path)
    paths[ch] <- out_path
  }

  list(output_dir = out_dir, chrom_paths = paths)
}

#' Run Beagle separately on each chromosome
#'
#' Splits `gt` into per-chromosome files, calls [run_beagle()] on each,
#' and returns one `beagler_run` per chromosome.
#'
#' Resume logic: if a per-chromosome `run_dir/<chrom>/run_record.json`
#' reports `status == 0` and the recorded input SHA-256 matches the
#' current per-chromosome `gt`, the chromosome is skipped and the
#' previous result is returned.
#'
#' @param gt Path to the target VCF.
#' @param out_dir Top-level output directory; per-chromosome subdirs are
#'   created inside it.
#' @param ref,map Passed verbatim to [run_beagle()] for each chromosome.
#' @param chrom_filter Optional character vector to restrict which
#'   chromosomes are run.
#' @param parallel `"none"` (default), `"future"`, or `"mclapply"`. When
#'   not `"none"`, `workers` controls parallelism.
#' @param workers Number of workers for parallel execution.
#' @param ... Extra args passed verbatim to each [run_beagle()] call
#'   (e.g. `mode`, `nthreads`, `memory_gb`).
#'
#' @return A named list of `beagler_run` objects (names = chromosome).
#' @export
run_beagle_by_chrom <- function(gt,
                                out_dir,
                                ref = NULL,
                                map = NULL,
                                chrom_filter = NULL,
                                parallel = c("none", "future", "mclapply"),
                                workers = 1L,
                                ...) {
  parallel <- match.arg(parallel)

  if (length(gt) != 1L || !is.character(gt) || is.na(gt) || !file.exists(gt)) {
    cli::cli_abort(
      "{.arg gt} must be a path to an existing VCF file.",
      class = "beagler_error_run_beagle_by_chrom_invalid_input"
    )
  }
  if (length(out_dir) != 1L || !is.character(out_dir) || is.na(out_dir)) {
    cli::cli_abort(
      "{.arg out_dir} must be a single non-NA character string.",
      class = "beagler_error_run_beagle_by_chrom_invalid_input"
    )
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  splits <- split_vcf_by_chrom(
    gt,
    file.path(out_dir, "_split"),
    chrom_filter = chrom_filter,
    overwrite = TRUE
  )

  chroms <- names(splits$chrom_paths)
  if (length(chroms) == 0L) {
    return(structure(list(), names = character()))
  }

  run_one <- function(ch) {
    chrom_gt <- splits$chrom_paths[[ch]]
    chrom_run_dir <- file.path(out_dir, safe_chrom_filename(ch))
    if (!dir.exists(chrom_run_dir)) {
      dir.create(chrom_run_dir, recursive = TRUE, showWarnings = FALSE)
    }
    out_prefix <- file.path(chrom_run_dir, "imputed")

    cached <- try_resume_chrom_run(chrom_gt, chrom_run_dir, out_prefix)
    if (!is.null(cached)) return(cached)

    run_beagle(
      gt = chrom_gt,
      out = out_prefix,
      ref = ref,
      map = map,
      run_dir = chrom_run_dir,
      overwrite = TRUE,
      ...
    )
  }

  results <- if (parallel == "future") {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      cli::cli_abort(
        "{.pkg future.apply} is required for {.arg parallel = \"future\"}.",
        class = "beagler_error_run_beagle_by_chrom_invalid_input"
      )
    }
    future.apply::future_lapply(chroms, run_one)
  } else if (parallel == "mclapply") {
    parallel::mclapply(chroms, run_one, mc.cores = as.integer(workers))
  } else {
    lapply(chroms, run_one)
  }

  stats::setNames(results, chroms)
}

#' Merge per-chromosome Beagle outputs into a single VCF
#'
#' Takes a list of `beagler_run` objects (typically the output of
#' [run_beagle_by_chrom()]) and concatenates their VCF outputs into
#' `output_vcf`. Header is taken from the first non-empty input;
#' subsequent inputs contribute body rows only.
#'
#' Pure R (uses `gzfile()` for `.vcf.gz` IO). A `bcftools concat`
#' fast path can be added in v0.2.x.
#'
#' @param run_results A named list of `beagler_run` objects.
#' @param output_vcf Path to the merged VCF (`.vcf` or `.vcf.gz`).
#' @param bcftools Optional path; ignored in v0.1.0 (kept for forward
#'   compatibility).
#' @param overwrite Refuse to clobber `output_vcf` unless `TRUE`.
#'
#' @return The output path, invisibly.
#' @export
merge_beagle_outputs <- function(run_results, output_vcf,
                                 bcftools = NULL, overwrite = FALSE) {
  if (!is.list(run_results) || length(run_results) == 0L) {
    cli::cli_abort(
      "{.arg run_results} must be a non-empty list of {.cls beagler_run} objects.",
      class = "beagler_error_merge_beagle_invalid_input"
    )
  }
  if (length(output_vcf) != 1L || !is.character(output_vcf) || is.na(output_vcf)) {
    cli::cli_abort(
      "{.arg output_vcf} must be a single non-NA character string.",
      class = "beagler_error_merge_beagle_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", output_vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg output_vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_merge_beagle_invalid_input"
    )
  }
  if (file.exists(output_vcf) && !overwrite) {
    cli::cli_abort(
      "{.arg output_vcf} exists: {.file {output_vcf}} (use {.code overwrite = TRUE}).",
      class = "beagler_error_merge_beagle_output_exists"
    )
  }

  paths <- unlist(lapply(run_results, function(r) {
    if (!inherits(r, "beagler_run")) {
      cli::cli_abort(
        "Each element of {.arg run_results} must be a {.cls beagler_run}.",
        class = "beagler_error_merge_beagle_invalid_input"
      )
    }
    r$output_paths
  }), use.names = FALSE)

  paths <- paths[file.exists(paths)]
  if (length(paths) == 0L) {
    cli::cli_abort(
      "None of the {.cls beagler_run} objects has an existing output file.",
      class = "beagler_error_merge_beagle_invalid_input"
    )
  }

  first_lines <- read_vcf_lines(paths[1])
  output_lines <- first_lines
  if (length(paths) > 1L) {
    for (p in paths[-1]) {
      lines <- read_vcf_lines(p)
      body_only <- lines[!startsWith(lines, "#") & nzchar(lines)]
      output_lines <- c(output_lines, body_only)
    }
  }

  if (grepl("\\.vcf\\.gz$", output_vcf, ignore.case = TRUE)) {
    con <- gzfile(output_vcf, "wt")
    on.exit(close(con), add = TRUE)
    writeLines(output_lines, con)
  } else {
    writeLines(output_lines, output_vcf)
  }
  invisible(output_vcf)
}

# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
safe_chrom_filename <- function(ch) {
  gsub("[^A-Za-z0-9_.-]", "_", ch)
}

#' If a previous successful run exists for this chromosome with a matching
#' input hash, return a synthesized `beagler_run` for it; otherwise NULL.
#'
#' @keywords internal
#' @noRd
try_resume_chrom_run <- function(chrom_gt, run_dir, out_prefix) {
  rec_path <- file.path(run_dir, "run_record.json")
  out_path <- paste0(out_prefix, ".vcf.gz")
  if (!file.exists(rec_path) || !file.exists(out_path)) {
    return(NULL)
  }
  rec <- tryCatch(
    jsonlite::fromJSON(rec_path, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (is.null(rec) || is.null(rec$status) || rec$status != 0L) {
    return(NULL)
  }
  expected_sha <- unname(tools::sha256sum(chrom_gt))
  recorded_sha <- if (length(rec$inputs) > 0L) rec$inputs[[1]]$sha256 else NULL
  if (is.null(recorded_sha) || recorded_sha != expected_sha) {
    return(NULL)
  }

  versions <- if (is.null(rec$versions)) list() else rec$versions
  command <- if (is.null(rec$command)) character() else as.character(unlist(rec$command))
  new_beagler_run(
    status = 0L,
    runtime_sec = if (is.null(rec$runtime_sec)) NA_real_ else as.numeric(rec$runtime_sec),
    command = command,
    version = versions,
    output_paths = out_path,
    stdout_log = file.path(run_dir, "stdout.log"),
    stderr_log = file.path(run_dir, "stderr.log"),
    run_record = rec_path
  )
}
