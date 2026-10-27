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
      class = "galick_error_split_vcf_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "galick_error_split_vcf_invalid_input"
    )
  }
  if (length(out_dir) != 1L || !is.character(out_dir) || is.na(out_dir)) {
    cli::cli_abort(
      "{.arg out_dir} must be a single non-NA character string.",
      class = "galick_error_split_vcf_invalid_input"
    )
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  parsed <- parse_vcf_raw(vcf)
  if (is.na(parsed$chrom_header_line)) {
    cli::cli_abort(
      "{.arg vcf} has no #CHROM header; cannot split.",
      class = "galick_error_split_vcf_invalid_input"
    )
  }

  header_lines <- c(parsed$meta, parsed$chrom_header_line)

  n_expected <- length(parsed$header_cols)
  body_chroms <- vapply(parsed$body_split, function(x) {
    if (length(x) >= 1L) x[1] else NA_character_
  }, character(1))
  body_col_ok <- parsed$body_col_counts == n_expected

  keepable <- !is.na(body_chroms) & body_col_ok
  if (any(!keepable)) {
    cli::cli_abort("Malformed VCF records cannot be silently omitted during splitting.",
      class = "galick_error_split_vcf_invalid_input"
    )
  }
  unique_chroms <- unique(body_chroms[keepable])
  if (!is.null(chrom_filter)) {
    unique_chroms <- intersect(unique_chroms, chrom_filter)
  }
  if (anyDuplicated(vapply(unique_chroms, safe_chrom_filename, character(1)))) {
    cli::cli_abort("Chromosome names collide after filename sanitization.",
      class = "galick_error_split_vcf_invalid_input"
    )
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
        class = "galick_error_split_vcf_output_exists"
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
#' and returns one `galick_run` per chromosome.
#'
#' Resume logic: if a per-chromosome `run_dir/<chrom>/run_record.json`
#' reports success, its input, reference, map, runtime and settings fingerprint
#' matches, and the output checksum and VCF validation succeed, the previous
#' result is returned. Older records without a fingerprint are rerun.
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
#' @return A named list of `galick_run` objects (names = chromosome).
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
  options <- list(...)
  if (length(options) && (is.null(names(options)) || any(!nzchar(names(options))) ||
                            anyDuplicated(names(options)) ||
                            any(names(options) %in% c("gt", "out", "ref", "map", "run_dir", "overwrite")))) {
    cli::cli_abort("Extra run options must be uniquely named, non-reserved arguments.",
      class = "galick_error_run_beagle_by_chrom_invalid_input"
    )
  }

  if (length(gt) != 1L || !is.character(gt) || is.na(gt) || !file.exists(gt)) {
    cli::cli_abort(
      "{.arg gt} must be a path to an existing VCF file.",
      class = "galick_error_run_beagle_by_chrom_invalid_input"
    )
  }
  if (length(out_dir) != 1L || !is.character(out_dir) || is.na(out_dir)) {
    cli::cli_abort(
      "{.arg out_dir} must be a single non-NA character string.",
      class = "galick_error_run_beagle_by_chrom_invalid_input"
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

    job <- do.call(
      prepare_beagle_job,
      c(list(gt = chrom_gt, out = out_prefix, ref = ref, map = map), options)
    )
    cached <- try_resume_chrom_run(
      chrom_gt, chrom_run_dir, out_prefix,
      job$fingerprint
    )
    if (!is.null(cached)) {
      return(cached)
    }

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
        class = "galick_error_run_beagle_by_chrom_invalid_input"
      )
    }
    future.apply::future_lapply(chroms, run_one)
  } else if (parallel == "mclapply") {
    parallel::mclapply(chroms, run_one, mc.cores = as.integer(workers))
  } else {
    lapply(chroms, run_one)
  }

  structure(results, names = chroms, expected_chromosomes = chroms)
}

#' Merge per-chromosome Beagle outputs into a single VCF
#'
#' Takes a list of `galick_run` objects (typically the output of
#' [run_beagle_by_chrom()]) and concatenates their VCF outputs into
#' `output_vcf`. All supplied runs must have succeeded. Sample columns are
#' aligned by name, compatible headers are combined, and FORMAT fields retained.
#'
#' Compressed output requires `bgzip` and is BGZF, not ordinary gzip.
#'
#' @param run_results A named list of `galick_run` objects.
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
      "{.arg run_results} must be a non-empty list of {.cls galick_run} objects.",
      class = "galick_error_merge_beagle_invalid_input"
    )
  }
  if (length(output_vcf) != 1L || !is.character(output_vcf) || is.na(output_vcf)) {
    cli::cli_abort(
      "{.arg output_vcf} must be a single non-NA character string.",
      class = "galick_error_merge_beagle_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", output_vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg output_vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "galick_error_merge_beagle_invalid_input"
    )
  }
  if (file.exists(output_vcf) && !overwrite) {
    cli::cli_abort(
      "{.arg output_vcf} exists: {.file {output_vcf}} (use {.code overwrite = TRUE}).",
      class = "galick_error_merge_beagle_output_exists"
    )
  }

  paths <- vapply(run_results, function(r) {
    if (!inherits(r, "galick_run") || !identical(r$status, 0L) ||
          length(r$output_paths) != 1L || is.na(r$output_paths) ||
          !file.exists(r$output_paths)) {
      cli::cli_abort(
        "Every expected run must have succeeded and have one existing VCF output.",
        class = "galick_error_merge_beagle_invalid_input"
      )
    }
    r$output_paths
  }, character(1), USE.NAMES = FALSE)
  for (path in paths) {
    checked <- validate_vcf(path)
    if (!checked$ok || checked$summary$n_records < 1L) {
      cli::cli_abort("Every chromosome output must contain valid VCF records.",
        class = "galick_error_merge_beagle_invalid_input"
      )
    }
  }
  inputs <- lapply(paths, read_genotype_vcf)
  expected <- attr(run_results, "expected_chromosomes")
  observed <- lapply(inputs, function(x) unique(x$fixed$CHROM))
  if (any(lengths(observed) != 1L) || anyDuplicated(unlist(observed)) ||
        (!is.null(expected) && !setequal(expected, unlist(observed)))) {
    cli::cli_abort("Chromosome outputs must cover the expected disjoint chromosomes.",
      class = "galick_error_merge_beagle_invalid_input"
    )
  }
  samples <- inputs[[1]]$samples
  if (any(!vapply(inputs, function(x) setequal(x$samples, samples), logical(1)))) {
    cli::cli_abort("Chromosome outputs must contain the same samples.",
      class = "galick_error_merge_beagle_invalid_input"
    )
  }
  for (i in seq_along(inputs)) {
    assert_resolvable_genotypes(inputs[[i]])
    index <- match(samples, inputs[[i]]$samples)
    inputs[[i]]$gt <- inputs[[i]]$gt[, index, drop = FALSE]
    inputs[[i]]$format_fields <- lapply(
      inputs[[i]]$format_fields,
      function(field) field[, index, drop = FALSE]
    )
    inputs[[i]]$samples <- samples
  }
  meta <- merge_vcf_metadata(inputs)
  inputs <- namespace_phase_sets(inputs, preserve_fields = TRUE)
  fields <- unique(unlist(lapply(inputs, function(x) names(x$format_fields))))
  values <- lapply(fields, function(field) {
    do.call(rbind, lapply(inputs, function(x) {
      if (is.null(x$format_fields[[field]])) {
        matrix(".", nrow(x$gt), ncol(x$gt))
      } else {
        x$format_fields[[field]]
      }
    }))
  })
  names(values) <- fields
  combined <- new_genotype_vcf(
    data.table::rbindlist(lapply(inputs, `[[`, "fixed")),
    do.call(rbind, lapply(inputs, `[[`, "gt")), samples, meta,
    format_fields = values
  )
  combined$formats <- unlist(lapply(inputs, `[[`, "formats"), use.names = FALSE)
  if ("PS" %in% fields) {
    combined$meta <- combined$meta[!grepl("^##FORMAT=<ID=PS[,>]", combined$meta)]
    combined$meta <- c(
      combined$meta,
      '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Input-specific phase set">'
    )
  }
  write_genotype_vcf(combined, output_vcf, overwrite)
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
#' input hash, return a synthesized `galick_run` for it; otherwise NULL.
#'
#' @keywords internal
#' @noRd
try_resume_chrom_run <- function(chrom_gt, run_dir, out_prefix,
                                 expected_fingerprint = NULL) {
  tryCatch(
    {
      rec_path <- file.path(run_dir, "run_record.json")
      out_path <- paste0(out_prefix, ".vcf.gz")
      if (!file.exists(rec_path) || !file.exists(out_path)) {
        return(NULL)
      }
      rec <- tryCatch(
        jsonlite::fromJSON(rec_path, simplifyVector = FALSE),
        error = function(e) NULL
      )
      if (is.null(rec) || !identical(rec$status, 0L) && !identical(rec$status, 0)) {
        return(NULL)
      }
      rec$command <- as.character(unlist(rec$command, use.names = FALSE))
      validate_run_record(rec)
      if (is.null(expected_fingerprint) ||
            !identical(rec$args$job_fingerprint, expected_fingerprint)) {
        return(NULL)
      }
      expected_sha <- file_sha256(chrom_gt)
      recorded_sha <- if (length(rec$inputs) > 0L) rec$inputs[[1]]$sha256 else NULL
      if (is.null(recorded_sha) || recorded_sha != expected_sha) {
        return(NULL)
      }
      if (length(rec$outputs) != 1L ||
        !identical(
          normalize_path_safe(rec$outputs[[1]]$path, must_exist = TRUE),
          normalize_path_safe(out_path, must_exist = TRUE)
        ) ||
        !identical(rec$outputs[[1]]$sha256, file_sha256(out_path))) {
        return(NULL)
      }
      valid <- tryCatch(validate_vcf(out_path, target = "beagle"), error = function(e) NULL)
      if (is.null(valid) || !valid$ok || valid$summary$n_records < 1L) {
        return(NULL)
      }

      versions <- if (is.null(rec$versions)) list() else rec$versions
      command <- if (is.null(rec$command)) character() else as.character(unlist(rec$command))
      new_galick_run(
        status = 0L,
        runtime_sec = if (is.null(rec$runtime_sec)) NA_real_ else as.numeric(rec$runtime_sec),
        command = command,
        version = versions,
        output_paths = out_path,
        stdout_log = file.path(run_dir, "stdout.log"),
        stderr_log = file.path(run_dir, "stderr.log"),
        run_record = rec_path
      )
    },
    error = function(e) NULL
  )
}
