#' Repair a VCF file
#'
#' Reads a (possibly broken) plain `.vcf` via the tolerant internal
#' `parse_vcf_raw()` and applies a sequence of repair steps. The output is a
#' minimal VCF v4.2 file containing only `CHROM`, `POS`, `ID`, `REF`, `ALT`
#' and `FORMAT=GT` plus the per-sample GT subfield extracted from the input
#' (any other FORMAT subfields are dropped).
#'
#' File-level problems that cannot be safely repaired (missing `#CHROM`
#' header, no sample columns, duplicate sample names) are reported as
#' `$errors` and **no VCF is written**. Row-level problems are addressed by
#' the requested `steps`. After all steps, the repaired output is itself
#' validated under `target = "beagle"`; if any errors remain (e.g. because
#' the user excluded a step that would have addressed them), the result
#' carries `repair_vcf_remaining_*` errors and no VCF is written.
#'
#' Repair actions surface as `$warnings` with per-step counts in
#' `$summary$step_summary`.
#'
#' @param vcf Path to a plain `.vcf` file.
#' @param output_vcf Path ending in `.vcf` to write the repaired file.
#' @param steps Character vector naming the steps to run, in order.
#'   Default runs all five. Choices: `"drop_invalid_records"`,
#'   `"filter_biallelic"`, `"fix_marker_ids"`, `"remove_duplicates"`,
#'   `"sort"`. The user fully controls the set; nothing is run secretly.
#' @param chrom_style `"auto"` leaves `CHROM` labels untouched; the other
#'   values rewrite labels in the output (and use the normalized form for
#'   duplicate detection and sorting).
#' @param overwrite Refuse to clobber `output_vcf` unless `TRUE`.
#'
#' @return A `galick_result` S3 object. On success,
#'   `summary$output_path` is the written file. Otherwise it is
#'   `NA_character_`.
#'
#' @examples
#' dirty <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
#' output <- tempfile(fileext = ".vcf")
#' result <- repair_vcf(dirty, output)
#' result$ok
#' result$summary$step_summary
#' unlink(output)
#' @export
repair_vcf <- function(vcf,
                       output_vcf,
                       steps = c(
                         "drop_invalid_records", "filter_biallelic",
                         "fix_marker_ids", "remove_duplicates", "sort"
                       ),
                       chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
                       overwrite = FALSE) {
  the_call <- match.call()
  chrom_style <- match.arg(chrom_style)

  valid_steps <- c(
    "drop_invalid_records", "filter_biallelic",
    "fix_marker_ids", "remove_duplicates", "sort"
  )
  bad_steps <- setdiff(steps, valid_steps)
  if (length(bad_steps) > 0L) {
    cli::cli_abort(
      "Unknown {.arg steps}: {.val {bad_steps}}",
      class = "galick_error_repair_vcf_invalid_input"
    )
  }

  if (length(output_vcf) != 1L || !is.character(output_vcf) || is.na(output_vcf)) {
    cli::cli_abort(
      "{.arg output_vcf} must be a single non-NA character string.",
      class = "galick_error_repair_vcf_invalid_input"
    )
  }
  if (!grepl("\\.vcf$", output_vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "Output path must end in {.val .vcf}; got {.file {output_vcf}}.",
      class = "galick_error_repair_vcf_unsupported_extension"
    )
  }
  if (file.exists(output_vcf) && !overwrite) {
    cli::cli_abort(
      "Output file exists: {.file {output_vcf}} (use {.code overwrite = TRUE}).",
      class = "galick_error_repair_vcf_output_exists"
    )
  }

  parsed <- parse_vcf_raw(vcf)

  blocker_errors <- check_repair_file_blockers(parsed)
  if (length(blocker_errors) > 0L) {
    return(repair_vcf_result(
      errors = blocker_errors,
      warnings = list(),
      summary_extras = list(
        n_records_in = length(parsed$body_lines),
        n_records_out = 0L,
        output_path = NA_character_,
        step_summary = empty_step_summary(),
        steps = steps,
        chrom_style = chrom_style
      ),
      call = the_call
    ))
  }

  record_set <- extract_record_set(parsed)
  if (chrom_style != "auto") {
    record_set$chrom <- normalize_chrom_label(record_set$chrom, target = chrom_style)
  }
  n_records_in <- length(record_set$chrom)

  warnings <- list()
  step_rows <- list()
  for (step in steps) {
    n_before <- length(record_set$chrom)
    res <- apply_repair_step(record_set, step)
    record_set <- res$record_set
    n_after <- length(record_set$chrom)
    n_changed <- res$n_changed

    step_rows[[length(step_rows) + 1L]] <- list(
      step = step,
      n_in = as.integer(n_before),
      n_out = as.integer(n_after),
      n_changed = as.integer(n_changed)
    )
    if (n_changed > 0L) {
      warnings[[length(warnings) + 1L]] <- list(
        code = sprintf("repair_vcf_%s_applied", step),
        message = sprintf(
          "Step '%s': %d record(s) in -> %d record(s) out (%d changed/dropped).",
          step, n_before, n_after, n_changed
        ),
        where = list(
          step = step, n_changed = as.integer(n_changed),
          n_in = as.integer(n_before), n_out = as.integer(n_after)
        )
      )
    }
  }

  step_summary <- step_summary_tibble(step_rows)
  n_records_out <- length(record_set$chrom)

  if (n_records_out == 0L) {
    errors <- list(list(
      code = "repair_vcf_no_records_left",
      message = "Every record was dropped; no VCF written.",
      where = list()
    ))
    return(repair_vcf_result(
      errors = errors,
      warnings = warnings,
      summary_extras = list(
        n_records_in = n_records_in,
        n_records_out = 0L,
        output_path = NA_character_,
        step_summary = step_summary,
        steps = steps,
        chrom_style = chrom_style
      ),
      call = the_call
    ))
  }

  tmp <- tempfile(fileext = ".vcf")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)

  write_vcf_v4_2(
    output_path = tmp,
    chrom = record_set$chrom,
    pos = as.integer(record_set$pos),
    id = record_set$id,
    ref = record_set$ref,
    alt = record_set$alt,
    gt_matrix = record_set$gt_matrix,
    samples = colnames(record_set$gt_matrix),
    overwrite = TRUE
  )

  validation <- validate_vcf(tmp, target = "beagle")
  if (!validation$ok) {
    errors <- lapply(seq_len(nrow(validation$errors)), function(i) {
      list(
        code = paste0("repair_vcf_remaining_", validation$errors$code[i]),
        message = sprintf(
          "Post-repair validation failed: %s",
          validation$errors$message[i]
        ),
        where = validation$errors$where[[i]]
      )
    })
    return(repair_vcf_result(
      errors = errors,
      warnings = warnings,
      summary_extras = list(
        n_records_in = n_records_in,
        n_records_out = n_records_out,
        output_path = NA_character_,
        step_summary = step_summary,
        steps = steps,
        chrom_style = chrom_style
      ),
      call = the_call
    ))
  }

  write_exact_lines(readLines(tmp, warn = FALSE), output_vcf)

  repair_vcf_result(
    errors = list(),
    warnings = warnings,
    summary_extras = list(
      n_records_in = n_records_in,
      n_records_out = n_records_out,
      output_path = output_vcf,
      step_summary = step_summary,
      steps = steps,
      chrom_style = chrom_style
    ),
    call = the_call
  )
}

# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
empty_step_summary <- function() {
  tibble::tibble(
    step = character(),
    n_in = integer(),
    n_out = integer(),
    n_changed = integer()
  )
}

#' @keywords internal
#' @noRd
step_summary_tibble <- function(rows) {
  if (length(rows) == 0L) {
    return(empty_step_summary())
  }
  tibble::tibble(
    step = vapply(rows, `[[`, character(1), "step"),
    n_in = vapply(rows, `[[`, integer(1), "n_in"),
    n_out = vapply(rows, `[[`, integer(1), "n_out"),
    n_changed = vapply(rows, `[[`, integer(1), "n_changed")
  )
}

#' @keywords internal
#' @noRd
repair_vcf_result <- function(errors, warnings, summary_extras, call) {
  summary_list <- c(
    list(target = "beagle"),
    summary_extras
  )
  new_galick_result(
    errors = issues_tibble(errors),
    warnings = issues_tibble(warnings),
    summary = summary_list,
    call = call
  )
}

#' Detect file-level repair blockers
#'
#' Missing `#CHROM`, no samples, and duplicate samples are not safely
#' repairable; the validator-style codes here surface as `$errors` and the
#' caller writes nothing.
#'
#' @keywords internal
#' @noRd
check_repair_file_blockers <- function(parsed) {
  errors <- list()
  if (is.na(parsed$chrom_header_line)) {
    errors[[length(errors) + 1L]] <- list(
      code = "repair_vcf_missing_chrom_header",
      message = "VCF has no #CHROM header; cannot repair.",
      where = list()
    )
  }
  if (!is.na(parsed$chrom_header_line) && length(parsed$samples) == 0L) {
    errors[[length(errors) + 1L]] <- list(
      code = "repair_vcf_no_samples",
      message = "VCF has no sample columns; not repairable.",
      where = list()
    )
  }
  if (length(parsed$samples) > 0L && any(duplicated(parsed$samples))) {
    errors[[length(errors) + 1L]] <- list(
      code = "repair_vcf_duplicated_samples",
      message = "VCF has duplicate sample names; not repairable.",
      where = list(
        samples = unique(parsed$samples[duplicated(parsed$samples)])
      )
    )
  }
  errors
}

#' Build the in-memory record set used by repair steps
#'
#' Wrong-column-count rows still appear in the record set with
#' `col_count_ok = FALSE` and missing fields, so the
#' `drop_invalid_records` step can drop them explicitly.
#'
#' @keywords internal
#' @noRd
extract_record_set <- function(parsed) {
  n_expected <- length(parsed$header_cols)
  n_samples <- length(parsed$samples)
  n_rows <- length(parsed$body_lines)

  if (n_rows == 0L) {
    return(list(
      chrom = character(),
      pos = numeric(),
      id = character(),
      ref = character(),
      alt = character(),
      format = character(),
      col_count_ok = logical(),
      gt_matrix = matrix(
        character(),
        nrow = 0L,
        ncol = n_samples,
        dimnames = list(NULL, parsed$samples)
      )
    ))
  }

  col_count_ok <- parsed$body_col_counts == n_expected
  pluck <- function(idx) {
    vapply(parsed$body_split, function(x) {
      if (length(x) >= idx) x[idx] else NA_character_
    }, character(1))
  }

  chrom <- pluck(1L)
  pos_raw <- pluck(2L)
  id <- pluck(3L)
  ref <- pluck(4L)
  alt <- pluck(5L)
  format <- pluck(9L)

  pos_num <- suppressWarnings(as.numeric(pos_raw))

  gt_matrix <- matrix(
    NA_character_,
    nrow = n_rows,
    ncol = n_samples,
    dimnames = list(NULL, parsed$samples)
  )

  if (n_samples > 0L) {
    for (i in seq_len(n_rows)) {
      if (!col_count_ok[i] || length(parsed$body_split[[i]]) < 10L) next
      sample_cells <- parsed$body_split[[i]][10:n_expected]
      fmt <- format[i]
      if (is.na(fmt) || !nzchar(fmt)) next
      subfields <- strsplit(fmt, ":", fixed = TRUE)[[1]]
      gt_pos <- which(subfields == "GT")
      if (length(gt_pos) == 0L) next
      idx <- gt_pos[1]
      cell_subfields <- strsplit(sample_cells, ":", fixed = TRUE)
      gt_tokens <- vapply(cell_subfields, function(s) {
        if (length(s) >= idx) s[idx] else NA_character_
      }, character(1))
      gt_matrix[i, ] <- gt_tokens
    }
  }

  list(
    chrom = chrom,
    pos = pos_num,
    id = id,
    ref = ref,
    alt = alt,
    format = format,
    col_count_ok = col_count_ok,
    gt_matrix = gt_matrix
  )
}

#' @keywords internal
#' @noRd
subset_record_set <- function(record_set, keep) {
  list(
    chrom = record_set$chrom[keep],
    pos = record_set$pos[keep],
    id = record_set$id[keep],
    ref = record_set$ref[keep],
    alt = record_set$alt[keep],
    format = record_set$format[keep],
    col_count_ok = record_set$col_count_ok[keep],
    gt_matrix = record_set$gt_matrix[keep, , drop = FALSE]
  )
}

#' @keywords internal
#' @noRd
apply_repair_step <- function(record_set, step) {
  switch(step,
    drop_invalid_records = step_drop_invalid_records(record_set),
    filter_biallelic = step_filter_biallelic(record_set),
    fix_marker_ids = step_fix_marker_ids(record_set),
    remove_duplicates = step_remove_duplicates(record_set),
    sort = step_sort(record_set),
    cli::cli_abort(
      "Unknown repair step: {.val {step}}",
      class = "galick_error_repair_vcf_invalid_input"
    )
  )
}

#' @keywords internal
#' @noRd
step_drop_invalid_records <- function(record_set) {
  n <- length(record_set$chrom)
  if (n == 0L) {
    return(list(record_set = record_set, n_changed = 0L))
  }

  bad_struct <- !record_set$col_count_ok

  bad_pos <- is.na(record_set$pos) |
    !is.finite(record_set$pos) |
    record_set$pos <= 0 |
    abs(record_set$pos - round(record_set$pos)) > .Machine$double.eps^0.5

  bad_chrom <- is.na(record_set$chrom) | !nzchar(trimws(record_set$chrom))

  bad_ref <- is.na(record_set$ref) |
    !nzchar(trimws(record_set$ref)) |
    grepl("[^ACGTNacgtn,*]", record_set$ref)

  bad_alt <- is.na(record_set$alt) |
    !nzchar(trimws(record_set$alt)) |
    vapply(strsplit(record_set$alt, ",", fixed = TRUE), function(a) {
      length(a) == 0L ||
        any(is.na(a)) ||
        any(!nzchar(trimws(a))) ||
        any(grepl("[^ACGTNacgtn*]", a))
    }, logical(1))

  ref_eq_alt <- !is.na(record_set$ref) &
    !is.na(record_set$alt) &
    !grepl(",", record_set$alt, fixed = TRUE) &
    toupper(record_set$ref) == toupper(record_set$alt)

  bad_gt <- gt_invalid_mask(record_set$gt_matrix)

  invalid <- bad_struct | bad_pos | bad_chrom | bad_ref | bad_alt |
    ref_eq_alt | bad_gt
  n_dropped <- sum(invalid)

  list(
    record_set = subset_record_set(record_set, !invalid),
    n_changed = as.integer(n_dropped)
  )
}

#' @keywords internal
#' @noRd
step_filter_biallelic <- function(record_set) {
  if (length(record_set$alt) == 0L) {
    return(list(record_set = record_set, n_changed = 0L))
  }
  multi <- !is.na(record_set$alt) & grepl(",", record_set$alt, fixed = TRUE)
  list(
    record_set = subset_record_set(record_set, !multi),
    n_changed = as.integer(sum(multi))
  )
}

#' @keywords internal
#' @noRd
step_fix_marker_ids <- function(record_set) {
  if (length(record_set$id) == 0L) {
    return(list(record_set = record_set, n_changed = 0L))
  }
  missing <- is.na(record_set$id) |
    record_set$id == "." |
    !nzchar(trimws(record_set$id))
  if (any(missing)) {
    record_set$id[missing] <- sprintf(
      "%s_%s_%s_%s",
      record_set$chrom[missing],
      format(record_set$pos[missing], scientific = FALSE, trim = TRUE),
      record_set$ref[missing],
      record_set$alt[missing]
    )
  }
  list(
    record_set = record_set,
    n_changed = as.integer(sum(missing))
  )
}

#' @keywords internal
#' @noRd
step_remove_duplicates <- function(record_set) {
  n <- length(record_set$chrom)
  if (n == 0L) {
    return(list(record_set = record_set, n_changed = 0L))
  }

  ids <- record_set$id
  non_missing_ids <- !is.na(ids) & ids != "." & nzchar(trimws(ids))
  dup_id <- non_missing_ids & duplicated(ids)

  pair_keys <- paste(
    record_set$chrom,
    format(record_set$pos, scientific = FALSE, trim = TRUE),
    sep = ":"
  )
  dup_pair <- duplicated(pair_keys)

  drop <- dup_id | dup_pair
  list(
    record_set = subset_record_set(record_set, !drop),
    n_changed = as.integer(sum(drop))
  )
}

#' @keywords internal
#' @noRd
step_sort <- function(record_set) {
  n <- length(record_set$chrom)
  if (n == 0L) {
    return(list(record_set = record_set, n_changed = 0L))
  }

  order_idx <- order(record_set$chrom, record_set$pos)
  was_unsorted <- !identical(order_idx, seq_along(order_idx))
  list(
    record_set = subset_record_set(record_set, order_idx),
    n_changed = if (was_unsorted) as.integer(sum(order_idx != seq_along(order_idx))) else 0L
  )
}

#' Per-row mask of GT-related invalidity
#'
#' Treats as invalid: any sample missing the GT subfield, malformed GT
#' tokens (per the strict allele-token regex), or non-diploid GT.
#'
#' @keywords internal
#' @noRd
gt_invalid_mask <- function(gt_matrix) {
  n <- nrow(gt_matrix)
  if (n == 0L || ncol(gt_matrix) == 0L) {
    return(rep(FALSE, n))
  }
  gt_token_re <- "^(\\.|[0-9]+)([/|](\\.|[0-9]+))*$"
  out <- logical(n)
  for (i in seq_len(n)) {
    tokens <- gt_matrix[i, ]
    if (any(is.na(tokens))) {
      out[i] <- TRUE
      next
    }
    if (!all(grepl(gt_token_re, tokens))) {
      out[i] <- TRUE
      next
    }
    ploidy <- vapply(strsplit(tokens, "[/|]"), length, integer(1))
    if (any(ploidy != 2L)) {
      out[i] <- TRUE
    }
  }
  out
}
