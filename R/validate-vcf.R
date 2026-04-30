#' Validate a VCF file
#'
#' Reads a plain `.vcf` file via the tolerant internal `parse_vcf_raw()` and
#' reports shape, structural, and semantic issues. Always-on checks cover the VCF
#' v4.2 essentials; `target = "beagle"` and `target = "plink"` add tool-
#' specific stricter checks.
#'
#' Programmer/shape problems (`NULL`, wrong type, unsupported extension)
#' raise immediately. Data-quality problems are returned in `$errors` /
#' `$warnings` so the result remains inspectable.
#'
#' @param vcf Path to a plain `.vcf` file. `.vcf.gz` is rejected in v0.1.0
#'   (M5 does not decompress; M7 handles BGZF).
#' @param target `"general"`, `"beagle"`, or `"plink"`. Each adds tool-
#'   specific checks; some codes (`vcf_no_samples`, `vcf_multiallelic`)
#'   shift severity by target.
#' @param max_records_to_scan When finite, the validator scans the header
#'   in full and only this many body records. Sortedness and duplicate
#'   `(chrom, pos)` checks are skipped in partial mode (would false-positive
#'   at chunk boundaries).
#'
#' @return A `beagler_result` S3 object.
#' @export
validate_vcf <- function(vcf,
                         target = c("general", "beagle", "plink"),
                         max_records_to_scan = Inf) {
  the_call <- match.call()
  target <- match.arg(target)

  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a single non-NA character string.",
      class = "beagler_error_validate_vcf_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_validate_vcf_unsupported_extension"
    )
  }
  has_vcf_gz_extension <- grepl("\\.vcf\\.gz$", vcf, ignore.case = TRUE)

  parsed <- parse_vcf_raw(vcf)

  errors <- list()
  warnings <- list()

  if (is.na(parsed$fileformat_line)) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_missing_fileformat",
      message = "##fileformat line not found.",
      where = list()
    )
  } else {
    fmt <- sub("^##fileformat=", "", parsed$fileformat_line)
    if (!grepl("^VCFv4\\.", fmt)) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_invalid_fileformat",
        message = sprintf("Unsupported fileformat: %s", fmt),
        where = list(fileformat = fmt)
      )
    }
  }

  if (is.na(parsed$chrom_header_line)) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_missing_chrom_header",
      message = "#CHROM header line not found; cannot validate body.",
      where = list()
    )
    return(new_beagler_result(
      errors = issues_tibble(errors),
      warnings = issues_tibble(warnings),
      summary = list(
        target = target,
        n_samples = 0L,
        n_records = NA_integer_,
        n_records_scanned = NA_integer_,
        scanned_partial = FALSE,
        n_chroms = 0L,
        chromosome_summary = empty_chromosome_summary(),
        inferred_chrom_style = "unknown",
        has_vcf_gz_extension = has_vcf_gz_extension
      ),
      call = the_call
    ))
  }

  if (length(parsed$header_cols) < 9L) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_invalid_column_count",
      message = sprintf(
        "Header has %d columns; VCF requires at least 9.",
        length(parsed$header_cols)
      ),
      where = list(header_cols = parsed$header_cols)
    )
  }

  n_samples <- length(parsed$samples)
  no_sample_severity <- if (target == "general") "warning" else "error"
  if (n_samples == 0L) {
    issue <- list(
      code = "vcf_no_samples",
      message = "VCF has no sample columns.",
      where = list()
    )
    if (no_sample_severity == "warning") {
      warnings[[length(warnings) + 1L]] <- issue
    } else {
      errors[[length(errors) + 1L]] <- issue
    }
  }
  if (n_samples > 0L && any(duplicated(parsed$samples))) {
    for (s in unique(parsed$samples[duplicated(parsed$samples)])) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_duplicated_samples",
        message = sprintf("Duplicate sample '%s' in header.", s),
        where = list(sample = s)
      )
    }
  }

  n_records_total <- length(parsed$body_lines)
  n_records_scanned <- if (is.finite(max_records_to_scan)) {
    min(n_records_total, as.integer(max_records_to_scan))
  } else {
    n_records_total
  }
  scanned_partial <- n_records_scanned < n_records_total
  scan_idx <- seq_len(n_records_scanned)

  n_expected_cols <- length(parsed$header_cols)
  body_col_counts_scanned <- parsed$body_col_counts[scan_idx]
  bad_col_rows <- which(body_col_counts_scanned != n_expected_cols)
  if (length(bad_col_rows) > 0L) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_invalid_column_count",
      message = sprintf(
        "%d body row(s) do not have %d columns.",
        length(bad_col_rows), n_expected_cols
      ),
      where = list(rows = bad_col_rows)
    )
  }

  good_rows <- scan_idx[body_col_counts_scanned == n_expected_cols]
  if (length(good_rows) > 0L) {
    body_mat <- do.call(rbind, parsed$body_split[good_rows])
    chrom <- body_mat[, 1]
    pos_raw <- body_mat[, 2]
    id <- body_mat[, 3]
    ref <- body_mat[, 4]
    alt <- body_mat[, 5]
    format_field <- if (n_expected_cols >= 9L) body_mat[, 9] else rep("", nrow(body_mat))
    sample_cells <- if (n_samples > 0L) {
      body_mat[, 10:n_expected_cols, drop = FALSE]
    } else {
      matrix(character(), nrow = nrow(body_mat), ncol = 0)
    }

    pos_num <- suppressWarnings(as.numeric(pos_raw))
    bad_pos_mask <- is.na(pos_num) | !is.finite(pos_num) | pos_num <= 0 |
      abs(pos_num - round(pos_num)) > .Machine$double.eps^0.5
    if (any(bad_pos_mask)) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_invalid_pos",
        message = sprintf(
          "%d POS value(s) are NA, non-finite, non-positive, or non-integer.",
          sum(bad_pos_mask)
        ),
        where = list(rows = good_rows[bad_pos_mask])
      )
    }

    bad_chrom_mask <- is.na(chrom) | !nzchar(trimws(chrom))
    if (any(bad_chrom_mask)) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_invalid_chrom",
        message = sprintf("%d CHROM value(s) are empty.", sum(bad_chrom_mask)),
        where = list(rows = good_rows[bad_chrom_mask])
      )
    }

    missing_id_mask <- id == "." | !nzchar(id)
    if (any(missing_id_mask)) {
      warnings[[length(warnings) + 1L]] <- list(
        code = "vcf_missing_marker_ids",
        message = sprintf("%d record(s) have missing marker IDs (.).", sum(missing_id_mask)),
        where = list(rows = good_rows[missing_id_mask])
      )
    }

    valid_ids <- id[!missing_id_mask]
    dup_ids <- unique(valid_ids[duplicated(valid_ids)])
    for (d in dup_ids) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_duplicated_marker_ids",
        message = sprintf("Marker ID '%s' appears more than once.", d),
        where = list(marker_id = d)
      )
    }

    valid_pair_mask <- !bad_chrom_mask & !bad_pos_mask
    if (!scanned_partial && any(valid_pair_mask)) {
      pair_keys <- paste(
        chrom[valid_pair_mask],
        format(pos_num[valid_pair_mask], scientific = FALSE, trim = TRUE),
        sep = ":"
      )
      dup_pairs <- unique(pair_keys[duplicated(pair_keys)])
      for (k in dup_pairs) {
        errors[[length(errors) + 1L]] <- list(
          code = "vcf_duplicated_chrom_pos",
          message = sprintf("(chrom, pos) pair '%s' appears more than once.", k),
          where = list(chrom_pos = k)
        )
      }
    }

    bad_ref_mask <- is.na(ref) | !nzchar(trimws(ref)) |
      grepl("[^ACGTNacgtn,*]", ref)
    if (any(bad_ref_mask)) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_invalid_ref_alt",
        message = sprintf(
          "%d record(s) have empty or non-DNA REF.",
          sum(bad_ref_mask)
        ),
        where = list(rows = good_rows[bad_ref_mask])
      )
    }

    alt_lists <- strsplit(alt, ",", fixed = TRUE)
    multi_mask <- vapply(alt_lists, function(a) length(a) > 1L, logical(1))
    bad_alt_mask <- vapply(alt_lists, function(a) {
      length(a) == 0L ||
        any(is.na(a)) ||
        any(!nzchar(trimws(a))) ||
        any(grepl("[^ACGTNacgtn*]", a))
    }, logical(1))
    if (any(bad_alt_mask)) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_invalid_ref_alt",
        message = sprintf(
          "%d record(s) have empty or non-DNA ALT alleles.",
          sum(bad_alt_mask)
        ),
        where = list(rows = good_rows[bad_alt_mask])
      )
    }
    if (any(multi_mask)) {
      multi_severity <- if (target == "beagle") "error" else "warning"
      issue <- list(
        code = "vcf_multiallelic",
        message = sprintf(
          "%d record(s) have multiallelic ALT (more than one alternate).",
          sum(multi_mask)
        ),
        where = list(rows = good_rows[multi_mask])
      )
      if (multi_severity == "error") {
        errors[[length(errors) + 1L]] <- issue
      } else {
        warnings[[length(warnings) + 1L]] <- issue
      }
    }

    if (!scanned_partial && any(valid_pair_mask)) {
      sortedness_issues <- check_sortedness(
        chrom[valid_pair_mask],
        pos_num[valid_pair_mask],
        good_rows[valid_pair_mask]
      )
      warnings <- c(warnings, sortedness_issues)
    }

    format_issues <- check_format_and_gt(
      format_field, sample_cells, n_samples, good_rows
    )
    errors <- c(errors, format_issues$errors)
    warnings <- c(warnings, format_issues$warnings)

    if (target == "beagle") {
      beagle_issues <- check_beagle_specific(
        ref, alt_lists, sample_cells, format_field, n_samples, good_rows
      )
      errors <- c(errors, beagle_issues$errors)
      warnings <- c(warnings, beagle_issues$warnings)
    }
    if (target == "plink") {
      plink_issues <- check_plink_specific(chrom, id, missing_id_mask, good_rows)
      errors <- c(errors, plink_issues$errors)
      warnings <- c(warnings, plink_issues$warnings)
    }

    chrom_for_summary <- chrom[!bad_chrom_mask]
    pos_for_summary <- pos_num[!bad_chrom_mask]
    inferred_chrom_style <- detect_chrom_style(chrom_for_summary)
    chromosome_summary <- build_chromosome_summary(
      chrom_for_summary, pos_for_summary, scanned_partial
    )
  } else {
    chromosome_summary <- empty_chromosome_summary()
    inferred_chrom_style <- "unknown"
  }

  if (target == "beagle" && !has_vcf_gz_extension) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "vcf_beagle_not_bgzf",
      message = "Plain .vcf is not BGZF; bgzip + tabix-index it before using as ref=.",
      where = list()
    )
  }

  summary_list <- list(
    target = target,
    n_samples = n_samples,
    n_records = if (scanned_partial) NA_integer_ else as.integer(n_records_total),
    n_records_scanned = as.integer(n_records_scanned),
    scanned_partial = scanned_partial,
    n_chroms = nrow(chromosome_summary),
    chromosome_summary = chromosome_summary,
    inferred_chrom_style = inferred_chrom_style,
    has_vcf_gz_extension = has_vcf_gz_extension
  )

  new_beagler_result(
    errors = issues_tibble(errors),
    warnings = issues_tibble(warnings),
    summary = summary_list,
    call = the_call
  )
}

# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
empty_chromosome_summary <- function() {
  tibble::tibble(
    chrom = character(),
    n_records = integer(),
    n_unsorted = integer()
  )
}

#' @keywords internal
#' @noRd
build_chromosome_summary <- function(chrom, pos, scanned_partial) {
  if (length(chrom) == 0L) return(empty_chromosome_summary())
  unique_chroms <- unique(chrom)
  rows <- lapply(unique_chroms, function(ch) {
    ch_pos <- pos[chrom == ch]
    n_unsorted <- if (length(ch_pos) <= 1L || scanned_partial) {
      0L
    } else {
      as.integer(sum(diff(ch_pos) < 0))
    }
    list(
      chrom = ch,
      n_records = length(ch_pos),
      n_unsorted = n_unsorted
    )
  })
  out <- tibble::tibble(
    chrom = vapply(rows, `[[`, character(1), "chrom"),
    n_records = vapply(rows, `[[`, integer(1), "n_records"),
    n_unsorted = vapply(rows, `[[`, integer(1), "n_unsorted")
  )
  out[order(out$chrom), ]
}

#' @keywords internal
#' @noRd
check_sortedness <- function(chrom, pos, row_numbers) {
  warnings <- list()
  for (ch in unique(chrom)) {
    ch_mask <- chrom == ch
    ch_pos <- pos[ch_mask]
    if (length(ch_pos) <= 1L) next
    n_break <- sum(diff(ch_pos) < 0)
    if (n_break > 0L) {
      warnings[[length(warnings) + 1L]] <- list(
        code = "vcf_unsorted_records",
        message = sprintf(
          "Chromosome '%s' has %d adjacent order break(s).",
          ch, n_break
        ),
        where = list(chrom = ch, n_break = as.integer(n_break))
      )
    }
  }
  warnings
}

#' Validate FORMAT and per-sample GT subfields
#'
#' Parses GT out of FORMAT (e.g. `GT:DP` -> GT is subfield 1) per record;
#' rejects malformed allele tokens via a strict regex (digits or `.`,
#' separated by `/` or `|`).
#' @keywords internal
#' @noRd
check_format_and_gt <- function(format_field, sample_cells, n_samples, row_numbers) {
  errors <- list()
  warnings <- list()
  if (n_samples == 0L) {
    return(list(errors = errors, warnings = warnings))
  }

  format_subfields <- strsplit(format_field, ":", fixed = TRUE)
  gt_idx <- vapply(format_subfields, function(f) {
    pos <- which(f == "GT")
    if (length(pos) == 0L) NA_integer_ else as.integer(pos[1])
  }, integer(1))

  bad_format_mask <- is.na(format_field) | !nzchar(format_field) | is.na(gt_idx)
  if (any(bad_format_mask)) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_invalid_format_field",
      message = sprintf(
        "%d record(s) have FORMAT missing or lacking GT subfield.",
        sum(bad_format_mask)
      ),
      where = list(rows = row_numbers[bad_format_mask])
    )
  }

  gt_token_re <- "^(\\.|[0-9]+)([/|](\\.|[0-9]+))*$"
  bad_gt_rows <- integer()
  for (i in which(!bad_format_mask)) {
    expected_n <- length(format_subfields[[i]])
    cells <- sample_cells[i, ]
    cell_subfields <- strsplit(cells, ":", fixed = TRUE)
    inconsistent <- vapply(cell_subfields, function(s) {
      length(s) != expected_n
    }, logical(1))
    if (any(inconsistent)) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_inconsistent_sample_subfields",
        message = sprintf(
          "Record %d: %d sample(s) do not match FORMAT subfield count.",
          row_numbers[i], sum(inconsistent)
        ),
        where = list(row = row_numbers[i])
      )
      next
    }
    gt_tokens <- vapply(cell_subfields, function(s) s[gt_idx[i]], character(1))
    bad_gt <- !grepl(gt_token_re, gt_tokens)
    if (any(bad_gt)) {
      bad_gt_rows <- c(bad_gt_rows, row_numbers[i])
    }
  }
  if (length(bad_gt_rows) > 0L) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_invalid_gt",
      message = sprintf(
        "%d record(s) have malformed GT field(s).",
        length(unique(bad_gt_rows))
      ),
      where = list(rows = unique(bad_gt_rows))
    )
  }

  list(errors = errors, warnings = warnings)
}

#' @keywords internal
#' @noRd
check_beagle_specific <- function(ref, alt_lists, sample_cells, format_field,
                                  n_samples, row_numbers) {
  errors <- list()
  warnings <- list()

  same_alleles <- vapply(seq_along(alt_lists), function(i) {
    a <- alt_lists[[i]]
    !is.na(ref[i]) && length(a) == 1L && !is.na(a) && toupper(ref[i]) == toupper(a)
  }, logical(1))
  if (any(same_alleles)) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_beagle_ref_eq_alt",
      message = sprintf("%d record(s) have REF == ALT.", sum(same_alleles)),
      where = list(rows = row_numbers[same_alleles])
    )
  }

  if (n_samples > 0L) {
    format_subfields <- strsplit(format_field, ":", fixed = TRUE)
    gt_idx <- vapply(format_subfields, function(f) {
      pos <- which(f == "GT")
      if (length(pos) == 0L) NA_integer_ else as.integer(pos[1])
    }, integer(1))
    non_diploid_rows <- integer()
    for (i in seq_along(format_subfields)) {
      if (is.na(gt_idx[i])) next
      cells <- sample_cells[i, ]
      cell_subfields <- strsplit(cells, ":", fixed = TRUE)
      gt_tokens <- vapply(cell_subfields, function(s) {
        if (length(s) >= gt_idx[i]) s[gt_idx[i]] else NA_character_
      }, character(1))
      ploidy <- vapply(strsplit(gt_tokens, "[/|]"), length, integer(1))
      if (any(!is.na(gt_tokens) & ploidy != 2L)) {
        non_diploid_rows <- c(non_diploid_rows, row_numbers[i])
      }
    }
    if (length(non_diploid_rows) > 0L) {
      errors[[length(errors) + 1L]] <- list(
        code = "vcf_beagle_non_diploid",
        message = sprintf(
          "%d record(s) have GT(s) not of ploidy 2.",
          length(unique(non_diploid_rows))
        ),
        where = list(rows = unique(non_diploid_rows))
      )
    }
  }

  list(errors = errors, warnings = warnings)
}

#' @keywords internal
#' @noRd
check_plink_specific <- function(chrom, id, missing_id_mask, row_numbers) {
  errors <- list()
  warnings <- list()

  plink_chrom_re <- "^(chr)?([0-9]+|X|Y|XY|MT|M)$"
  unsupported_chrom_mask <- !grepl(plink_chrom_re, chrom, ignore.case = TRUE)
  if (any(unsupported_chrom_mask)) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "vcf_plink_unsupported_chrom",
      message = sprintf(
        "%d record(s) use chromosome labels PLINK may not accept.",
        sum(unsupported_chrom_mask)
      ),
      where = list(rows = row_numbers[unsupported_chrom_mask])
    )
  }

  if (any(missing_id_mask)) {
    errors[[length(errors) + 1L]] <- list(
      code = "vcf_plink_missing_ids",
      message = sprintf(
        "%d record(s) have missing marker IDs (PLINK requires unique non-empty IDs).",
        sum(missing_id_mask)
      ),
      where = list(rows = row_numbers[missing_id_mask])
    )
  }

  list(errors = errors, warnings = warnings)
}
