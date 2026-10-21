#' Validate compatibility between a target VCF and a reference panel VCF
#'
#' Compares two VCFs for chromosome-label consistency, shared marker
#' positions, and REF/ALT allele agreement at those shared positions.
#' Useful as a pre-flight check before passing a reference panel to
#' [run_beagle()].
#'
#' Both inputs may be plain `.vcf` or `.vcf.gz`.
#'
#' @param target_vcf Path to the target VCF.
#' @param reference_vcf Path to the reference-panel VCF.
#' @param chrom_style How to interpret chromosome labels for comparison.
#'   `"auto"` keeps labels as-is and warns on mismatched style; `"numeric"`
#'   / `"chr_prefixed"` normalize before comparison; `"as_is"` is identical
#'   to `"auto"` minus the warning.
#' @param require_phased Require a nonempty diploid, phased, fully called
#'   reference and disjoint sample IDs, as used by [run_beagle()].
#'
#' @return A `beagler_result` S3 object.
#' @export
validate_reference_match <- function(target_vcf,
                                     reference_vcf,
                                     chrom_style = c(
                                       "auto", "numeric",
                                       "chr_prefixed", "as_is"
                                     ),
                                     require_phased = FALSE) {
  the_call <- match.call()
  chrom_style <- match.arg(chrom_style)

  validate_vcf_path(target_vcf, "target_vcf")
  validate_vcf_path(reference_vcf, "reference_vcf")

  parsed_t <- parse_vcf_raw(target_vcf)
  parsed_r <- parse_vcf_raw(reference_vcf)

  errors <- list()
  warnings <- list()
  for (path in c(target_vcf, reference_vcf)) {
    checked <- validate_vcf(path, target = if (require_phased) "beagle" else "general")
    for (i in seq_len(nrow(checked$errors))) {
      errors[[length(errors) + 1L]] <- list(
        code = paste0("reference_match_", checked$errors$code[i]),
        message = checked$errors$message[i], where = list(path = path)
      )
    }
  }
  header_error <- tryCatch(
    {
      merge_vcf_metadata(list(list(meta = parsed_t$meta), list(meta = parsed_r$meta)))
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(header_error)) {
    errors[[length(errors) + 1L]] <- list(
      code = "reference_match_header_conflict",
      message = header_error, where = list()
    )
  }
  declared <- vapply(list(parsed_t, parsed_r), function(x) {
    any(grepl("^##reference=|^##contig=<.*assembly=", x$meta))
  }, logical(1))
  if (!all(declared)) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "reference_match_assembly_unverified",
      message = "Assembly identity is unverified without matching assembly metadata; no FASTA was checked.",
      where = list()
    )
  }
  overlap <- intersect(parsed_t$samples, parsed_r$samples)
  if (length(overlap)) {
    issue <- list(
      code = "reference_match_sample_overlap",
      message = "Target and reference share sample IDs; exclude target truth from the panel.",
      where = list(samples = overlap)
    )
    if (require_phased) {
      errors[[length(errors) + 1L]] <- issue
    } else {
      warnings[[length(warnings) + 1L]] <- issue
    }
  }
  if (require_phased && !length(errors)) {
    reference <- read_genotype_vcf(reference_vcf)
    if (!nrow(reference$gt) || !ncol(reference$gt) ||
          any(!grepl("^[01]\\|[01]$", reference$gt))) {
      errors[[length(errors) + 1L]] <- list(
        code = "reference_match_incomplete_or_unphased",
        message = "The reference must contain complete, phased diploid genotypes.",
        where = list()
      )
    }
  }
  if (length(errors)) {
    return(new_beagler_result(
      errors = issues_tibble(errors),
      warnings = issues_tibble(warnings), summary = empty_reference_match_summary(),
      call = the_call
    ))
  }

  if (is.na(parsed_t$chrom_header_line) || is.na(parsed_r$chrom_header_line)) {
    errors[[length(errors) + 1L]] <- list(
      code = "reference_match_invalid_vcf",
      message = "One or both VCFs are missing a #CHROM header.",
      where = list()
    )
    return(new_beagler_result(
      errors = issues_tibble(errors),
      warnings = issues_tibble(warnings),
      summary = empty_reference_match_summary(),
      call = the_call
    ))
  }

  rec_t <- extract_match_records(parsed_t)
  rec_r <- extract_match_records(parsed_r)

  norm_t <- maybe_normalize_chrom(rec_t$chrom, chrom_style)
  norm_r <- maybe_normalize_chrom(rec_r$chrom, chrom_style)

  style_t <- detect_chrom_style(norm_t)
  style_r <- detect_chrom_style(norm_r)

  if (chrom_style == "auto" &&
        style_t != style_r &&
        !(style_t %in% c("unknown")) &&
        !(style_r %in% c("unknown"))) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "reference_match_chrom_style_mismatch",
      message = sprintf(
        paste0(
          "Chromosome label styles differ: target = '%s', reference = '%s'. ",
          "Pass chrom_style = 'numeric' or 'chr_prefixed' to normalize before comparison."
        ),
        style_t, style_r
      ),
      where = list(target_style = style_t, reference_style = style_r)
    )
  }

  shared_chroms <- intersect(unique(norm_t), unique(norm_r))
  if (length(shared_chroms) == 0L && length(norm_t) > 0L && length(norm_r) > 0L) {
    errors[[length(errors) + 1L]] <- list(
      code = "reference_match_no_shared_chromosomes",
      message = "Target and reference share no chromosomes after normalization.",
      where = list(
        target_chroms = unique(norm_t),
        reference_chroms = unique(norm_r)
      )
    )
  }

  key_t <- paste(norm_t, rec_t$pos, sep = ":")
  key_r <- paste(norm_r, rec_r$pos, sep = ":")
  if (anyDuplicated(key_t) || anyDuplicated(key_r)) {
    return(new_beagler_result(
      errors = issues_tibble(list(list(
        code = "reference_match_ambiguous_positions",
        message = "Repeated positions after chromosome normalization are ambiguous.",
        where = list()
      ))), warnings = issues_tibble(warnings),
      summary = empty_reference_match_summary(), call = the_call
    ))
  }
  shared_keys <- intersect(key_t, key_r)
  n_shared <- length(shared_keys)

  if (n_shared == 0L && length(shared_chroms) > 0L) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "reference_match_no_shared_markers",
      message = "Target and reference share chromosomes but no (chrom, pos) markers.",
      where = list()
    )
  }

  n_ref_mismatch <- 0L
  n_alt_mismatch <- 0L
  if (n_shared > 0L) {
    idx_t <- match(shared_keys, key_t)
    idx_r <- match(shared_keys, key_r)
    ref_mm <- !is.na(rec_t$ref[idx_t]) & !is.na(rec_r$ref[idx_r]) &
      toupper(rec_t$ref[idx_t]) != toupper(rec_r$ref[idx_r])
    alt_mm <- !is.na(rec_t$alt[idx_t]) & !is.na(rec_r$alt[idx_r]) &
      toupper(rec_t$alt[idx_t]) != toupper(rec_r$alt[idx_r])
    n_ref_mismatch <- sum(ref_mm)
    n_alt_mismatch <- sum(alt_mm)
    if (n_ref_mismatch > 0L) {
      errors[[length(errors) + 1L]] <- list(
        code = "reference_match_ref_mismatch",
        message = sprintf(
          "%d shared marker(s) have mismatched REF alleles.", n_ref_mismatch
        ),
        where = list(n = n_ref_mismatch)
      )
    }
    if (n_alt_mismatch > 0L) {
      errors[[length(errors) + 1L]] <- list(
        code = "reference_match_alt_mismatch",
        message = sprintf(
          "%d shared marker(s) have mismatched ALT alleles.", n_alt_mismatch
        ),
        where = list(n = n_alt_mismatch)
      )
    }
  }

  summary_list <- list(
    n_chroms_target = length(unique(norm_t)),
    n_chroms_reference = length(unique(norm_r)),
    n_chroms_shared = length(shared_chroms),
    n_markers_target = length(key_t),
    n_markers_reference = length(key_r),
    n_markers_shared = n_shared,
    n_markers_target_only = length(setdiff(key_t, key_r)),
    n_markers_reference_only = length(setdiff(key_r, key_t)),
    n_ref_mismatch = n_ref_mismatch,
    n_alt_mismatch = n_alt_mismatch,
    inferred_chrom_style_target = style_t,
    inferred_chrom_style_reference = style_r,
    n_samples_target = length(parsed_t$samples),
    n_samples_reference = length(parsed_r$samples),
    shared_samples = overlap,
    assembly_metadata_present = all(declared)
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
validate_vcf_path <- function(path, arg) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    cli::cli_abort(
      "{.arg {arg}} must be a single non-NA character string.",
      class = "beagler_error_reference_match_invalid_input"
    )
  }
  if (!file.exists(path)) {
    cli::cli_abort(
      "{.arg {arg}} not found: {.file {path}}",
      class = "beagler_error_reference_match_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", path, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg {arg}} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_reference_match_invalid_input"
    )
  }
  invisible(path)
}

#' @keywords internal
#' @noRd
extract_match_records <- function(parsed) {
  n_expected <- length(parsed$header_cols)
  if (n_expected < 5L) {
    return(list(
      chrom = character(), pos = integer(),
      ref = character(), alt = character()
    ))
  }
  good_rows <- parsed$body_col_counts == n_expected
  if (!any(good_rows)) {
    return(list(
      chrom = character(), pos = integer(),
      ref = character(), alt = character()
    ))
  }
  body_mat <- do.call(rbind, parsed$body_split[good_rows])
  list(
    chrom = body_mat[, 1],
    pos = suppressWarnings(as.numeric(body_mat[, 2])),
    ref = body_mat[, 4],
    alt = body_mat[, 5]
  )
}

#' @keywords internal
#' @noRd
maybe_normalize_chrom <- function(chrom, chrom_style) {
  if (chrom_style %in% c("auto", "as_is")) {
    chrom
  } else {
    normalize_chrom_label(chrom, target = chrom_style)
  }
}

#' @keywords internal
#' @noRd
empty_reference_match_summary <- function() {
  list(
    n_chroms_target = 0L,
    n_chroms_reference = 0L,
    n_chroms_shared = 0L,
    n_markers_target = 0L,
    n_markers_reference = 0L,
    n_markers_shared = 0L,
    n_markers_target_only = 0L,
    n_markers_reference_only = 0L,
    n_ref_mismatch = 0L,
    n_alt_mismatch = 0L,
    inferred_chrom_style_target = "unknown",
    inferred_chrom_style_reference = "unknown",
    n_samples_target = 0L,
    n_samples_reference = 0L
  )
}
