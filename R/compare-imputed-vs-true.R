#' Compare imputed calls or expected dosage with verified genotypes
#'
#' Aligns records by CHROM, POS, REF and ALT without guessing allele flips.
#' Missing predictions and absent output records remain in the eligible-cell
#' denominator. Truth is always GT; `field` selects prediction values.
#' Correlations are computed per marker; undefined cases are reported explicitly.
#'
#' @param imputed_vcf Path to the imputed VCF.
#' @param true_vcf Path to the truth VCF.
#' @param masked_vcf Optional complete, allele-matched truth copy with held-out
#'   GT cells marked missing. Only those cells are eligible.
#' @param field `"GT"` for ALT copy counts or `"DS"` for expected ALT dosage.
#'   DS never falls back to GT. Evaluation requires diploid biallelic calls.
#' @return A `beagler_result` with eligible, evaluated and missing counts,
#'   completion, unconditional and conditional GT concordance, RMSE,
#'   mean defined per-marker `dosage_r2`, `pooled_dosage_r2`, and `per_marker`.
#'   GT concordance is never obtained by rounding DS.
#' @export
compare_imputed_vs_true <- function(imputed_vcf, true_vcf, masked_vcf = NULL,
                                    field = c("GT", "DS")) {
  the_call <- match.call()
  field <- match.arg(field)
  for (arg in c("imputed_vcf", "true_vcf", "masked_vcf")) {
    value <- get(arg)
    if (arg == "masked_vcf" && is.null(value)) next
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
          !file.exists(value)) {
      cli::cli_abort("{.arg {arg}} must name an existing VCF.",
        class = "beagler_error_compare_imputed_invalid_input"
      )
    }
  }
  summary <- empty_compare_summary()
  summary$field <- field
  summary$used_mask_filter <- !is.null(masked_vcf)
  fail <- function(code, message) {
    new_beagler_result(
      errors = issues_tibble(list(list(
        code = code, message = message,
        where = list()
      ))), summary = summary, call = the_call
    )
  }
  truth <- read_genotype_vcf(true_vcf)
  imputed <- read_genotype_vcf(imputed_vcf)
  invalid <- comparison_contract_error(truth, imputed)
  if (!is.null(invalid)) {
    return(fail(invalid$code, invalid$message))
  }
  rows <- match(vcf_variant_keys(truth), vcf_variant_keys(imputed))
  columns <- match(truth$samples, imputed$samples)
  truth_values <- genotype_counts(truth)
  eligible <- !is.na(truth_values)
  if (!is.null(masked_vcf)) {
    mask <- read_genotype_vcf(masked_vcf)
    invalid <- comparison_contract_error(truth, mask)
    if (!is.null(invalid)) {
      return(fail(invalid$code, invalid$message))
    }
    mask_rows <- match(vcf_variant_keys(truth), vcf_variant_keys(mask))
    mask_columns <- match(truth$samples, mask$samples)
    if (anyNA(mask_rows) || anyNA(mask_columns)) {
      return(fail(
        "compare_imputed_incomplete_mask",
        "The mask must contain every truth variant and sample."
      ))
    }
    mask_counts <- genotype_counts(mask)[mask_rows, mask_columns, drop = FALSE]
    eligible <- eligible & is.na(mask_counts)
  }
  prediction_gt <- matrix(NA_real_, nrow(truth$gt), ncol(truth$gt))
  present_rows <- which(!is.na(rows))
  present_columns <- which(!is.na(columns))
  prediction_gt[present_rows, present_columns] <- genotype_counts(imputed)[
    rows[present_rows], columns[present_columns],
    drop = FALSE
  ]
  predictions <- prediction_gt
  if (field == "DS") {
    if (is.null(imputed$format_fields$DS)) {
      return(fail("compare_imputed_missing_ds", "Requested DS is absent."))
    }
    tokens <- imputed$format_fields$DS
    ds <- suppressWarnings(matrix(as.numeric(tokens), nrow(tokens), ncol(tokens)))
    if (any(tokens != "." & (!is.finite(ds) | ds < 0 | ds > 2))) {
      return(fail("compare_imputed_invalid_ds", "DS must be missing or between 0 and 2."))
    }
    predictions[, ] <- NA_real_
    predictions[present_rows, present_columns] <- ds[
      rows[present_rows], columns[present_columns],
      drop = FALSE
    ]
  }
  evaluated <- eligible & !is.na(predictions)
  gt_evaluated <- eligible & !is.na(prediction_gt)
  concordant <- gt_evaluated & prediction_gt == truth_values
  concordant[is.na(concordant)] <- FALSE
  summary$n_eligible <- sum(eligible)
  summary$n_compared <- sum(evaluated)
  summary$n_gt_compared <- sum(gt_evaluated)
  summary$n_missing_predictions <- sum(eligible & !evaluated)
  summary$n_concordant <- sum(concordant)
  summary$n_samples_compared <- length(present_columns)
  summary$n_markers_compared <- length(present_rows)
  summary$n_missing_markers <- sum(is.na(rows))
  summary$n_missing_samples <- sum(is.na(columns))
  if (!summary$n_eligible) {
    return(fail("compare_imputed_no_cells", "There are no eligible truth cells."))
  }
  summary$completion_rate <- summary$n_compared / summary$n_eligible
  summary$concordance_rate <- summary$n_concordant / summary$n_eligible
  summary$conditional_concordance_rate <- if (summary$n_gt_compared) {
    summary$n_concordant / summary$n_gt_compared
  } else {
    NA_real_
  }
  summary$rmse <- comparison_rmse(predictions[evaluated], truth_values[evaluated])
  summary$pooled_dosage_r2 <- comparison_r2(
    predictions[evaluated], truth_values[evaluated]
  )
  summary$per_marker <- tibble::tibble(
    chrom = truth$fixed$CHROM, pos = truth$fixed$POS, id = truth$fixed$ID,
    ref = truth$fixed$REF, alt = truth$fixed$ALT,
    n_eligible = rowSums(eligible), n_compared = rowSums(evaluated),
    n_missing_predictions = rowSums(eligible & !evaluated),
    dosage_r2 = vapply(seq_len(nrow(truth$gt)), function(i) {
      use <- evaluated[i, ]
      comparison_r2(predictions[i, use], truth_values[i, use])
    }, numeric(1)),
    rmse = vapply(seq_len(nrow(truth$gt)), function(i) {
      use <- evaluated[i, ]
      comparison_rmse(predictions[i, use], truth_values[i, use])
    }, numeric(1))
  )
  defined <- is.finite(summary$per_marker$dosage_r2)
  summary$n_undefined_r2 <- sum(!defined)
  summary$dosage_r2 <- if (any(defined)) {
    mean(summary$per_marker$dosage_r2[defined])
  } else {
    NA_real_
  }
  warnings <- if (summary$n_missing_predictions) {
    list(list(
      code = "compare_imputed_incomplete_predictions",
      message = "Some eligible truth cells have no prediction; inspect completion_rate.",
      where = list(n = summary$n_missing_predictions)
    ))
  } else {
    list()
  }
  new_beagler_result(
    warnings = issues_tibble(warnings),
    summary = summary, call = the_call
  )
}

#' @keywords internal
#' @noRd
vcf_variant_keys <- function(x) {
  paste(x$fixed$CHROM, format_vcf_pos(x$fixed$POS),
    x$fixed$REF, x$fixed$ALT,
    sep = "\r"
  )
}

#' @keywords internal
#' @noRd
comparison_contract_error <- function(truth, other) {
  metadata_error <- tryCatch(
    {
      merge_vcf_metadata(list(truth, other))
      NULL
    },
    error = function(e) conditionMessage(e)
  )
  if (!is.null(metadata_error)) {
    return(list(
      code = "compare_imputed_variant_mismatch",
      message = metadata_error
    ))
  }
  for (x in list(truth, other)) {
    assert_resolvable_genotypes(x)
    slots <- lengths(strsplit(as.character(x$gt), "[/|]"))
    if (any(grepl(",", x$fixed$ALT, fixed = TRUE)) ||
          any(slots != 2L & as.character(x$gt) != ".")) {
      return(list(
        code = "compare_imputed_unsupported_genotypes",
        message = "Evaluation requires diploid biallelic genotypes."
      ))
    }
    ids <- x$fixed$ID[!x$fixed$ID %in% c(".", "")]
    positions <- paste(x$fixed$CHROM, x$fixed$POS)
    if (anyDuplicated(positions) || anyDuplicated(ids)) {
      return(list(
        code = "compare_imputed_ambiguous_variants",
        message = "Duplicated variant identities or IDs make comparison ambiguous."
      ))
    }
  }
  by_id <- match(truth$fixed$ID, other$fixed$ID)
  named <- !truth$fixed$ID %in% c(".", "") & !is.na(by_id)
  by_position <- match(
    paste(truth$fixed$CHROM, truth$fixed$POS),
    paste(other$fixed$CHROM, other$fixed$POS)
  )
  shared <- !is.na(by_position)
  if (any(vcf_variant_keys(truth)[named] != vcf_variant_keys(other)[by_id[named]]) ||
        any(vcf_variant_keys(truth)[shared] != vcf_variant_keys(other)[by_position[shared]])) {
    return(list(
      code = "compare_imputed_variant_mismatch",
      message = "Shared IDs or positions disagree in coordinates or REF/ALT; no allele flips were inferred."
    ))
  }
  NULL
}

#' @keywords internal
#' @noRd
genotype_counts <- function(x) {
  matrix(vapply(as.character(x$gt), gt_to_dosage, numeric(1)),
    nrow(x$gt), ncol(x$gt),
    dimnames = dimnames(x$gt)
  )
}

#' @keywords internal
#' @noRd
comparison_r2 <- function(predicted, truth) {
  if (length(truth) < 2L || stats::sd(truth) == 0 || stats::sd(predicted) == 0) {
    return(NA_real_)
  }
  unname(stats::cor(predicted, truth)^2)
}

#' @keywords internal
#' @noRd
comparison_rmse <- function(predicted, truth) {
  if (!length(truth)) {
    return(NA_real_)
  }
  sqrt(mean((predicted - truth)^2))
}

#' @keywords internal
#' @noRd
empty_compare_summary <- function() {
  list(
    n_eligible = 0L, n_compared = 0L, n_gt_compared = 0L,
    n_concordant = 0L, n_missing_predictions = 0L,
    concordance_rate = NA_real_, conditional_concordance_rate = NA_real_,
    completion_rate = NA_real_, dosage_r2 = NA_real_, pooled_dosage_r2 = NA_real_,
    n_undefined_r2 = 0L, rmse = NA_real_, n_samples_compared = 0L,
    n_markers_compared = 0L, field = "GT", used_mask_filter = FALSE,
    per_marker = tibble::tibble()
  )
}
