#' Compare an imputed VCF to a truth VCF
#'
#' Computes per-cell accuracy statistics on shared `(sample, marker)`
#' positions. Cells where either side is missing are excluded.
#' If `masked_vcf` is provided, the comparison is restricted to cells that
#' were missing in the masked VCF (i.e. the cells imputation actually had
#' to fill in).
#'
#' @param imputed_vcf Path to the imputed VCF.
#' @param true_vcf Path to the reference / truth VCF (typically the
#'   original VCF before [mask_genotypes()]).
#' @param masked_vcf Optional path to the masked VCF — if provided, only
#'   cells that were `./.` here AND non-missing in `true_vcf` count
#'   toward the comparison.
#' @return A `beagler_result` with these in `$summary`:
#'   `n_compared`, `n_concordant`, `concordance_rate`, `dosage_r2`,
#'   `rmse`.
#' @export
compare_imputed_vs_true <- function(imputed_vcf, true_vcf, masked_vcf = NULL) {
  the_call <- match.call()

  for (arg in c("imputed_vcf", "true_vcf")) {
    val <- get(arg)
    if (length(val) != 1L || !is.character(val) || is.na(val) || !file.exists(val)) {
      cli::cli_abort(
        "{.arg {arg}} must be a path to an existing VCF file.",
        class = "beagler_error_compare_imputed_invalid_input"
      )
    }
  }
  if (!is.null(masked_vcf)) {
    if (length(masked_vcf) != 1L || !is.character(masked_vcf) ||
        is.na(masked_vcf) || !file.exists(masked_vcf)) {
      cli::cli_abort(
        "{.arg masked_vcf} must be NULL or a path to an existing VCF file.",
        class = "beagler_error_compare_imputed_invalid_input"
      )
    }
  }

  imputed_d <- vcf_to_dosage(imputed_vcf)
  true_d <- vcf_to_dosage(true_vcf)

  common_samples <- intersect(rownames(imputed_d), rownames(true_d))
  common_markers <- intersect(colnames(imputed_d), colnames(true_d))
  if (length(common_samples) == 0L || length(common_markers) == 0L) {
    return(new_beagler_result(
      errors = issues_tibble(list(list(
        code = "compare_imputed_no_overlap",
        message = "No shared (sample, marker) cells between imputed and truth.",
        where = list()
      ))),
      warnings = issues_tibble(list()),
      summary = empty_compare_summary(),
      call = the_call
    ))
  }
  imputed_d <- imputed_d[common_samples, common_markers, drop = FALSE]
  true_d <- true_d[common_samples, common_markers, drop = FALSE]

  use_mask <- !is.null(masked_vcf)
  cell_in_mask <- if (use_mask) {
    masked_d <- vcf_to_dosage(masked_vcf)
    masked_d <- masked_d[
      intersect(common_samples, rownames(masked_d)),
      intersect(common_markers, colnames(masked_d)),
      drop = FALSE
    ]
    masked_aligned <- matrix(
      NA_real_, nrow = length(common_samples), ncol = length(common_markers),
      dimnames = list(common_samples, common_markers)
    )
    masked_aligned[rownames(masked_d), colnames(masked_d)] <- masked_d
    is.na(masked_aligned)
  } else {
    matrix(TRUE, nrow = length(common_samples), ncol = length(common_markers))
  }

  has_truth <- !is.na(true_d)
  has_imputed <- !is.na(imputed_d)
  considered <- cell_in_mask & has_truth & has_imputed

  n_compared <- sum(considered)
  if (n_compared == 0L) {
    return(new_beagler_result(
      errors = issues_tibble(list(list(
        code = "compare_imputed_no_cells",
        message = "No cells available for comparison.",
        where = list()
      ))),
      warnings = issues_tibble(list()),
      summary = empty_compare_summary(),
      call = the_call
    ))
  }

  imp_v <- imputed_d[considered]
  true_v <- true_d[considered]

  imp_hardcall <- as.integer(round(imp_v))
  n_concordant <- sum(imp_hardcall == as.integer(round(true_v)))
  concordance_rate <- n_concordant / n_compared
  rmse <- sqrt(mean((imp_v - true_v)^2))
  dosage_r2 <- if (stats::sd(true_v) > 0 && stats::sd(imp_v) > 0) {
    suppressWarnings(stats::cor(imp_v, true_v))^2
  } else {
    NA_real_
  }

  new_beagler_result(
    errors = issues_tibble(list()),
    warnings = issues_tibble(list()),
    summary = list(
      n_compared = n_compared,
      n_concordant = n_concordant,
      concordance_rate = concordance_rate,
      dosage_r2 = dosage_r2,
      rmse = rmse,
      n_samples_compared = length(common_samples),
      n_markers_compared = length(common_markers),
      used_mask_filter = use_mask
    ),
    call = the_call
  )
}

#' @keywords internal
#' @noRd
empty_compare_summary <- function() {
  list(
    n_compared = 0L,
    n_concordant = 0L,
    concordance_rate = NA_real_,
    dosage_r2 = NA_real_,
    rmse = NA_real_,
    n_samples_compared = 0L,
    n_markers_compared = 0L,
    used_mask_filter = FALSE
  )
}
