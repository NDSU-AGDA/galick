#' Randomly mask a fraction of non-missing genotypes
#'
#' Reads `vcf`, masks (sets to `./.`) a random fraction of the non-missing
#' GT cells, and writes the result to `output_vcf` as a minimal VCF v4.2.
#' Useful for benchmarking imputation accuracy: keep the original VCF as
#' the truth, run Beagle on the masked file, then compare with
#' [compare_imputed_vs_true()].
#'
#' Already-missing cells are not counted toward the fraction.
#'
#' @param vcf Path to the input VCF (`.vcf` or `.vcf.gz`).
#' @param output_vcf Path ending in `.vcf` for the masked VCF.
#' @param fraction Fraction of non-missing GTs to mask, in `(0, 1)`.
#' @param seed Optional integer seed for reproducible masking.
#' @param samples Optional character vector restricting which samples
#'   are eligible to be masked.
#' @param markers Optional character vector restricting which markers
#'   are eligible to be masked.
#' @param overwrite Refuse to clobber `output_vcf` unless `TRUE`.
#'
#' @return A `galick_result` with `summary$n_masked`,
#'   `summary$n_eligible`, and `summary$output_path`.
#' @examples
#' vcf <- system.file("extdata", "tiny_founders.vcf", package = "galick")
#' output <- tempfile(fileext = ".vcf")
#' result <- mask_genotypes(vcf, output, fraction = 0.2, seed = 42L)
#' result$summary$n_masked
#' unlink(output)
#' @export
mask_genotypes <- function(vcf,
                           output_vcf,
                           fraction = 0.1,
                           seed = NULL,
                           samples = NULL,
                           markers = NULL,
                           overwrite = FALSE) {
  the_call <- match.call()

  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf) || !file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a path to an existing VCF file.",
      class = "galick_error_mask_genotypes_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "galick_error_mask_genotypes_invalid_input"
    )
  }
  if (length(output_vcf) != 1L || !is.character(output_vcf) || is.na(output_vcf)) {
    cli::cli_abort(
      "{.arg output_vcf} must be a single non-NA character string.",
      class = "galick_error_mask_genotypes_invalid_input"
    )
  }
  if (!grepl("\\.vcf$", output_vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg output_vcf} must end in {.val .vcf}.",
      class = "galick_error_mask_genotypes_invalid_input"
    )
  }
  if (file.exists(output_vcf) && !overwrite) {
    cli::cli_abort(
      "{.arg output_vcf} exists: {.file {output_vcf}} (use {.code overwrite = TRUE}).",
      class = "galick_error_mask_genotypes_output_exists"
    )
  }
  if (!is.numeric(fraction) || length(fraction) != 1L ||
        is.na(fraction) || fraction <= 0 || fraction >= 1) {
    cli::cli_abort(
      "{.arg fraction} must be a single number strictly between 0 and 1.",
      class = "galick_error_mask_genotypes_invalid_input"
    )
  }

  rec <- read_genotype_records(vcf)

  if (nrow(rec$gt_matrix) == 0L) {
    cli::cli_abort(
      "{.arg vcf} has no records to mask.",
      class = "galick_error_mask_genotypes_invalid_input"
    )
  }

  sample_eligible <- if (is.null(samples)) {
    rep(TRUE, length(rec$samples))
  } else {
    rec$samples %in% samples
  }
  marker_names <- marker_names_from_metadata(rec$metadata)
  marker_eligible <- if (is.null(markers)) {
    rep(TRUE, length(marker_names))
  } else {
    marker_names %in% markers
  }

  is_missing <- is.na(rec$gt_matrix) |
    grepl(".", rec$gt_matrix, fixed = TRUE) |
    !nzchar(rec$gt_matrix)
  dim(is_missing) <- dim(rec$gt_matrix)

  eligibility_mask <- outer(marker_eligible, sample_eligible, FUN = "&") &
    !is_missing
  eligible_idx <- which(eligibility_mask, arr.ind = TRUE)
  n_eligible <- nrow(eligible_idx)

  if (!is.null(seed)) {
    withr::local_seed(as.integer(seed))
  }
  n_to_mask <- as.integer(round(fraction * n_eligible))
  n_to_mask <- max(0L, min(n_to_mask, n_eligible))

  if (n_to_mask > 0L) {
    sampled <- eligible_idx[sample.int(n_eligible, n_to_mask), , drop = FALSE]
    rec$gt_matrix[sampled] <- "./."
  }

  write_vcf_v4_2(
    output_path = output_vcf,
    chrom = rec$metadata$chrom,
    pos = rec$metadata$pos,
    id = marker_names,
    ref = rec$metadata$ref,
    alt = rec$metadata$alt,
    gt_matrix = rec$gt_matrix,
    samples = rec$samples,
    overwrite = overwrite
  )

  new_galick_result(
    errors = issues_tibble(list()),
    warnings = issues_tibble(list()),
    summary = list(
      n_eligible = n_eligible,
      n_masked = n_to_mask,
      fraction_requested = fraction,
      fraction_actual = if (n_eligible > 0L) n_to_mask / n_eligible else 0,
      output_path = output_vcf
    ),
    call = the_call
  )
}
