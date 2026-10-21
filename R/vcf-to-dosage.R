#' Extract a dosage matrix from a VCF
#'
#' Returns a numeric matrix of ALT-allele counts. Rows are samples, columns
#' are markers. Values are `0`, `1`, `2`, or `NA` (missing GT). Multiallelic
#' GTs (e.g. `1/2`) are counted as 2 ALT alleles — the conservative
#' interpretation when REF/ALT disambiguation isn't possible from GT alone.
#'
#' Marker column names use `ID` when present; markers with `ID = "."` get
#' `chrom_pos_ref_alt` fallbacks so the matrix has unique column names.
#'
#' Accepts both `.vcf` and `.vcf.gz`.
#'
#' @param vcf Path to a VCF file.
#' @param samples Optional character vector to subset rows (samples).
#' @param markers Optional character vector to subset columns (marker IDs).
#' @return A numeric matrix `samples x markers`.
#' @export
vcf_to_dosage <- function(vcf, samples = NULL, markers = NULL) {
  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a single non-NA character string.",
      class = "beagler_error_vcf_to_dosage_invalid_input"
    )
  }
  if (!file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} not found: {.file {vcf}}",
      class = "beagler_error_vcf_to_dosage_invalid_input"
    )
  }

  rec <- read_genotype_records(vcf)

  n_markers <- nrow(rec$gt_matrix)
  n_samples <- length(rec$samples)
  out <- matrix(NA_real_, nrow = n_samples, ncol = n_markers)
  if (n_markers > 0L && n_samples > 0L) {
    for (i in seq_len(n_markers)) {
      out[, i] <- vapply(rec$gt_matrix[i, ], gt_to_dosage, numeric(1))
    }
  }
  rownames(out) <- rec$samples
  colnames(out) <- marker_names_from_metadata(rec$metadata)

  if (!is.null(samples)) {
    keep <- intersect(samples, rownames(out))
    out <- out[keep, , drop = FALSE]
  }
  if (!is.null(markers)) {
    keep <- intersect(markers, colnames(out))
    out <- out[, keep, drop = FALSE]
  }
  out
}
