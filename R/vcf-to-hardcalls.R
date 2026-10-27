#' Extract a hard-call genotype matrix from a VCF
#'
#' Returns a character matrix of hard-call genotypes under the chosen
#' encoding. Rows are samples, columns are markers. `NA` for missing GT.
#'
#' Accepts both `.vcf` and `.vcf.gz`.
#'
#' @param vcf Path to a VCF file.
#' @param samples Optional character vector to subset rows.
#' @param markers Optional character vector to subset columns.
#' @param encoding `"012"` (default), `"AHB"`, or `"AABB"`. Mapping uses
#'   the dosage from the GT field: 0 alleles -> `"0"` / `"A"` / `"AA"`,
#'   1 -> `"1"` / `"H"` / `"AB"`, 2 -> `"2"` / `"B"` / `"BB"`.
#' @return A character matrix `samples x markers`.
#' @export
vcf_to_hardcalls <- function(vcf,
                             samples = NULL,
                             markers = NULL,
                             encoding = c("012", "AHB", "AABB")) {
  encoding <- match.arg(encoding)
  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a single non-NA character string.",
      class = "galick_error_vcf_to_hardcalls_invalid_input"
    )
  }
  if (!file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} not found: {.file {vcf}}",
      class = "galick_error_vcf_to_hardcalls_invalid_input"
    )
  }

  rec <- read_genotype_records(vcf)

  n_markers <- nrow(rec$gt_matrix)
  n_samples <- length(rec$samples)
  out <- matrix(NA_character_, nrow = n_samples, ncol = n_markers)
  if (n_markers > 0L && n_samples > 0L) {
    for (i in seq_len(n_markers)) {
      out[, i] <- vapply(
        rec$gt_matrix[i, ],
        function(gt) gt_to_hardcall(gt, encoding),
        character(1)
      )
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
