#' Write a minimal VCF v4.2 file (internal)
#'
#' Used by [csv_to_vcf()] and [hapmap_to_vcf()] after validation has passed.
#' No compression; output is plain `.vcf`. M7 handles BGZIP / tabix.
#'
#' Header structure:
#' - `##fileformat=VCFv4.2`
#' - `##fileDate=YYYYMMDD`
#' - `##source=<source string>` (defaults to `galick/<version>`)
#' - `##contig=<ID=<chrom>>` per unique chromosome — length deliberately
#'   omitted because galick does not require a reference assembly.
#' - `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
#' - `#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples...>`
#'
#' Body: one record per marker; `QUAL=.`, `FILTER=PASS`, `INFO=.`, `FORMAT=GT`.
#'
#' @keywords internal
#' @noRd
write_vcf_v4_2 <- function(output_path,
                           chrom,
                           pos,
                           id,
                           ref,
                           alt,
                           gt_matrix,
                           samples = colnames(gt_matrix),
                           source = NULL,
                           reference = NULL,
                           file_date = format(Sys.Date(), "%Y%m%d"),
                           overwrite = FALSE) {
  if (length(output_path) != 1L || !is.character(output_path) || is.na(output_path)) {
    cli::cli_abort(
      "{.arg output_path} must be a single non-NA character string.",
      class = "galick_error_invalid_argument"
    )
  }
  if (!grepl("\\.vcf$", output_path, ignore.case = TRUE)) {
    cli::cli_abort(
      "Output path must end in {.val .vcf}; got {.file {output_path}}.",
      class = "galick_error_unsupported_extension"
    )
  }
  if (file.exists(output_path) && !overwrite) {
    cli::cli_abort(
      "Output file exists: {.file {output_path}} (use {.code overwrite = TRUE}).",
      class = "galick_error_output_exists"
    )
  }

  n_markers <- length(chrom)
  if (length(pos) != n_markers || length(id) != n_markers ||
        length(ref) != n_markers || length(alt) != n_markers) {
    cli::cli_abort(
      "Marker metadata vectors must all have the same length.",
      class = "galick_error_invalid_argument"
    )
  }
  if (!is.matrix(gt_matrix) || nrow(gt_matrix) != n_markers) {
    cli::cli_abort(
      "{.arg gt_matrix} must be a character matrix with one row per marker.",
      class = "galick_error_invalid_argument"
    )
  }
  if (length(samples) != ncol(gt_matrix)) {
    cli::cli_abort(
      "{.arg samples} length must match {.arg gt_matrix} column count.",
      class = "galick_error_invalid_argument"
    )
  }

  if (is.null(source)) {
    source <- sprintf("galick/%s", utils::packageVersion("galick"))
  }
  unique_chroms <- unique(chrom)

  header <- c(
    "##fileformat=VCFv4.2",
    sprintf("##fileDate=%s", file_date),
    sprintf("##source=%s", source),
    if (!is.null(reference) && !is.na(reference) && nzchar(reference)) {
      sprintf("##reference=%s", reference)
    },
    sprintf("##contig=<ID=%s>", unique_chroms),
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(
      c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
      collapse = "\t"
    )
  )

  records <- vapply(seq_len(n_markers), function(i) {
    paste(
      chrom[i], format_vcf_pos(pos[i]), id[i], ref[i], alt[i],
      ".", "PASS", ".", "GT",
      paste(gt_matrix[i, ], collapse = "\t"),
      sep = "\t"
    )
  }, character(1))

  writeLines(c(header, records), output_path)
  invisible(output_path)
}
