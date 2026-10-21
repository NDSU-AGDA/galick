#' Read a minimal VCF (v4.2-ish, GT field only)
#'
#' Internal helper used by M4 round-trip tests; M5 will use it as the seed
#' for `validate_vcf()`. Parses only header metadata, sample columns, and
#' per-record fields. Does not interpret INFO; does not parse FORMAT fields
#' other than GT.
#'
#' @param path Path to a plain `.vcf` file.
#' @return A list with `meta` (character vector of `##` lines), `chrom`,
#'   `pos`, `id`, `ref`, `alt`, `qual`, `filter`, `info`, `format`,
#'   `samples`, and `gt_matrix` (n_markers x n_samples character matrix).
#' @keywords internal
#' @noRd
read_vcf_minimal <- function(path) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    cli::cli_abort(
      "{.arg path} must be a single non-NA character string.",
      class = "beagler_error_invalid_argument"
    )
  }
  if (!file.exists(path)) {
    cli::cli_abort(
      "File not found: {.file {path}}",
      class = "beagler_error_invalid_argument"
    )
  }

  lines <- readLines(path, warn = FALSE)
  is_meta <- startsWith(lines, "##")
  is_header <- startsWith(lines, "#CHROM")
  is_blank <- !nzchar(lines)

  meta <- lines[is_meta]
  header_lines <- lines[is_header]
  if (length(header_lines) != 1L) {
    cli::cli_abort(
      "Expected exactly one {.code #CHROM} header line; got {length(header_lines)}.",
      class = "beagler_error_invalid_argument"
    )
  }

  header_cols <- strsplit(sub("^#", "", header_lines), "\t", fixed = TRUE)[[1]]
  if (length(header_cols) < 9L) {
    cli::cli_abort(
      "Header line has {length(header_cols)} columns; minimum is 9.",
      class = "beagler_error_invalid_argument"
    )
  }
  samples <- header_cols[-(1:9)]

  body_idx <- which(!is_meta & !is_header & !is_blank)
  if (length(body_idx) == 0L) {
    return(list(
      meta = meta,
      chrom = character(),
      pos = integer(),
      id = character(),
      ref = character(),
      alt = character(),
      qual = character(),
      filter = character(),
      info = character(),
      format = character(),
      samples = samples,
      gt_matrix = matrix(
        character(),
        nrow = 0L, ncol = length(samples),
        dimnames = list(NULL, samples)
      )
    ))
  }

  body <- strsplit(lines[body_idx], "\t", fixed = TRUE)
  expected_cols <- length(samples) + 9L
  body_mat <- do.call(rbind, body)
  if (ncol(body_mat) != expected_cols) {
    cli::cli_abort(
      "Body row column count ({ncol(body_mat)}) does not match header ({expected_cols}).",
      class = "beagler_error_invalid_argument"
    )
  }

  gt <- body_mat[, 10:expected_cols, drop = FALSE]
  colnames(gt) <- samples

  list(
    meta = meta,
    chrom = body_mat[, 1],
    pos = as.integer(body_mat[, 2]),
    id = body_mat[, 3],
    ref = body_mat[, 4],
    alt = body_mat[, 5],
    qual = body_mat[, 6],
    filter = body_mat[, 7],
    info = body_mat[, 8],
    format = body_mat[, 9],
    samples = samples,
    gt_matrix = gt
  )
}
