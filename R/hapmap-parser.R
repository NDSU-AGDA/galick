#' Read a TASSEL-format HapMap file (internal)
#'
#' Parses the standard 11-column TASSEL HapMap header (`rs#`, `alleles`,
#' `chrom`, `pos`, `strand`, `assembly#`, `center`, `protLSID`, `assayLSID`,
#' `panel`, `QCcode`) followed by per-sample columns. Cell values are
#' single-character IUPAC codes.
#'
#' Returns a list with:
#' - `metadata`: data.frame with `marker_id`, `chrom`, `pos`, `ref`, `alt`
#'   columns (alleles parsed from the `alleles` field on `/`).
#' - `samples`: character vector of sample IDs.
#' - `gt_matrix`: character matrix, n_markers x n_samples, of raw IUPAC
#'   tokens (no decoding here; that's the caller's job).
#'
#' @keywords internal
#' @noRd
read_hapmap <- function(path) {
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

  raw <- as.data.frame(
    vroom::vroom(
      path,
      delim = "\t",
      col_types = vroom::cols(.default = vroom::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    ),
    stringsAsFactors = FALSE
  )

  required <- c("rs#", "alleles", "chrom", "pos")
  missing_required <- setdiff(required, names(raw))
  if (length(missing_required) > 0L) {
    cli::cli_abort(
      "HapMap file missing required column(s): {.field {missing_required}}",
      class = "beagler_error_hapmap_invalid_input"
    )
  }

  standard_meta <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#",
                     "center", "protLSID", "assayLSID", "panel", "QCcode")
  sample_cols <- setdiff(names(raw), standard_meta)
  if (length(sample_cols) == 0L) {
    cli::cli_abort(
      "HapMap file has no sample columns (after stripping standard metadata).",
      class = "beagler_error_hapmap_invalid_input"
    )
  }

  alleles_split <- strsplit(as.character(raw$alleles), "/", fixed = TRUE)
  ref_alt <- do.call(rbind, lapply(alleles_split, function(x) {
    if (length(x) >= 2L) c(x[1], x[2]) else c(NA_character_, NA_character_)
  }))

  gt_matrix <- as.matrix(raw[, sample_cols, drop = FALSE])
  rownames(gt_matrix) <- as.character(raw[["rs#"]])
  colnames(gt_matrix) <- sample_cols

  list(
    metadata = data.frame(
      marker_id = as.character(raw[["rs#"]]),
      chrom = as.character(raw$chrom),
      pos = as.integer(suppressWarnings(as.numeric(as.character(raw$pos)))),
      ref = ref_alt[, 1],
      alt = ref_alt[, 2],
      stringsAsFactors = FALSE
    ),
    samples = sample_cols,
    gt_matrix = gt_matrix
  )
}
