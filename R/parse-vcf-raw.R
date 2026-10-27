#' Tolerant raw VCF parser (internal)
#'
#' Reads a VCF file and returns the structural pieces without aborting on
#' shape problems (missing `#CHROM`, malformed column counts, short headers,
#' invalid POS). The validator [validate_vcf()] consumes this and turns the
#' diagnostics into `$errors` / `$warnings` codes.
#'
#' Aborts only on truly catastrophic conditions (path arg invalid, file not
#' found). Everything else is preserved as raw text + structural metadata
#' for the validator to interpret.
#'
#' @return A list with:
#'   - `path`: original path (normalized).
#'   - `meta`: character vector of `##` metadata lines (zero-length if none).
#'   - `fileformat_line`: the first `##fileformat=` line, or `NA_character_`.
#'   - `chrom_header_line`: the `#CHROM` line (raw), or `NA_character_`.
#'   - `header_cols`: split header columns (after stripping leading `#`),
#'     or `character(0)` if no header.
#'   - `samples`: header_cols beyond column 9, or `character(0)`.
#'   - `body_lines`: non-meta non-header non-blank lines, raw.
#'   - `body_split`: list of split-on-tab character vectors per body line.
#'   - `body_col_counts`: integer vector, one per body line.
#'
#' @keywords internal
#' @noRd
parse_vcf_raw <- function(path) {
  if (length(path) != 1L || !is.character(path) || is.na(path)) {
    cli::cli_abort(
      "{.arg path} must be a single non-NA character string.",
      class = "galick_error_invalid_argument"
    )
  }
  if (!file.exists(path)) {
    cli::cli_abort(
      "File not found: {.file {path}}",
      class = "galick_error_invalid_argument"
    )
  }

  lines <- read_vcf_lines(path)
  is_meta <- startsWith(lines, "##")
  is_header <- startsWith(lines, "#CHROM")
  is_blank <- !nzchar(lines)
  is_body <- !is_meta & !is_header & !is_blank

  meta <- lines[is_meta]
  fileformat_idx <- grep("^##fileformat=", meta)
  fileformat_line <- if (length(fileformat_idx) > 0L) meta[fileformat_idx[1]] else NA_character_

  chrom_header_lines <- lines[is_header]
  chrom_header_line <- if (length(chrom_header_lines) >= 1L) chrom_header_lines[1] else NA_character_

  header_cols <- if (!is.na(chrom_header_line)) {
    strsplit(sub("^#", "", chrom_header_line), "\t", fixed = TRUE)[[1]]
  } else {
    character(0)
  }

  samples <- if (length(header_cols) > 9L) header_cols[-(1:9)] else character(0)

  body_lines <- lines[is_body]
  body_split <- if (length(body_lines) > 0L) {
    strsplit(body_lines, "\t", fixed = TRUE)
  } else {
    list()
  }
  body_col_counts <- vapply(body_split, length, integer(1))

  list(
    path = path,
    meta = meta,
    fileformat_line = fileformat_line,
    chrom_header_line = chrom_header_line,
    header_cols = header_cols,
    samples = samples,
    body_lines = body_lines,
    body_split = body_split,
    body_col_counts = body_col_counts
  )
}

#' Read all lines from a `.vcf` or `.vcf.gz` file
#'
#' BGZF (Beagle's output format) is a gzip variant; standard `gzfile()` reads
#' it as a stream just fine. Random access (tabix) requires BGZF-aware
#' tools but is not needed for whole-file parsing.
#'
#' @keywords internal
#' @noRd
read_vcf_lines <- function(path) {
  if (grepl("\\.vcf\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, "rt")
    on.exit(close(con))
    readLines(con, warn = FALSE)
  } else {
    readLines(path, warn = FALSE)
  }
}
