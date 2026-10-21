#' Read verified GT-only marker metadata and a call matrix
#' @keywords internal
#' @noRd
read_genotype_records <- function(path) {
  x <- read_genotype_vcf(path)
  assert_resolvable_genotypes(x)
  list(
    metadata = data.frame(
      chrom = x$fixed$CHROM, pos = x$fixed$POS, id = x$fixed$ID,
      ref = x$fixed$REF, alt = x$fixed$ALT, stringsAsFactors = FALSE
    ),
    samples = x$samples,
    gt_matrix = x$gt
  )
}

#' Extract GT-only marker metadata + matrix from a parsed VCF (internal)
#'
#' Used by `vcf_to_dosage()`, `vcf_to_hardcalls()`, and other read-back
#' utilities. Filters to rows whose column count matches the header
#' (i.e. structurally valid rows). Returns the GT subfield only — any
#' other FORMAT subfields (DP, GQ, ...) are dropped.
#'
#' @keywords internal
#' @noRd
extract_vcf_records <- function(parsed) {
  n_expected <- length(parsed$header_cols)
  n_samples <- length(parsed$samples)

  empty_metadata <- data.frame(
    chrom = character(), pos = integer(), id = character(),
    ref = character(), alt = character(),
    stringsAsFactors = FALSE
  )
  if (n_expected < 9L || length(parsed$body_split) == 0L) {
    return(list(
      metadata = empty_metadata,
      samples = parsed$samples,
      gt_matrix = matrix(
        character(),
        nrow = 0L, ncol = n_samples,
        dimnames = list(NULL, parsed$samples)
      )
    ))
  }

  good_rows <- parsed$body_col_counts == n_expected
  if (!any(good_rows)) {
    return(list(
      metadata = empty_metadata,
      samples = parsed$samples,
      gt_matrix = matrix(
        character(),
        nrow = 0L, ncol = n_samples,
        dimnames = list(NULL, parsed$samples)
      )
    ))
  }

  body_mat <- do.call(rbind, parsed$body_split[good_rows])
  format_field <- body_mat[, 9]

  gt_matrix <- matrix(NA_character_, nrow = nrow(body_mat), ncol = n_samples)
  if (n_samples > 0L) {
    for (i in seq_len(nrow(body_mat))) {
      fmt <- format_field[i]
      if (is.na(fmt) || !nzchar(fmt)) next
      subfields <- strsplit(fmt, ":", fixed = TRUE)[[1]]
      gt_pos <- which(subfields == "GT")
      if (length(gt_pos) == 0L) next
      idx <- gt_pos[1]
      sample_cells <- body_mat[i, 10:n_expected]
      cell_subfields <- strsplit(sample_cells, ":", fixed = TRUE)
      gt_tokens <- vapply(cell_subfields, function(s) {
        if (length(s) >= idx) s[idx] else NA_character_
      }, character(1))
      gt_matrix[i, ] <- gt_tokens
    }
  }
  colnames(gt_matrix) <- parsed$samples

  list(
    metadata = data.frame(
      chrom = body_mat[, 1],
      pos = suppressWarnings(as.integer(body_mat[, 2])),
      id = body_mat[, 3],
      ref = body_mat[, 4],
      alt = body_mat[, 5],
      stringsAsFactors = FALSE
    ),
    samples = parsed$samples,
    gt_matrix = gt_matrix
  )
}

#' Pick a stable per-marker name (uses ID; falls back to chrom_pos_ref_alt)
#' @keywords internal
#' @noRd
marker_names_from_metadata <- function(metadata) {
  ifelse(
    is.na(metadata$id) | metadata$id == "." | !nzchar(trimws(metadata$id)),
    sprintf(
      "%s_%s_%s_%s",
      metadata$chrom,
      format(metadata$pos, scientific = FALSE, trim = TRUE),
      metadata$ref, metadata$alt
    ),
    metadata$id
  )
}

#' Convert a single GT token to a numeric ALT-allele count.
#'
#' Counts non-zero, non-missing allele indices. NA on missing or malformed
#' GT. Multiallelic GTs (e.g. "1/2") count both as ALT alleles, giving a
#' dosage of 2; this is the conservative interpretation when
#' biallelic/multiallelic disambiguation isn't possible from GT alone.
#'
#' @keywords internal
#' @noRd
gt_to_dosage <- function(gt) {
  if (is.na(gt) || !nzchar(gt)) {
    return(NA_real_)
  }
  alleles <- strsplit(gt, "[/|]")[[1]]
  if (any(alleles == ".")) {
    return(NA_real_)
  }
  als <- suppressWarnings(as.integer(alleles))
  if (any(is.na(als))) {
    return(NA_real_)
  }
  as.numeric(sum(als > 0L))
}

#' Convert a single GT token to a hardcall under the given encoding.
#' @keywords internal
#' @noRd
gt_to_hardcall <- function(gt, encoding) {
  d <- gt_to_dosage(gt)
  if (is.na(d)) {
    return(NA_character_)
  }
  d <- as.integer(d)
  switch(encoding,
    `012` = as.character(d),
    AHB = c("A", "H", "B")[d + 1L],
    AABB = c("AA", "AB", "BB")[d + 1L],
    NA_character_
  )
}
