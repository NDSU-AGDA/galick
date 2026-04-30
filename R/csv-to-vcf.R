#' Convert a CSV genotype table to a VCF v4.2 file
#'
#' Runs [validate_genotype_table()] and [validate_marker_map()] internally,
#' applies M4-specific stricter checks, and only writes a VCF if zero errors
#' were collected. Otherwise returns a `beagler_result(ok = FALSE)` with all
#' errors aggregated (no partial output).
#'
#' M4-specific checks applied on top of M2 and M3:
#' - Markers in the CSV that are not in the marker map → error (escalated
#'   from M3's warning, since they have no `ref`/`alt`).
#' - `ref` / `alt` outside `{A, C, G, T}` → error (v0.1.0 supports biallelic
#'   SNPs only).
#' - `output_vcf` extension must be `.vcf` (M4 does not compress; M7 handles
#'   BGZIP).
#' - `output_vcf` exists and `overwrite = FALSE` → error.
#'
#' @param csv Path to a CSV file or a `data.frame` already in M2 wide format.
#' @param marker_map Path (TSV/CSV) or `data.frame` with `marker_id`, `chrom`,
#'   `pos`, `ref`, `alt` columns.
#' @param output_vcf Path ending in `.vcf` to write the converted VCF.
#' @param sample_col Sample-id column name. `NULL` (default) auto-detects per
#'   M2's rule: first character/factor column with no NA, 100% unique, and
#'   tokens not recognized as a genotype encoding. Numeric sample IDs require
#'   explicit `sample_col`.
#' @param encoding,ploidy,missing_codes Passed through to
#'   [validate_genotype_table()].
#' @param chrom_style Passed through to [validate_marker_map()].
#' @param ahb_ref REF-mapping convention for AHB encoding: `"A"` (default)
#'   means data `A` codes homo-REF; pass `"B"` to swap.
#' @param aabb_ref_homo Mirror of `ahb_ref` for AABB encoding.
#' @param overwrite Refuse to clobber `output_vcf` unless `TRUE`.
#'
#' @return A `beagler_result` S3 object. On success, `summary$output_path`
#'   is the written file. On failure, `output_path` is `NA_character_` and
#'   no VCF is written.
#' @export
csv_to_vcf <- function(csv,
                       marker_map,
                       output_vcf,
                       sample_col = NULL,
                       encoding = c("auto", "AHB", "012", "centered",
                                    "AABB", "VCF", "REFHET", "IUPAC"),
                       ploidy = 2L,
                       chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
                       ahb_ref = c("A", "B"),
                       aabb_ref_homo = c("AA", "BB"),
                       missing_codes = c(NA, ".", "--", "", "N", "NN", "./."),
                       overwrite = FALSE) {
  the_call <- match.call()
  encoding <- match.arg(encoding)
  chrom_style <- match.arg(chrom_style)
  ahb_ref <- match.arg(ahb_ref)
  aabb_ref_homo <- match.arg(aabb_ref_homo)

  if (length(output_vcf) != 1L || !is.character(output_vcf) || is.na(output_vcf)) {
    cli::cli_abort(
      "{.arg output_vcf} must be a single non-NA character string.",
      class = "beagler_error_csv_to_vcf_invalid_input"
    )
  }

  csv_data <- read_csv_input(csv)
  map_data <- read_marker_map_input(marker_map)

  gt_result <- validate_genotype_table(
    csv_data,
    sample_col = sample_col,
    encoding = encoding,
    ploidy = ploidy,
    missing_codes = missing_codes
  )
  map_result <- validate_marker_map(
    map_data,
    genotype_table = csv_data,
    require_alleles = TRUE,
    chrom_style = chrom_style
  )

  resolved_sample_col <- if (is.null(sample_col)) {
    detect_sample_col(csv_data, missing_codes = missing_codes)
  } else {
    sample_col
  }
  csv_marker_names <- setdiff(names(csv_data), resolved_sample_col)
  map_marker_ids <- as.character(map_data$marker_id)

  m4_errors <- list()

  if (!grepl("\\.vcf$", output_vcf, ignore.case = TRUE)) {
    m4_errors[[length(m4_errors) + 1L]] <- list(
      code = "csv_to_vcf_unsupported_extension",
      message = sprintf("output_vcf must end in .vcf; got '%s'.", output_vcf),
      where = list(output_vcf = output_vcf)
    )
  }
  if (file.exists(output_vcf) && !overwrite) {
    m4_errors[[length(m4_errors) + 1L]] <- list(
      code = "csv_to_vcf_output_exists",
      message = sprintf("Output exists: '%s' (use overwrite = TRUE).", output_vcf),
      where = list(output_vcf = output_vcf)
    )
  }

  missing_in_map <- setdiff(csv_marker_names, map_marker_ids)
  if (length(missing_in_map) > 0L) {
    m4_errors[[length(m4_errors) + 1L]] <- list(
      code = "csv_to_vcf_unmappable_markers",
      message = sprintf(
        "%d marker(s) in CSV are missing from the marker map.",
        length(missing_in_map)
      ),
      where = list(missing = missing_in_map)
    )
  }

  if (all(c("ref", "alt") %in% names(map_data))) {
    ref_norm <- toupper(trimws(as.character(map_data$ref)))
    alt_norm <- toupper(trimws(as.character(map_data$alt)))
    bad_alleles <- !(ref_norm %in% c("A", "C", "G", "T")) |
      !(alt_norm %in% c("A", "C", "G", "T"))
    if (any(bad_alleles)) {
      m4_errors[[length(m4_errors) + 1L]] <- list(
        code = "csv_to_vcf_invalid_alleles",
        message = sprintf(
          "%d marker(s) have non-A/C/G/T alleles (biallelic SNPs only in v0.1.0).",
          sum(bad_alleles)
        ),
        where = list(rows = which(bad_alleles))
      )
    }
  }

  combined_errors <- rbind_issues(gt_result$errors, map_result$errors, issues_tibble(m4_errors))
  combined_warnings <- rbind_issues(gt_result$warnings, map_result$warnings)

  base_summary <- list(
    n_samples = gt_result$summary$n_samples,
    n_markers = gt_result$summary$n_markers,
    inferred_encoding = gt_result$summary$inferred_encoding,
    inferred_chrom_style = map_result$summary$inferred_chrom_style,
    output_path = NA_character_
  )

  if (nrow(combined_errors) > 0L) {
    return(new_beagler_result(
      errors = combined_errors,
      warnings = combined_warnings,
      summary = base_summary,
      call = the_call
    ))
  }

  inferred_encoding <- gt_result$summary$inferred_encoding
  ordered_idx <- match(csv_marker_names, map_marker_ids)
  ordered_chrom <- as.character(map_data$chrom)[ordered_idx]
  ordered_pos <- as.integer(suppressWarnings(as.numeric(as.character(map_data$pos))))[ordered_idx]
  ordered_ref <- as.character(map_data$ref)[ordered_idx]
  ordered_alt <- as.character(map_data$alt)[ordered_idx]

  gt_matrix <- build_gt_matrix(
    csv_data = csv_data,
    sample_col = resolved_sample_col,
    marker_names = csv_marker_names,
    refs = ordered_ref,
    alts = ordered_alt,
    encoding = inferred_encoding,
    ahb_ref = ahb_ref,
    aabb_ref_homo = aabb_ref_homo,
    missing_codes = missing_codes
  )

  write_vcf_v4_2(
    output_path = output_vcf,
    chrom = ordered_chrom,
    pos = ordered_pos,
    id = csv_marker_names,
    ref = ordered_ref,
    alt = ordered_alt,
    gt_matrix = gt_matrix,
    samples = as.character(csv_data[[resolved_sample_col]]),
    overwrite = overwrite
  )

  base_summary$output_path <- output_vcf
  new_beagler_result(
    errors = combined_errors,
    warnings = combined_warnings,
    summary = base_summary,
    call = the_call
  )
}

# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
read_csv_input <- function(csv) {
  if (inherits(csv, "data.frame")) {
    return(csv)
  }
  if (is.character(csv) && length(csv) == 1L && !is.na(csv) && file.exists(csv)) {
    return(as.data.frame(
      vroom::vroom(
        csv,
        col_types = vroom::cols(.default = vroom::col_character()),
        show_col_types = FALSE,
        progress = FALSE
      ),
      stringsAsFactors = FALSE
    ))
  }
  cli::cli_abort(
    "{.arg csv} must be a path to an existing CSV file or a data.frame.",
    class = "beagler_error_csv_to_vcf_invalid_input"
  )
}

#' @keywords internal
#' @noRd
read_marker_map_input <- function(marker_map) {
  if (inherits(marker_map, "data.frame")) {
    return(marker_map)
  }
  if (is.character(marker_map) && length(marker_map) == 1L && !is.na(marker_map) &&
      file.exists(marker_map)) {
    return(as.data.frame(
      vroom::vroom(
        marker_map,
        col_types = vroom::cols(.default = vroom::col_character()),
        show_col_types = FALSE,
        progress = FALSE
      ),
      stringsAsFactors = FALSE
    ))
  }
  cli::cli_abort(
    "{.arg marker_map} must be a path to an existing TSV/CSV file or a data.frame.",
    class = "beagler_error_csv_to_vcf_invalid_input"
  )
}

#' Recode each marker column with its per-marker REF/ALT (for IUPAC) and
#' assemble the n_markers x n_samples GT matrix.
#' @keywords internal
#' @noRd
build_gt_matrix <- function(csv_data, sample_col, marker_names,
                            refs, alts, encoding,
                            ahb_ref, aabb_ref_homo, missing_codes) {
  n_samples <- nrow(csv_data)
  n_markers <- length(marker_names)
  gt <- matrix(character(), nrow = n_markers, ncol = n_samples)
  for (i in seq_along(marker_names)) {
    col_values <- csv_data[[marker_names[i]]]
    gt[i, ] <- recode_to_vcf_gt(
      col_values,
      encoding = encoding,
      ref = refs[i],
      alt = alts[i],
      ahb_ref = ahb_ref,
      aabb_ref_homo = aabb_ref_homo,
      missing_codes = missing_codes
    )
  }
  colnames(gt) <- as.character(csv_data[[sample_col]])
  gt
}

#' rbind issues tibbles, treating any NULL or zero-row argument as empty.
#' @keywords internal
#' @noRd
rbind_issues <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, logical(1))]
  parts <- parts[vapply(parts, nrow, integer(1)) > 0L]
  if (length(parts) == 0L) {
    return(tibble::tibble(code = character(), message = character(), where = list()))
  }
  do.call(rbind, parts)
}
