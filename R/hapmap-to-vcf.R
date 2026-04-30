#' Convert a TASSEL HapMap file to a VCF v4.2 file
#'
#' Parses the HapMap file, validates the embedded marker map (via
#' [validate_marker_map()]), applies M4-specific stricter checks
#' (biallelic A/C/G/T alleles only; per-cell IUPAC consistency with each
#' marker's REF/ALT), and only writes a VCF if zero errors were collected.
#'
#' @param hmp Path to a TASSEL-format HapMap file.
#' @param output_vcf Path ending in `.vcf` to write the converted VCF.
#' @param build Optional reference build / assembly name; written as
#'   `##reference=<build>` in the VCF header.
#' @param chrom_style Passed through to [validate_marker_map()].
#' @param missing_codes Tokens treated as missing.
#' @param overwrite Refuse to clobber `output_vcf` unless `TRUE`.
#'
#' @return A `beagler_result` S3 object.
#' @export
hapmap_to_vcf <- function(hmp,
                          output_vcf,
                          build = NA_character_,
                          chrom_style = c("auto", "numeric", "chr_prefixed", "as_is"),
                          missing_codes = c(NA, ".", "--", "", "N", "NN"),
                          overwrite = FALSE) {
  the_call <- match.call()
  chrom_style <- match.arg(chrom_style)

  if (length(output_vcf) != 1L || !is.character(output_vcf) || is.na(output_vcf)) {
    cli::cli_abort(
      "{.arg output_vcf} must be a single non-NA character string.",
      class = "beagler_error_hapmap_to_vcf_invalid_input"
    )
  }

  hmp_data <- read_hapmap(hmp)

  map_result <- validate_marker_map(
    hmp_data$metadata,
    require_alleles = TRUE,
    chrom_style = chrom_style
  )

  m4_errors <- list()

  if (!grepl("\\.vcf$", output_vcf, ignore.case = TRUE)) {
    m4_errors[[length(m4_errors) + 1L]] <- list(
      code = "hapmap_to_vcf_unsupported_extension",
      message = sprintf("output_vcf must end in .vcf; got '%s'.", output_vcf),
      where = list(output_vcf = output_vcf)
    )
  }
  if (file.exists(output_vcf) && !overwrite) {
    m4_errors[[length(m4_errors) + 1L]] <- list(
      code = "hapmap_to_vcf_output_exists",
      message = sprintf("Output exists: '%s' (use overwrite = TRUE).", output_vcf),
      where = list(output_vcf = output_vcf)
    )
  }

  ref_norm <- toupper(trimws(as.character(hmp_data$metadata$ref)))
  alt_norm <- toupper(trimws(as.character(hmp_data$metadata$alt)))
  bad_alleles <- !(ref_norm %in% c("A", "C", "G", "T")) |
    !(alt_norm %in% c("A", "C", "G", "T")) |
    is.na(ref_norm) | is.na(alt_norm) |
    (!is.na(ref_norm) & !is.na(alt_norm) & ref_norm == alt_norm)
  if (any(bad_alleles)) {
    m4_errors[[length(m4_errors) + 1L]] <- list(
      code = "hapmap_to_vcf_invalid_alleles",
      message = sprintf(
        "%d marker(s) have non-A/C/G/T or duplicate alleles (biallelic SNPs only).",
        sum(bad_alleles)
      ),
      where = list(rows = which(bad_alleles))
    )
  }

  iupac_errors <- list()
  if (!any(bad_alleles)) {
    iupac_errors <- check_iupac_consistency(
      hmp_data$gt_matrix,
      ref_norm,
      alt_norm,
      hmp_data$metadata$marker_id,
      missing_codes
    )
  }

  combined_errors <- rbind_issues(
    map_result$errors,
    issues_tibble(c(m4_errors, iupac_errors))
  )
  combined_warnings <- rbind_issues(map_result$warnings)

  base_summary <- list(
    n_samples = length(hmp_data$samples),
    n_markers = nrow(hmp_data$metadata),
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

  decoded <- matrix(
    character(),
    nrow = nrow(hmp_data$gt_matrix),
    ncol = ncol(hmp_data$gt_matrix)
  )
  for (i in seq_len(nrow(hmp_data$gt_matrix))) {
    decoded[i, ] <- recode_to_vcf_gt(
      hmp_data$gt_matrix[i, ],
      encoding = "IUPAC",
      ref = ref_norm[i],
      alt = alt_norm[i],
      missing_codes = missing_codes
    )
  }
  colnames(decoded) <- hmp_data$samples

  write_vcf_v4_2(
    output_path = output_vcf,
    chrom = as.character(hmp_data$metadata$chrom),
    pos = as.integer(hmp_data$metadata$pos),
    id = as.character(hmp_data$metadata$marker_id),
    ref = ref_norm,
    alt = alt_norm,
    gt_matrix = decoded,
    samples = hmp_data$samples,
    reference = if (is.na(build)) NULL else build,
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

#' Check that every non-missing cell in gt_matrix is consistent with its
#' marker's REF / ALT under IUPAC encoding.
#' @keywords internal
#' @noRd
check_iupac_consistency <- function(gt_matrix, ref_norm, alt_norm,
                                    marker_ids, missing_codes) {
  miss_norm <- toupper(trimws(as.character(missing_codes)))
  iupac_het <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")

  issues <- list()
  for (i in seq_len(nrow(gt_matrix))) {
    ref <- ref_norm[i]
    alt <- alt_norm[i]
    if (is.na(ref) || is.na(alt)) next
    het_pair <- paste0(sort(c(ref, alt)), collapse = "")
    het_code <- unname(iupac_het[het_pair])
    expected <- c(ref, alt, het_code)

    cells <- toupper(trimws(as.character(gt_matrix[i, ])))
    is_missing <- is.na(cells) | !nzchar(cells) | cells %in% miss_norm
    bad <- !is_missing & !(cells %in% expected)

    if (any(bad)) {
      issues[[length(issues) + 1L]] <- list(
        code = "hapmap_to_vcf_invalid_iupac_token",
        message = sprintf(
          "Marker %s (ref=%s, alt=%s) has %d cell(s) inconsistent with IUPAC alphabet.",
          marker_ids[i], ref, alt, sum(bad)
        ),
        where = list(
          marker_id = marker_ids[i],
          rows = which(bad)
        )
      )
    }
  }
  issues
}
