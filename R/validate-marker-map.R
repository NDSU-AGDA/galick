#' Validate a marker map
#'
#' Checks shape, required columns, marker-id integrity, position and
#' chromosome validity, duplicates, sortedness within chromosomes, optional
#' allele columns, and (optionally) cross-check coverage against a genotype
#' table or a bare list of marker IDs.
#'
#' Programmer/shape problems (`NULL`, wrong class, missing required columns)
#' raise immediately. Data-quality problems are returned as `$errors`;
#' advisory issues as `$warnings`.
#'
#' @param map A `data.frame` / `tibble` with at minimum `marker_id`, `chrom`,
#'   and `pos` columns. Optional: `ref`, `alt`, `cm`.
#' @param genotype_table Optional `data.frame` in M2 wide format. The
#'   sample-id column is auto-detected the same way as in
#'   [validate_genotype_table()]; remaining column names are the marker set
#'   to cross-check against `map$marker_id`.
#' @param marker_ids Optional character vector of marker IDs to cross-check
#'   against `map$marker_id`. Mutually exclusive with `genotype_table`.
#' @param require_alleles When `TRUE`, require `ref` and `alt` columns to be
#'   present, non-empty, and distinct (`ref != alt` per marker).
#' @param chrom_style How to interpret `map$chrom` for duplicate/sortedness
#'   checks. `"auto"` infers via the M1 detector; `"numeric"` and
#'   `"chr_prefixed"` normalize first; `"as_is"` uses raw labels.
#'
#' @return A `galick_result` S3 object.
#' @examples
#' map <- data.frame(marker_id = c("rs1", "rs2"), chrom = "1",
#'                   pos = c(100L, 200L), ref = "A", alt = "G")
#' validate_marker_map(map, marker_ids = c("rs1", "rs2"), require_alleles = TRUE)
#' @export
validate_marker_map <- function(map,
                                genotype_table = NULL,
                                marker_ids = NULL,
                                require_alleles = FALSE,
                                chrom_style = c(
                                  "auto", "numeric",
                                  "chr_prefixed", "as_is"
                                )) {
  the_call <- match.call()
  chrom_style <- match.arg(chrom_style)

  if (is.null(map) || !inherits(map, "data.frame")) {
    cli::cli_abort(
      "{.arg map} must be a data.frame or tibble; got {.cls {class(map)[1]}}.",
      class = "galick_error_marker_map_invalid_input"
    )
  }
  if (nrow(map) == 0L) {
    cli::cli_abort(
      "{.arg map} has zero rows.",
      class = "galick_error_marker_map_invalid_input"
    )
  }
  required_cols <- c("marker_id", "chrom", "pos")
  missing_cols <- setdiff(required_cols, names(map))
  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "{.arg map} missing required columns: {.field {missing_cols}}",
      class = "galick_error_marker_map_missing_required_columns"
    )
  }
  if (!is.null(genotype_table) && !is.null(marker_ids)) {
    cli::cli_abort(
      "Pass either {.arg genotype_table} or {.arg marker_ids}, not both.",
      class = "galick_error_marker_map_invalid_input"
    )
  }

  cross_check_markers <- resolve_cross_check_markers(genotype_table, marker_ids)

  errors <- list()
  warnings <- list()

  marker_ids_in_map <- as.character(map$marker_id)
  chrom_raw <- as.character(map$chrom)
  pos_raw <- map$pos

  bad_id_mask <- is.na(marker_ids_in_map) | nchar(trimws(marker_ids_in_map)) == 0L
  if (any(bad_id_mask)) {
    errors[[length(errors) + 1L]] <- list(
      code = "marker_map_invalid_marker_id",
      message = sprintf(
        "%d marker_id value(s) are NA, empty, or whitespace-only.",
        sum(bad_id_mask)
      ),
      where = list(rows = which(bad_id_mask))
    )
  }

  pos_num <- suppressWarnings(as.numeric(as.character(pos_raw)))
  bad_pos_mask <- is.na(pos_num) | !is.finite(pos_num) | pos_num <= 0 |
    abs(pos_num - round(pos_num)) > .Machine$double.eps^0.5
  if (any(bad_pos_mask)) {
    errors[[length(errors) + 1L]] <- list(
      code = "marker_map_invalid_pos",
      message = sprintf(
        "%d pos value(s) are NA, non-finite, non-positive, or non-integer.",
        sum(bad_pos_mask)
      ),
      where = list(rows = which(bad_pos_mask))
    )
  }

  bad_chrom_mask <- is.na(chrom_raw) | nchar(trimws(chrom_raw)) == 0L
  if (any(bad_chrom_mask)) {
    errors[[length(errors) + 1L]] <- list(
      code = "marker_map_invalid_chrom",
      message = sprintf("%d chrom value(s) are NA or empty.", sum(bad_chrom_mask)),
      where = list(rows = which(bad_chrom_mask))
    )
  }

  inferred_style <- detect_chrom_style(chrom_raw[!bad_chrom_mask])

  if (chrom_style == "auto") {
    if (inferred_style == "unknown") {
      warnings[[length(warnings) + 1L]] <- list(
        code = "marker_map_unknown_chrom_style",
        message = "Chromosome labels do not match any recognized style.",
        where = list()
      )
    } else if (inferred_style == "mixed") {
      warnings[[length(warnings) + 1L]] <- list(
        code = "marker_map_mixed_chrom_style",
        message = "Chromosome labels mix multiple recognized styles.",
        where = list()
      )
    }
    chrom_used <- chrom_raw
  } else if (chrom_style == "as_is") {
    chrom_used <- chrom_raw
  } else {
    chrom_used <- normalize_chrom_label(chrom_raw, target = chrom_style)
  }

  valid_ids_mask <- !bad_id_mask
  valid_ids <- marker_ids_in_map[valid_ids_mask]
  dup_ids <- unique(valid_ids[duplicated(valid_ids)])
  for (id in dup_ids) {
    errors[[length(errors) + 1L]] <- list(
      code = "marker_map_duplicated_marker_ids",
      message = sprintf("marker_id '%s' appears more than once.", id),
      where = list(marker_id = id)
    )
  }

  valid_pair_mask <- !bad_chrom_mask & !bad_pos_mask
  dup_pairs <- character()
  if (any(valid_pair_mask)) {
    pair_keys <- paste(
      chrom_used[valid_pair_mask],
      format(pos_num[valid_pair_mask], scientific = FALSE, trim = TRUE),
      sep = ":"
    )
    dup_pairs <- unique(pair_keys[duplicated(pair_keys)])
    for (k in dup_pairs) {
      errors[[length(errors) + 1L]] <- list(
        code = "marker_map_duplicated_chrom_pos",
        message = sprintf("(chrom, pos) pair '%s' appears more than once.", k),
        where = list(chrom_pos = k)
      )
    }
  }

  chromosome_summary_rows <- list()
  total_unsorted <- 0L
  if (any(valid_pair_mask)) {
    chrom_for_sort <- chrom_used[valid_pair_mask]
    pos_for_sort <- pos_num[valid_pair_mask]
    unique_chroms <- unique(chrom_for_sort)
    for (ch in unique_chroms) {
      ch_mask <- chrom_for_sort == ch
      pos_in_ch <- pos_for_sort[ch_mask]
      n_unsorted <- if (length(pos_in_ch) <= 1L) 0L else sum(diff(pos_in_ch) < 0)
      chromosome_summary_rows[[length(chromosome_summary_rows) + 1L]] <- list(
        chrom = ch,
        n_markers = length(pos_in_ch),
        n_unsorted = as.integer(n_unsorted)
      )
      total_unsorted <- total_unsorted + as.integer(n_unsorted)
      if (n_unsorted > 0L) {
        warnings[[length(warnings) + 1L]] <- list(
          code = "marker_map_unsorted_within_chrom",
          message = sprintf(
            "Chromosome '%s' has %d adjacent order break(s).",
            ch, n_unsorted
          ),
          where = list(chrom = ch, n_unsorted = as.integer(n_unsorted))
        )
      }
    }
  }

  if (require_alleles) {
    allele_issues <- check_alleles(map)
    errors <- c(errors, allele_issues)
  }

  cross_check_results <- list(
    in_table_not_in_map = NULL,
    in_map_not_in_table = NULL,
    n_in_table = NULL
  )
  if (!is.null(cross_check_markers)) {
    cross_check_results$n_in_table <- length(cross_check_markers)
    cross_check_results$in_table_not_in_map <- setdiff(cross_check_markers, valid_ids)
    cross_check_results$in_map_not_in_table <- setdiff(valid_ids, cross_check_markers)
    if (length(cross_check_results$in_table_not_in_map) > 0L) {
      warnings[[length(warnings) + 1L]] <- list(
        code = "marker_map_markers_in_table_not_in_map",
        message = sprintf(
          "%d marker(s) in genotype table are missing from the map.",
          length(cross_check_results$in_table_not_in_map)
        ),
        where = list(missing = cross_check_results$in_table_not_in_map)
      )
    }
    if (length(cross_check_results$in_map_not_in_table) > 0L) {
      warnings[[length(warnings) + 1L]] <- list(
        code = "marker_map_markers_in_map_not_in_table",
        message = sprintf(
          "%d marker(s) in the map are missing from the genotype table.",
          length(cross_check_results$in_map_not_in_table)
        ),
        where = list(missing = cross_check_results$in_map_not_in_table)
      )
    }
  }

  chrom_summary_tbl <- if (length(chromosome_summary_rows) > 0L) {
    tbl <- tibble::tibble(
      chrom = vapply(chromosome_summary_rows, `[[`, character(1), "chrom"),
      n_markers = vapply(chromosome_summary_rows, `[[`, integer(1), "n_markers"),
      n_unsorted = vapply(chromosome_summary_rows, `[[`, integer(1), "n_unsorted")
    )
    tbl[order(tbl$chrom), ]
  } else {
    tibble::tibble(chrom = character(), n_markers = integer(), n_unsorted = integer())
  }

  summary_list <- list(
    n_markers_in_map = nrow(map),
    inferred_chrom_style = inferred_style,
    chromosome_summary = chrom_summary_tbl,
    n_duplicated_marker_ids = length(dup_ids),
    n_duplicated_chrom_pos = length(dup_pairs),
    n_unsorted_within_chrom = total_unsorted,
    markers_in_table_not_in_map = cross_check_results$in_table_not_in_map,
    markers_in_map_not_in_table = cross_check_results$in_map_not_in_table,
    n_markers_in_table = cross_check_results$n_in_table
  )

  new_galick_result(
    errors = issues_tibble(errors),
    warnings = issues_tibble(warnings),
    summary = summary_list,
    call = the_call
  )
}

# ---------------------------------------------------------------------------
# Helpers (internal)
# ---------------------------------------------------------------------------

#' Resolve the cross-check marker set from genotype_table and/or marker_ids
#' @keywords internal
#' @noRd
resolve_cross_check_markers <- function(genotype_table, marker_ids) {
  if (!is.null(marker_ids)) {
    if (!is.character(marker_ids)) {
      cli::cli_abort(
        "{.arg marker_ids} must be a character vector.",
        class = "galick_error_marker_map_invalid_input"
      )
    }
    return(marker_ids)
  }
  if (!is.null(genotype_table)) {
    if (!inherits(genotype_table, "data.frame")) {
      cli::cli_abort(
        "{.arg genotype_table} must be a data.frame or tibble.",
        class = "galick_error_marker_map_invalid_input"
      )
    }
    if (ncol(genotype_table) < 2L) {
      cli::cli_abort(
        "{.arg genotype_table} needs at least one sample column and one marker column.",
        class = "galick_error_marker_map_invalid_input"
      )
    }
    sample_col <- detect_sample_col(
      genotype_table,
      missing_codes = c(NA, ".", "--", "", "N", "NN")
    )
    if (is.null(sample_col)) {
      sample_col <- names(genotype_table)[1]
    }
    return(setdiff(names(genotype_table), sample_col))
  }
  NULL
}

#' Allele-validity issues when `require_alleles = TRUE`
#' @keywords internal
#' @noRd
check_alleles <- function(map) {
  issues <- list()
  if (!"ref" %in% names(map) || !"alt" %in% names(map)) {
    issues[[length(issues) + 1L]] <- list(
      code = "marker_map_missing_alleles",
      message = "require_alleles is TRUE but ref/alt columns are missing.",
      where = list()
    )
    return(issues)
  }

  ref_vals <- as.character(map$ref)
  alt_vals <- as.character(map$alt)
  bad_ref <- is.na(ref_vals) | nchar(trimws(ref_vals)) == 0L
  bad_alt <- is.na(alt_vals) | nchar(trimws(alt_vals)) == 0L
  bad_either <- bad_ref | bad_alt
  if (any(bad_either)) {
    issues[[length(issues) + 1L]] <- list(
      code = "marker_map_missing_alleles",
      message = sprintf(
        "%d marker(s) have missing or empty ref/alt.",
        sum(bad_either)
      ),
      where = list(rows = which(bad_either))
    )
  }

  same_alleles <- !bad_ref & !bad_alt &
    toupper(trimws(ref_vals)) == toupper(trimws(alt_vals))
  if (any(same_alleles)) {
    issues[[length(issues) + 1L]] <- list(
      code = "marker_map_missing_alleles",
      message = sprintf("%d marker(s) have ref == alt.", sum(same_alleles)),
      where = list(rows = which(same_alleles))
    )
  }

  issues
}
