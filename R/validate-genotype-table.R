#' Validate a wide genotype table
#'
#' Checks shape, encoding, sample-id integrity, marker-name uniqueness,
#' missingness, and monomorphism. Returns a `galick_result` capturing
#' errors, warnings, and a structured summary.
#'
#' Programmer/input-shape problems (`NULL` data, wrong class, missing
#' `sample_col`, unsupported `ploidy`) raise immediately. Data-quality
#' problems (unsupported / mixed encoding, duplicate sample IDs or marker
#' names) are returned as `$errors` so the result remains inspectable.
#'
#' @param data A `data.frame` / `tibble` in wide format (rows = samples,
#'   columns = markers + a sample-id column).
#' @param sample_col Name of the sample-id column. If `NULL` (default),
#'   auto-detected as the first character/factor column with no missing
#'   values, 100% unique entries, and tokens that do not look like a
#'   genotype encoding. Numeric sample IDs require explicit `sample_col`.
#' @param encoding One of `"auto"`, `"AHB"`, `"012"`, `"centered"`, `"AABB"`,
#'   `"VCF"`, `"REFHET"`, `"IUPAC"`. `"centered"` is the `-1`/`0`/`1`
#'   marker-centered encoding common in plant breeding.
#' @param ploidy Integer; only `2L` is supported in v0.1.0.
#' @param missing_codes Tokens treated as missing before encoding detection.
#'
#' @return A `galick_result` S3 object.
#' @examples
#' calls <- data.frame(sample_id = c("plant_1", "plant_2", "plant_3"),
#'                     rs1 = c(0, 1, 2), rs2 = c(2, NA, 0))
#' validate_genotype_table(calls, sample_col = "sample_id", encoding = "012")
#' @export
validate_genotype_table <- function(data,
                                    sample_col = NULL,
                                    encoding = c(
                                      "auto", "AHB", "012", "centered",
                                      "AABB", "VCF", "REFHET", "IUPAC"
                                    ),
                                    ploidy = 2L,
                                    missing_codes = c(NA, ".", "--", "", "N", "NN", "./.")) {
  the_call <- match.call()
  encoding <- match.arg(encoding)

  if (is.null(data) || !inherits(data, "data.frame")) {
    cli::cli_abort(
      "{.arg data} must be a data.frame or tibble; got {.cls {class(data)[1]}}.",
      class = "galick_error_invalid_input"
    )
  }
  if (nrow(data) == 0L || ncol(data) < 2L) {
    cli::cli_abort(
      "{.arg data} must have at least one row and two columns (sample column + at least one marker).",
      class = "galick_error_invalid_input"
    )
  }
  if (length(ploidy) != 1L || !is.numeric(ploidy) || !identical(as.integer(ploidy), 2L)) {
    cli::cli_abort(
      "Only diploid (`ploidy = 2L`) is supported in v0.1.0; got {.val {ploidy}}.",
      class = "galick_error_unsupported_ploidy"
    )
  }

  sample_col <- resolve_sample_col(data, sample_col, missing_codes)

  parts <- split_sample_and_markers(data, sample_col)
  sample_ids <- parts$sample_ids
  marker_data <- parts$marker_data
  marker_names <- parts$marker_names

  errors <- list()
  warnings <- list()

  encoding_used <- if (encoding == "auto") {
    determine_encoding_auto(marker_data, missing_codes, errors)
  } else {
    determine_encoding_user(marker_data, encoding, missing_codes, errors)
  }
  errors <- attr(encoding_used, "errors_accum") %||% errors
  encoding_label <- as.character(encoding_used)

  dup_sample_issues <- duplicate_issues(
    sample_ids,
    code = "genotype_table_duplicated_samples",
    label = "Sample ID",
    field = "sample_id"
  )
  errors <- c(errors, dup_sample_issues)

  dup_marker_issues <- duplicate_issues(
    marker_names,
    code = "genotype_table_duplicated_markers",
    label = "Marker name",
    field = "marker_id"
  )
  errors <- c(errors, dup_marker_issues)

  missingness <- compute_missingness(marker_data, missing_codes)
  high_miss_threshold <- 0.5

  high_sample_issues <- high_missingness_issues(
    rates = missingness$by_sample,
    ids = sample_ids,
    threshold = high_miss_threshold,
    code = "genotype_table_high_missingness_sample",
    field = "sample_id"
  )
  warnings <- c(warnings, high_sample_issues)

  high_marker_issues <- high_missingness_issues(
    rates = missingness$by_marker,
    ids = marker_names,
    threshold = high_miss_threshold,
    code = "genotype_table_high_missingness_marker",
    field = "marker_id"
  )
  warnings <- c(warnings, high_marker_issues)

  monomorphic_idx <- which(count_distinct_per_marker(marker_data, missing_codes) <= 1L)
  for (j in monomorphic_idx) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "genotype_table_monomorphic_markers",
      message = sprintf("Marker '%s' is monomorphic.", marker_names[j]),
      where = list(marker_id = marker_names[j])
    )
  }

  summary_list <- list(
    n_samples = nrow(marker_data),
    n_markers = ncol(marker_data),
    inferred_encoding = encoding_label,
    inferred_ploidy = 2L,
    missing_rate_overall = missingness$overall,
    missing_rate_by_sample = stats::setNames(missingness$by_sample, sample_ids),
    missing_rate_by_marker = stats::setNames(missingness$by_marker, marker_names),
    n_monomorphic_markers = length(monomorphic_idx),
    n_duplicated_samples = length(unique(sample_ids[duplicated(sample_ids)])),
    n_duplicated_markers = length(unique(marker_names[duplicated(marker_names)]))
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

#' @keywords internal
#' @noRd
resolve_sample_col <- function(data, sample_col, missing_codes) {
  if (is.null(sample_col)) {
    detected <- detect_sample_col(data, missing_codes)
    if (is.null(detected)) {
      cli::cli_abort(
        c(
          "Could not auto-detect a sample-id column.",
          i = paste0(
            "Pass {.arg sample_col} explicitly. Auto-detection requires a ",
            "character/factor column with no missing values, 100% unique ",
            "entries, and tokens that do not look like a genotype encoding."
          )
        ),
        class = "galick_error_invalid_sample_col"
      )
    }
    return(detected)
  }
  if (!is.character(sample_col) || length(sample_col) != 1L || is.na(sample_col)) {
    cli::cli_abort(
      "{.arg sample_col} must be a single column name.",
      class = "galick_error_invalid_sample_col"
    )
  }
  if (!sample_col %in% names(data)) {
    cli::cli_abort(
      "{.arg sample_col} {.val {sample_col}} not found in data.",
      class = "galick_error_invalid_sample_col"
    )
  }
  sample_col
}

#' @keywords internal
#' @noRd
detect_sample_col <- function(data, missing_codes) {
  for (nm in names(data)) {
    col <- data[[nm]]
    if (!(is.character(col) || is.factor(col))) next
    vals <- as.character(col)
    if (any(is.na(vals)) || !all(nzchar(vals))) next
    if (length(unique(vals)) != length(vals)) next
    enc <- detect_genotype_encoding(vals, missing_codes = missing_codes)
    if (enc != "unknown") next
    return(nm)
  }
  NULL
}

#' Split a wide table into sample IDs and a marker-only data frame
#' Position-based to preserve duplicate marker names if any.
#' @keywords internal
#' @noRd
split_sample_and_markers <- function(data, sample_col) {
  sample_idx <- which(names(data) == sample_col)[1]
  marker_idx <- setdiff(seq_along(data), sample_idx)
  list(
    sample_ids = as.character(data[[sample_idx]]),
    marker_data = data[, marker_idx, drop = FALSE],
    marker_names = names(data)[marker_idx]
  )
}

#' Auto-detect encoding (encoding = "auto" branch)
#' @keywords internal
#' @noRd
determine_encoding_auto <- function(marker_data, missing_codes, errors) {
  all_tokens <- unlist(lapply(marker_data, as.character), use.names = FALSE)
  inferred <- detect_genotype_encoding(all_tokens, missing_codes = missing_codes)
  if (inferred == "unknown") {
    errors <- c(errors, list(list(
      code = "genotype_table_unsupported_encoding",
      message = "No supported genotype encoding contains the observed values.",
      where = list()
    )))
  } else if (inferred == "mixed") {
    errors <- c(errors, list(list(
      code = "genotype_table_mixed_encoding",
      message = "Observed values span multiple encodings; no single encoding fits.",
      where = list()
    )))
  }
  structure(inferred, errors_accum = errors)
}

#' Check user-specified encoding against observed tokens
#' @keywords internal
#' @noRd
determine_encoding_user <- function(marker_data, encoding, missing_codes, errors) {
  alphabet <- encoding_alphabets[[encoding]]
  all_tokens <- unlist(lapply(marker_data, as.character), use.names = FALSE)
  observed <- normalize_observed_tokens(all_tokens, missing_codes)
  out_of_alphabet <- setdiff(observed, alphabet)
  if (length(out_of_alphabet) > 0L) {
    sample_unexpected <- utils::head(out_of_alphabet, 5L)
    errors <- c(errors, list(list(
      code = "genotype_table_mixed_encoding",
      message = sprintf(
        "Specified encoding '%s' does not contain values: %s%s",
        encoding,
        paste(sample_unexpected, collapse = ", "),
        if (length(out_of_alphabet) > 5L) ", ..." else ""
      ),
      where = list(unexpected = list(out_of_alphabet))
    )))
  }
  structure(encoding, errors_accum = errors)
}

#' Build issues for repeated values
#' @keywords internal
#' @noRd
duplicate_issues <- function(values, code, label, field) {
  dups <- unique(values[duplicated(values)])
  lapply(dups, function(v) {
    issue <- list(
      code = code,
      message = sprintf("%s '%s' appears more than once.", label, v),
      where = list()
    )
    issue$where[[field]] <- v
    issue
  })
}

#' Compute overall, by-sample, by-marker missingness rates
#' @keywords internal
#' @noRd
compute_missingness <- function(marker_data, missing_codes) {
  miss_norm <- toupper(trimws(as.character(missing_codes)))
  is_cell_missing <- function(values) {
    v <- toupper(trimws(as.character(values)))
    is.na(v) | !nzchar(v) | v %in% miss_norm
  }
  per_col <- lapply(marker_data, is_cell_missing)
  m <- do.call(cbind, per_col)
  if (!is.matrix(m)) m <- matrix(m, nrow = nrow(marker_data))
  list(
    overall = if (length(m) == 0L) 0 else mean(m),
    by_sample = if (nrow(m) == 0L) numeric() else rowMeans(m),
    by_marker = if (ncol(m) == 0L) numeric() else colMeans(m)
  )
}

#' Build high-missingness warning issues
#' @keywords internal
#' @noRd
high_missingness_issues <- function(rates, ids, threshold, code, field) {
  flagged <- which(rates >= threshold)
  lapply(flagged, function(i) {
    issue <- list(
      code = code,
      message = sprintf(
        "%s '%s' is %.1f%% missing.",
        if (field == "sample_id") "Sample" else "Marker",
        ids[i],
        rates[i] * 100
      ),
      where = list()
    )
    issue$where[[field]] <- ids[i]
    issue$where$rate <- rates[i]
    issue
  })
}

#' Count distinct non-missing values per marker
#' @keywords internal
#' @noRd
count_distinct_per_marker <- function(marker_data, missing_codes) {
  miss_norm <- toupper(trimws(as.character(missing_codes)))
  vapply(marker_data, function(col) {
    v <- toupper(trimws(as.character(col)))
    v <- v[!(is.na(v) | !nzchar(v) | v %in% miss_norm)]
    length(unique(v))
  }, integer(1))
}

#' Coerce a list-of-issues into the canonical 3-column tibble
#' @keywords internal
#' @noRd
issues_tibble <- function(lst) {
  if (length(lst) == 0L) {
    return(tibble::tibble(code = character(), message = character(), where = list()))
  }
  tibble::tibble(
    code = vapply(lst, `[[`, character(1), "code"),
    message = vapply(lst, `[[`, character(1), "message"),
    where = lapply(lst, `[[`, "where")
  )
}

#' Local null-coalesce
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
