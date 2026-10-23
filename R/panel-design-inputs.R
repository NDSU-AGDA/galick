#' @keywords internal
#' @noRd
panel_abort <- function(message) {
  cli::cli_abort(message, class = "beagler_error_panel_design_invalid_input")
}

#' @keywords internal
#' @noRd
panel_number <- function(x, name, lower = 0, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower ||
        (integer && (x != floor(x) || x > .Machine$integer.max))) {
    panel_abort(paste0("Invalid ", name, "."))
  }
}

#' Count cells without materializing the VCF genotype representation.
#' @keywords internal
#' @noRd
panel_cell_count <- function(path, limit) {
  validate_vcf_path(path)
  con <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  samples <- NULL
  records <- 0
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "#CHROM\t")) samples <- length(strsplit(line, "\t")[[1]]) - 9L
    if (nzchar(line) && !startsWith(line, "#")) {
      if (is.null(samples) || samples < 1L) panel_abort("VCF records require a sample header.")
      records <- records + 1
      if (records * samples > limit) panel_abort("max_genotype_cells exceeded before genotype allocation.")
    }
  }
  if (is.null(samples) || !records) panel_abort("A nonempty sample VCF is required.")
  records * samples
}

#' @keywords internal
#' @noRd
panel_marker_rows <- function(data, truth, name, empty = FALSE) {
  if (is.null(data) && empty) {
    return(integer())
  }
  cols <- c("chrom", "pos", "ref", "alt")
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(cols %in% names(data)) ||
        (!empty && !nrow(data)) || anyNA(data[cols]) || !is.numeric(data$pos) ||
        any(!is.finite(data$pos) | data$pos < 1 | data$pos != floor(data$pos))) {
    panel_abort(paste0(name, " requires unique exact chrom, pos, ref, alt identities."))
  }
  keys <- paste(data$chrom, format_vcf_pos(data$pos), data$ref, data$alt, sep = "\r")
  rows <- match(keys, vcf_variant_keys(truth))
  if (anyDuplicated(keys) || anyNA(rows)) panel_abort(paste0(name, " contains duplicate or unmatched markers."))
  sort(rows)
}

#' @keywords internal
#' @noRd
panel_samples <- function(data, truth) {
  cols <- c("sample_id", "partition")
  if (!is.data.frame(data) || anyDuplicated(names(data)) || !all(cols %in% names(data)) || !nrow(data)) {
    panel_abort("sample_data requires sample_id and partition.")
  }
  data <- as.data.frame(data[intersect(c(cols, "group", "block"), names(data))])
  data[] <- lapply(data, as.character)
  if (anyNA(data) || any(vapply(data, function(x) any(!nzchar(trimws(x))), logical(1))) ||
        anyDuplicated(data$sample_id) || any(!data$sample_id %in% truth$samples) ||
        !setequal(data$partition, c("design", "tune", "test"))) {
    panel_abort("Invalid sample partitions or IDs.")
  }
  if (!"group" %in% names(data)) data$group <- "all"
  if ("block" %in% names(data) && any(vapply(
    split(data$partition, data$block),
    function(x) length(unique(x)) != 1L, logical(1)
  ))) {
    panel_abort("A block cannot cross partitions.")
  }
  counts <- table(factor(data$partition, levels = c("design", "tune", "test")), data$group)
  if (any(counts["design", ] < 20L) || any(counts[c("tune", "test"), , drop = FALSE] < 10L)) {
    panel_abort("Every group needs at least 20 design, 10 tuning and 10 test samples.")
  }
  data[order(match(data$sample_id, truth$samples)), , drop = FALSE]
}

#' @keywords internal
#' @noRd
panel_inputs <- function(truth_vcf, reference_vcf, sample_data, candidate_markers, scoring_markers,
                         mandatory_markers, excluded_markers, max_genotype_cells) {
  used <- panel_cell_count(truth_vcf, max_genotype_cells)
  panel_cell_count(reference_vcf, max_genotype_cells - used)
  checked <- validate_reference_match(truth_vcf, reference_vcf, require_phased = TRUE)
  if (!checked$ok) panel_abort("Reference preflight failed; require phased, complete, disjoint reference calls.")
  truth <- read_genotype_vcf(truth_vcf)
  reference <- read_genotype_vcf(reference_vcf)
  for (x in list(truth, reference)) {
    if (any(!x$fixed$REF %in% c("A", "C", "G", "T")) ||
          any(!x$fixed$ALT %in% c("A", "C", "G", "T")) ||
          anyDuplicated(paste(x$fixed$CHROM, x$fixed$POS))) {
      panel_abort("Unique biallelic SNP positions required.")
    }
  }
  index <- order_vcf_rows(truth)
  truth$fixed <- data.table::copy(truth$fixed[index])
  truth$gt <- truth$gt[index, , drop = FALSE]
  truth$formats <- truth$formats[index]
  truth$format_fields <- lapply(truth$format_fields, function(x) x[index, , drop = FALSE])
  samples <- panel_samples(sample_data, truth)
  candidates <- panel_marker_rows(candidate_markers, truth, "candidate_markers")
  scoring <- panel_marker_rows(scoring_markers, truth, "scoring_markers")
  mandatory <- panel_marker_rows(mandatory_markers, truth, "mandatory_markers", TRUE)
  excluded <- panel_marker_rows(excluded_markers, truth, "excluded_markers", TRUE)
  if (length(intersect(candidates, scoring)) || any(!mandatory %in% candidates) ||
        length(intersect(mandatory, excluded))) {
    panel_abort("Conflicting candidate/scoring/mandatory/excluded sets.")
  }
  candidates <- setdiff(candidates, excluded)
  requested <- union(candidates, scoring)
  if (any(!vcf_variant_keys(truth)[requested] %in% vcf_variant_keys(reference))) {
    panel_abort("Every selectable and scoring marker must occur in the reference.")
  }
  counts <- genotype_counts(truth)
  design <- match(samples$sample_id[samples$partition == "design"], truth$samples)
  calls <- counts[candidates, design, drop = FALSE]
  af <- rowMeans(calls, na.rm = TRUE) / 2
  eligible <- rowMeans(!is.na(calls)) >= 0.95 & is.finite(af) & pmin(af, 1 - af) >= 0.01
  removed <- candidates[!eligible]
  if (any(mandatory %in% removed)) panel_abort("Mandatory markers fail design-only call-rate/MAF eligibility.")
  candidates <- candidates[eligible]
  if (!setequal(truth$fixed$CHROM[candidates], truth$fixed$CHROM[scoring])) {
    panel_abort("Candidate and scoring sets must cover the same chromosomes after eligibility filtering.")
  }
  for (partition in c("design", "tune", "test")) {
    for (group in unique(samples$group)) {
      cols <- match(samples$sample_id[samples$partition == partition & samples$group == group], truth$samples)
      if (!any(!is.na(counts[scoring, cols, drop = FALSE]))) panel_abort("Every partition/group needs scoring truth.")
    }
  }
  list(
    truth = truth, counts = counts, samples = samples, candidates = candidates, scoring = scoring,
    mandatory = mandatory, removed = removed, warnings = checked$warnings
  )
}

#' @keywords internal
#' @noRd
panel_markers <- function(input, rows) {
  x <- input$truth$fixed[sort(rows)]
  tibble::tibble(chrom = x$CHROM, pos = x$POS, ref = x$REF, alt = x$ALT)
}
