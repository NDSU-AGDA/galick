#' @keywords internal
#' @noRd
vcf_diagnostic_tables <- function() {
  list(
    allele_usage = tibble::tibble(
      record = integer(), allele_index = integer(), allele = character(), role = character(),
      supported = logical(), n_copies = double(), n_genotype_cells = integer()
    ),
    genotype_counts = tibble::tibble(
      record = integer(), genotype = character(), category = character(), n = integer()
    ),
    genotype_summary = tibble::tibble(
      genotype = character(), category = character(), n_records = integer(), n_genotype_cells = double()
    ),
    unsupported_alleles = tibble::tibble(
      allele = character(), n_records = integer(), n_genotype_cells = double(),
      n_copies = double(), records = list()
    ),
    quality_fields = tibble::tibble(
      field = character(), number = character(), type = character(), description = character(),
      profiled = logical()
    ),
    quality_summary = tibble::tibble(
      field = character(), category = character(), n = integer(), n_finite = integer(),
      n_missing = integer(), n_nonfinite = integer(), n_nonnumeric = integer(),
      min = double(), q25 = double(), median = double(), q75 = double(), max = double()
    ),
    quality_by_record = tibble::tibble(
      record = integer(), field = character(), category = character(), n = integer(),
      n_finite = integer(), median = double()
    ),
    quality_comparisons = tibble::tibble(
      field = character(), n_comparable_records = integer(), n_not_comparable_records = integer(),
      n_ratio_records = integer(), n_lower = integer(),
      n_equal = integer(), n_higher = integer(), median_difference = double(),
      median_positive_baseline_ratio = double()
    )
  )
}

#' @keywords internal
#' @noRd
vcf_diagnostic_bind <- function(rows, empty) {
  rows <- Filter(function(x) !is.null(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(empty)
  out <- tibble::as_tibble(data.table::rbindlist(rows, use.names = TRUE))
  if (length(names(empty))) out <- out[, names(empty)]
  out
}

#' @keywords internal
#' @noRd
vcf_diagnostic_fields <- function(meta) {
  lines <- meta[startsWith(meta, "##FORMAT=<") & endsWith(meta, ">")]
  rows <- lapply(lines, function(line) {
    contents <- substr(line, 11L, nchar(line) - 1L)
    parts <- suppressWarnings(tryCatch(utils::read.csv(
      text = contents, header = FALSE, colClasses = "character", check.names = FALSE,
      comment.char = "", strip.white = TRUE
    ), error = function(e) NULL))
    if (is.null(parts) || nrow(parts) != 1L) return(NULL)
    parts <- unlist(parts, use.names = FALSE)
    keys <- sub("=.*$", "", parts)
    values <- sub("^[^=]*=", "", parts)
    get <- function(key) {
      index <- match(key, keys)
      if (is.na(index)) NA_character_ else values[index]
    }
    id <- get("ID")
    if (is.na(id) || !nzchar(id)) return(NULL)
    tibble::tibble(
      field = id, number = get("Number"), type = get("Type"), description = get("Description"),
      profiled = identical(get("Number"), "1") && get("Type") %in% c("Integer", "Float")
    )
  })
  out <- vcf_diagnostic_bind(rows, vcf_diagnostic_tables()$quality_fields)
  # Ambiguous duplicate declarations are not used to interpret numeric fields.
  out$profiled[duplicated(out$field) | duplicated(out$field, fromLast = TRUE)] <- FALSE
  out
}

#' @keywords internal
#' @noRd
new_vcf_diagnostic_state <- function(parsed, mode, n_scanned) {
  state <- new.env(parent = emptyenv())
  state$mode <- mode
  state$fields <- vcf_diagnostic_fields(parsed$meta)
  state$fields$profiled <- state$fields$profiled & mode == "extended"
  state$alleles <- vector("list", n_scanned)
  state$genotypes <- vector("list", n_scanned)
  state$quality <- vector("list", n_scanned)
  state$unsupported <- vector("list", n_scanned)
  state$non_diploid_rows <- integer()
  state$flagged_samples <- logical(length(parsed$samples))
  state$n_flagged_cells <- 0
  state$n_inspected_cells <- 0
  state$n_unsupported_cells <- 0
  state$unsupported_samples <- logical(length(parsed$samples))
  state
}

#' @keywords internal
#' @noRd
collect_vcf_diagnostic_row <- function(state, row, tokens, subfields, fields, alleles, supported,
                                       valid, inconsistent) {
  genotypes <- unique(tokens)
  index <- match(tokens, genotypes)
  slots <- strsplit(genotypes, "[/|]")
  has_unsupported <- vapply(seq_along(genotypes), function(j) {
    if (!valid[match(genotypes[j], tokens)]) return(FALSE)
    called <- slots[[j]][slots[[j]] != "."]
    any(!supported[match(as.numeric(called), seq_along(alleles) - 1L)])
  }, logical(1))
  missing <- vapply(slots, function(x) length(x) > 0L && all(x == "."), logical(1))
  partial <- vapply(slots, function(x) any(x == "."), logical(1))
  categories <- ifelse(
    has_unsupported, "unsupported_allele",
    ifelse(missing, "missing", ifelse(partial, "partially_missing", "called_supported"))
  )
  categories[!valid[match(genotypes, tokens)]] <- "invalid"
  non_diploid <- state$target == "beagle" & !is.na(genotypes) & lengths(slots) != 2L
  categories[non_diploid & !has_unsupported & valid[match(genotypes, tokens)]] <- "non_diploid"
  unsupported <- has_unsupported[index] & !inconsistent
  flagged <- unsupported | !valid | inconsistent | non_diploid[index]
  state$flagged_samples <- state$flagged_samples | flagged
  state$unsupported_samples <- state$unsupported_samples | unsupported
  state$n_flagged_cells <- state$n_flagged_cells + sum(flagged)
  state$n_unsupported_cells <- state$n_unsupported_cells + sum(unsupported)
  state$n_inspected_cells <- state$n_inspected_cells + sum(!inconsistent)
  if (any(!is.na(genotypes) & lengths(slots) != 2L)) {
    state$non_diploid_rows <- c(state$non_diploid_rows, row)
  }

  interesting <- row %in% state$problem_rows || length(alleles) > 2L || any(!supported) || any(flagged) ||
    (state$target == "beagle" && any(lengths(slots) != 2L))
  if (!interesting) return(invisible(NULL))
  copies <- cells <- numeric(length(alleles))
  for (j in seq_along(genotypes)) {
    eligible <- index == j & valid & !inconsistent
    n <- sum(eligible)
    if (!n) next
    called <- slots[[j]][slots[[j]] != "."]
    allele_index <- match(as.numeric(called), seq_along(alleles) - 1L)
    copies <- copies + tabulate(allele_index, length(alleles)) * n
    cells[unique(allele_index)] <- cells[unique(allele_index)] + n
  }
  state$alleles[[row]] <- tibble::tibble(
    record = row, allele_index = seq_along(alleles) - 1L, allele = alleles,
    role = c("REF", rep("ALT", length(alleles) - 1L)), supported = supported,
    n_copies = copies, n_genotype_cells = as.integer(cells)
  )
  token_rows <- lapply(unique(alleles[!supported]), function(token) {
    selected <- if (is.na(token)) {
      which(is.na(alleles) & !supported)
    } else {
      which(!is.na(alleles) & alleles == token & !supported)
    }
    touches <- vapply(seq_along(genotypes), function(j) {
      if (!valid[match(genotypes[j], tokens)]) return(FALSE)
      any(suppressWarnings(as.numeric(slots[[j]])) %in% (selected - 1L))
    }, logical(1))
    tibble::tibble(
      record = row, allele = token, n_genotype_cells = sum(touches[index] & !inconsistent),
      n_copies = sum(copies[selected])
    )
  })
  if (length(token_rows)) state$unsupported[[row]] <- vcf_diagnostic_bind(token_rows, NULL)
  # Keep FORMAT-shape failures separate, even if their GT token looks valid.
  cell_categories <- categories[index]
  cell_categories[inconsistent] <- "uninterpretable"
  combinations <- unique(data.frame(genotype = tokens, category = cell_categories))
  state$genotypes[[row]] <- tibble::tibble(
    record = row, genotype = combinations$genotype, category = combinations$category,
    n = vapply(seq_len(nrow(combinations)), function(j) {
      same <- is.na(tokens) & is.na(combinations$genotype[j]) |
        !is.na(tokens) & !is.na(combinations$genotype[j]) & tokens == combinations$genotype[j]
      sum(same & cell_categories == combinations$category[j])
    }, integer(1))
  )
  if (state$mode != "extended") return(invisible(NULL))
  numeric_fields <- state$fields$field[state$fields$profiled]
  rows <- list()
  for (field in intersect(numeric_fields, fields)) {
    position <- match(field, fields)
    raw <- vapply(subfields, function(x) {
      if (length(x) >= position) x[position] else NA_character_
    }, character(1))
    number <- suppressWarnings(as.numeric(raw))
    for (category in unique(cell_categories[!inconsistent])) {
      selected <- cell_categories == category & !inconsistent
      absent <- is.na(raw[selected]) | raw[selected] %in% c(".", "")
      x <- number[selected]
      rows[[length(rows) + 1L]] <- tibble::tibble(
        record = row, field = field, category = category, n = sum(selected),
        n_missing = sum(absent), n_nonfinite = sum(!absent & (is.nan(x) | is.infinite(x))),
        n_nonnumeric = sum(!absent & is.na(x) & !is.nan(x)), values = list(x[is.finite(x)])
      )
    }
  }
  if (length(rows)) state$quality[[row]] <- vcf_diagnostic_bind(rows, NULL)
  invisible(NULL)
}

#' @keywords internal
#' @noRd
summarize_vcf_quality <- function(rows, tables) {
  if (!nrow(rows)) return(tables)
  by_record <- lapply(seq_len(nrow(rows)), function(i) {
    values <- rows$values[[i]]
    tibble::tibble(
      record = rows$record[i], field = rows$field[i], category = rows$category[i],
      n = rows$n[i], n_finite = length(values),
      median = if (length(values)) stats::median(values) else NA_real_
    )
  })
  tables$quality_by_record <- vcf_diagnostic_bind(by_record, tables$quality_by_record)
  groups <- split(seq_len(nrow(rows)), interaction(rows$field, rows$category, drop = TRUE))
  summary <- lapply(groups, function(indices) {
    values <- unlist(rows$values[indices], use.names = FALSE)
    quantiles <- if (length(values)) {
      stats::quantile(values, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
    } else {
      rep(NA_real_, 5L)
    }
    tibble::tibble(
      field = rows$field[indices[1]], category = rows$category[indices[1]], n = sum(rows$n[indices]),
      n_finite = length(values), n_missing = sum(rows$n_missing[indices]),
      n_nonfinite = sum(rows$n_nonfinite[indices]), n_nonnumeric = sum(rows$n_nonnumeric[indices]),
      min = quantiles[1], q25 = quantiles[2], median = quantiles[3], q75 = quantiles[4], max = quantiles[5]
    )
  })
  tables$quality_summary <- vcf_diagnostic_bind(summary, tables$quality_summary)
  priority <- match(tables$quality_summary$category, c(
    "unsupported_allele", "missing", "partially_missing", "invalid", "non_diploid", "called_supported"
  ))
  tables$quality_summary <- tables$quality_summary[order(priority, tables$quality_summary$field), ]
  comparisons <- lapply(unique(rows$field), function(field) {
    x <- tables$quality_by_record
    baseline <- x[x$field == field & x$category == "called_supported", ]
    unusual <- x[x$field == field & x$category == "unsupported_allele", ]
    matched <- match(unusual$record, baseline$record)
    a <- unusual$median
    b <- baseline$median[matched]
    good <- is.finite(a) & is.finite(b)
    a <- a[good]
    b <- b[good]
    ratio <- a[b > 0] / b[b > 0]
    tibble::tibble(
      field = field, n_comparable_records = length(a), n_lower = sum(a < b),
      n_not_comparable_records = sum(!good), n_ratio_records = length(ratio),
      n_equal = sum(a == b), n_higher = sum(a > b),
      median_difference = if (length(a)) stats::median(a - b) else NA_real_,
      median_positive_baseline_ratio = if (length(ratio)) stats::median(ratio) else NA_real_
    )
  })
  tables$quality_comparisons <- vcf_diagnostic_bind(comparisons, tables$quality_comparisons)
  tables
}

#' @keywords internal
#' @noRd
vcf_diagnostic_records <- function(parsed, n_scanned) {
  pieces <- parsed$body_split[seq_len(n_scanned)]
  column <- function(j) vapply(pieces, function(x) if (length(x) >= j) x[j] else NA_character_, character(1))
  tibble::tibble(record = seq_len(n_scanned), chrom = column(1), pos = column(2), marker_id = column(3))
}

#' @keywords internal
#' @noRd
vcf_diagnostic_issue_rows <- function(issue, records) {
  where <- issue$where
  rows <- where$rows
  if (!is.null(where$marker_id)) rows <- which(records$marker_id %in% where$marker_id)
  if (!is.null(where$chrom_pos)) {
    position <- suppressWarnings(as.numeric(records$pos))
    keys <- paste(records$chrom, format(position, scientific = FALSE, trim = TRUE), sep = ":")
    rows <- which(keys %in% where$chrom_pos)
  }
  if (identical(issue$code, "vcf_unsorted_records")) {
    position <- suppressWarnings(as.numeric(records$pos))
    selected <- which(
      records$chrom %in% where$chrom & is.finite(position) & position > 0 &
        abs(position - round(position)) <= .Machine$double.eps^0.5
    )
    breaks <- which(diff(position[selected]) < 0)
    rows <- selected[unique(c(breaks, breaks + 1L))]
  }
  sort(unique(as.integer(rows[rows %in% records$record])))
}

#' @keywords internal
#' @noRd
finish_vcf_diagnostics <- function(state, parsed, errors, warnings, n_scanned, partial) {
  records <- state$records
  if (is.null(records)) records <- vcf_diagnostic_records(parsed, n_scanned)
  tables <- vcf_diagnostic_tables()
  tables$quality_fields <- state$fields
  tables$allele_usage <- vcf_diagnostic_bind(state$alleles, tables$allele_usage)
  tables$genotype_counts <- vcf_diagnostic_bind(state$genotypes, tables$genotype_counts)
  counts <- tables$genotype_counts
  combinations <- unique(counts[, c("genotype", "category")])
  summaries <- lapply(seq_len(nrow(combinations)), function(i) {
    same <- is.na(counts$genotype) & is.na(combinations$genotype[i]) |
      !is.na(counts$genotype) & !is.na(combinations$genotype[i]) & counts$genotype == combinations$genotype[i]
    rows <- counts[same & counts$category == combinations$category[i], ]
    tibble::tibble(
      genotype = combinations$genotype[i], category = combinations$category[i],
      n_records = length(unique(rows$record)), n_genotype_cells = as.double(sum(rows$n))
    )
  })
  tables$genotype_summary <- vcf_diagnostic_bind(summaries, tables$genotype_summary)
  tokens <- vcf_diagnostic_bind(state$unsupported, tibble::tibble())
  if (nrow(tokens)) {
    tokens <- lapply(unique(tokens$allele), function(token) {
      x <- tokens[if (is.na(token)) is.na(tokens$allele) else tokens$allele %in% token, ]
      tibble::tibble(
        allele = token, n_records = length(unique(x$record)),
        n_genotype_cells = sum(x$n_genotype_cells), n_copies = sum(x$n_copies), records = list(x$record)
      )
    })
    tables$unsupported_alleles <- vcf_diagnostic_bind(tokens, tables$unsupported_alleles)
  }
  quality <- vcf_diagnostic_bind(state$quality, tibble::tibble())
  tables <- summarize_vcf_quality(quality, tables)
  tables$quality_fields$profiled <- tables$quality_fields$profiled &
    tables$quality_fields$field %in% tables$quality_summary$field
  issues <- c(errors, warnings)
  severity <- c(rep("error", length(errors)), rep("warning", length(warnings)))
  locations <- lapply(issues, vcf_diagnostic_issue_rows, records = records)
  record_codes <- vector("list", n_scanned)
  record_errors <- logical(n_scanned)
  for (i in seq_along(issues)) {
    for (row in locations[[i]]) record_codes[[row]] <- unique(c(record_codes[[row]], issues[[i]]$code))
    if (severity[i] == "error") record_errors[locations[[i]]] <- TRUE
  }
  selected <- which(lengths(record_codes) > 0L)
  record_issues <- records[selected, ]
  record_issues$codes <- lapply(record_codes[selected], sort)
  keys <- vapply(record_issues$codes, paste, collapse = "; ", FUN.VALUE = character(1))
  groups <- lapply(unique(keys), function(key) {
    rows <- record_issues$record[keys == key]
    tibble::tibble(
      codes = list(record_codes[[rows[1]]]), severity = if (any(record_errors[rows])) "error" else "warning",
      n_records = length(rows), records = list(rows)
    )
  })
  issue_groups <- vcf_diagnostic_bind(groups, tibble::tibble(
    codes = list(), severity = character(), n_records = integer(), records = list()
  ))
  issue_keys <- vapply(seq_along(issues), function(i) paste(severity[i], issues[[i]]$code), character(1))
  impact <- lapply(unique(issue_keys), function(key) {
    indices <- which(issue_keys == key)
    rows <- unique(unlist(locations[indices], use.names = FALSE))
    samples <- unique(unlist(lapply(issues[indices], function(x) x$where$sample), use.names = FALSE))
    scope <- if (length(rows)) "record" else if (length(samples)) "sample" else "file_or_unlocated"
    tibble::tibble(
      severity = severity[indices[1]], code = issues[[indices[1]]]$code,
      n_entries = length(indices), n_records = if (length(rows)) length(rows) else NA_integer_,
      n_samples = if (length(samples)) sum(parsed$samples %in% samples) else NA_integer_, scope = scope
    )
  })
  issue_impact <- vcf_diagnostic_bind(impact, tibble::tibble(
    severity = character(), code = character(), n_entries = integer(),
    n_records = integer(), n_samples = integer(), scope = character()
  ))
  codes <- unique(vapply(issues, `[[`, character(1), "code"))
  actions <- lapply(codes, function(code) {
    action <- if (code == "vcf_unsorted_records") {
      "Sort records after resolving invalid coordinates; do not change genotype calls."
    } else if (code %in% c("vcf_duplicated_chrom_pos", "vcf_duplicated_marker_ids")) {
      "Compare alleles and calls before choosing, merging or discarding duplicate records."
    } else if (code == "vcf_invalid_ref_alt") {
      "Confirm allele semantics. Review genotype indices and allele-numbered annotations before changing alleles."
    } else if (code == "vcf_multiallelic") {
      "Review target-specific representation; splitting or filtering requires genotype and annotation handling."
    } else if (code == "vcf_beagle_not_bgzf") {
      "Plain VCF is usable as gt=; reference-panel input needs verified BGZF compression and indexing."
    } else {
      "Inspect the reported locations and original source before changing data."
    }
    risk <- if (code == "vcf_unsorted_records") {
      "ordering_only"
    } else if (code == "vcf_beagle_not_bgzf") {
      "storage_only"
    } else {
      "review_required"
    }
    tibble::tibble(code = code, action = action, risk = risk)
  })
  header_available <- !is.na(parsed$chrom_header_line)
  coverage <- list(
    mode = state$mode,
    scope = if (!header_available) "header_unavailable" else if (partial) "scanned_prefix" else "full_file",
    n_records_total = length(parsed$body_lines), n_records_scanned = n_scanned,
    n_samples = length(parsed$samples),
    duplicate_coordinates_checked = isTRUE(state$coordinate_checks),
    sortedness_checked = isTRUE(state$coordinate_checks),
    quality_scope = if (state$mode == "extended") {
      "interpretable problem records; declared scalar numeric FORMAT fields"
    } else {
      "not_requested"
    }
  )
  structure(c(list(
    coverage = coverage,
    impact = list(
      n_records_with_findings = length(selected), n_records_with_errors = sum(record_errors),
      n_issue_entries_without_record_locations = sum(lengths(locations) == 0L),
      n_genotype_cells_inspected = state$n_inspected_cells,
      n_genotype_cells_not_inspected = as.double(n_scanned) * length(parsed$samples) - state$n_inspected_cells,
      n_genotype_cells_flagged = state$n_flagged_cells, n_samples_with_flagged_calls = sum(state$flagged_samples),
      n_genotype_cells_using_unsupported_alleles = state$n_unsupported_cells,
      n_samples_using_unsupported_alleles = sum(state$unsupported_samples)
    ),
    issue_impact = issue_impact, issue_groups = issue_groups, record_issues = record_issues,
    records = records,
    actions = vcf_diagnostic_bind(
      actions, tibble::tibble(code = character(), action = character(), risk = character())
    ),
    interpretation = paste(
      "Unsupported allele meaning is unresolved. Numeric FORMAT summaries are descriptive,",
      "not probabilities, universal quality thresholds, or a missing-call/deletion diagnosis."
    )
  ), tables), class = c("galick_vcf_diagnostics", "list"))
}
