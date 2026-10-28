#' @keywords internal
#' @noRd
print_vcf_validation <- function(x, full = FALSE) {
  d <- x$summary$diagnostics
  coverage <- d$coverage
  cli::cli_h1("galick_result: {if (x$ok) 'OK' else 'NOT OK'}")
  cli::cli_text("Validation entries: {nrow(x$errors)} errors; {nrow(x$warnings)} warnings.")
  cli::cli_text("Entry counts are not counts of affected records or genotype calls.")
  scope <- switch(coverage$scope,
    full_file = "Full scan", scanned_prefix = "PARTIAL scan", header_unavailable = "Header unavailable"
  )
  cli::cli_text(sprintf(
    "%s: %d of %d records; %d samples.", scope,
    coverage$n_records_scanned, coverage$n_records_total, coverage$n_samples
  ))
  if (coverage$scope != "full_file") {
    cli::cli_alert_warning("No full-file clearance: coordinate duplicates and sortedness were not checked.")
  }
  cli::cli_h3("Impact in scanned records")
  cli::cli_text(
    "Records with located findings: {d$impact$n_records_with_findings}; with errors: {d$impact$n_records_with_errors}."
  )
  cli::cli_text("File/sample or unlocated entries: {d$impact$n_issue_entries_without_record_locations}.")
  cli::cli_text(paste0(
    "Flagged genotype cells: {d$impact$n_genotype_cells_flagged}; ",
    "affected sample columns: {d$impact$n_samples_with_flagged_calls}."
  ))
  if (d$impact$n_genotype_cells_not_inspected > 0) {
    cli::cli_alert_warning("{d$impact$n_genotype_cells_not_inspected} genotype cells could not be inspected.")
  }
  cli::cli_text("Cell counts cover genotype/FORMAT findings, not every call in a structurally flagged record.")
  if (nrow(d$issue_impact)) {
    cli::cli_h3("Rule impact")
    print(d$issue_impact, n = if (full) Inf else 6L)
  }
  if (nrow(d$issue_groups)) {
    cli::cli_h3("Overlapping record findings")
    groups <- tibble::tibble(
      n_records = d$issue_groups$n_records,
      codes = vapply(d$issue_groups$codes, paste, collapse = "; ", FUN.VALUE = character(1))
    )
    print(groups, n = if (full) Inf else 6L)
  }
  if (nrow(d$unsupported_alleles)) {
    cli::cli_h3("Unsupported allele dependencies")
    tokens <- d$unsupported_alleles
    for (i in seq_len(min(nrow(tokens), if (full) nrow(tokens) else 3L))) {
      cli::cli_text(paste0(
        "Allele {.val {tokens$allele[i]}}: {tokens$n_records[i]} records; ",
        "{tokens$n_genotype_cells[i]} dependent genotype cells."
      ))
      example <- d$allele_usage[d$allele_usage$record == tokens$records[[i]][1], ]
      labels <- paste0(example$allele_index, "=", example$allele)
      if (length(labels) > 6L) labels <- c(utils::head(labels, 6L), "...")
      cli::cli_verbatim(sprintf(
        "Example allele indices (record %d): %s", example$record[1], paste(labels, collapse = ", ")
      ))
    }
    cli::cli_alert_warning("{d$interpretation}")
    cli::cli_text("Removing an allele requires review of dependent calls and allele-numbered INFO/FORMAT annotations.")
    genotypes <- d$genotype_summary[d$genotype_summary$category == "unsupported_allele", ]
    if (nrow(genotypes)) print(genotypes, n = if (full) Inf else 3L)
  }
  if (nrow(d$quality_summary)) {
    cli::cli_h3("Numeric FORMAT evidence on problem records")
    fields <- c("field", "category", "n_finite", "median", "n_missing", "n_nonfinite", "n_nonnumeric")
    print(d$quality_summary[, fields], n = if (full) Inf else 6L)
    cli::cli_text("These pooled summaries are descriptive; within-record comparisons are reported separately.")
  } else if (coverage$mode == "extended") {
    cli::cli_text("No eligible numeric FORMAT evidence was available on interpretable problem records.")
  } else {
    cli::cli_text('Quality profiling not requested; use diagnostics = "extended" for numeric FORMAT evidence.')
  }
  if (full) {
    if (nrow(d$quality_fields)) {
      cli::cli_h3("FORMAT declarations")
      print(d$quality_fields, n = Inf)
      cli::cli_text("Only unambiguous Number=1 Integer/Float declarations are profiled in extended mode.")
    }
    if (nrow(d$quality_comparisons)) {
      cli::cli_h3("Within-record comparisons")
      cli::cli_text("Unsupported-allele calls versus fully called supported alleles at the same record.")
      columns <- c("field", "n_comparable_records", "n_not_comparable_records", "n_lower", "n_equal", "n_higher")
      print(d$quality_comparisons[, columns], n = Inf)
      for (i in seq_len(nrow(d$quality_comparisons))) {
        row <- d$quality_comparisons[i, ]
        row$median_difference <- signif(row$median_difference, 4)
        row$median_positive_baseline_ratio <- signif(row$median_positive_baseline_ratio, 4)
        cli::cli_text(paste0(
          "{row$field}: median difference {row$median_difference}; ",
          "median ratio {row$median_positive_baseline_ratio} from {row$n_ratio_records} positive baselines."
        ))
      }
      cli::cli_text(paste(
        "Ratios use positive baseline medians only;",
        "numeric differences do not establish biological meaning."
      ))
    }
    if (nrow(d$genotype_counts)) {
      cli::cli_h3("Genotype totals on problem records")
      print(d$genotype_summary, n = Inf)
      cli::cli_h3("Genotype categories on problem records")
      print(d$genotype_counts, n = 20L)
    }
    if (nrow(d$allele_usage)) {
      cli::cli_h3("Allele-index usage on problem records")
      print(d$allele_usage, n = 20L)
    }
    if (nrow(d$record_issues)) {
      cli::cli_h3("Record locations")
      print(d$record_issues[, c("record", "chrom", "pos", "marker_id")], n = 10L)
    }
    if (nrow(x$summary$chromosome_summary)) {
      cli::cli_h3("Chromosomes")
      print(x$summary$chromosome_summary, n = Inf)
    }
  }
  if (nrow(d$actions)) {
    cli::cli_h3("Action assessment")
    for (i in seq_len(min(nrow(d$actions), if (full) nrow(d$actions) else 3L))) {
      cli::cli_text("{d$actions$code[i]} [{d$actions$risk[i]}]: {d$actions$action[i]}")
    }
  }
  invisible(x)
}
