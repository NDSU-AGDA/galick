diagnostic_vcf <- function() {
  header <- c(
    "##fileformat=VCFv4.2",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=Score,Number=1,Type=Float,Description="Custom score, not a probability">',
    '##FORMAT=<ID=Depth,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
    '##FORMAT=<ID=Label,Number=1,Type=String,Description="Text">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\te"
  )
  rows <- list(
    c("1", "100", "m1", "A", "G,-", ".", "PASS", ".", "GT:Score:Depth:AD:Label:Extra",
      "2/2:0.8:2:0,0,2:x:4", "0/0:0.9:20:20,0,0:x:4", "./.:.:.:.:x:4",
      "0/2:Inf:bad:1,0,1:x:4", "2|.:NaN:0:0,0,1:x:4"),
    c("1", "100", "m2", "C", "T,-", ".", "PASS", ".", "Depth:GT:Score:AD:Label:Extra",
      "2:2/2:0.7:0,0,2:x:4", "30:1/1:0.6:0,2,0:x:4", ".:./.:oops:.:x:4",
      "5:0/1:-0.1:2,3,0:x:4", "9:1|1:0.6:0,2,0:x:4"),
    c("1", "90", "m3", "G", "A,-", ".", "PASS", ".", "GT:Score:Depth:AD:Label:Extra",
      rep("0/0:0.9:20:20,0,0:x:4", 5)),
    c("2", "50", "m4", "T", "C", ".", "PASS", ".", "GT:Score:Depth:AD:Label:Extra",
      rep("0/0:0.9:20:20,0:x:4", 5))
  )
  path <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  writeLines(c(header, vapply(rows, paste, collapse = "\t", FUN.VALUE = character(1))), path)
  path
}

capture_vcf_report <- function(x, full = FALSE) {
  capture.output(invisible(capture.output(if (full) summary(x) else print(x))), type = "message")
}

test_that("diagnostics distinguish issue entries, unique records and dependent calls", {
  result <- validate_vcf(diagnostic_vcf(), target = "beagle", progress = FALSE)
  d <- result$summary$diagnostics
  storage <- d$actions$action[d$actions$code == "vcf_beagle_not_bgzf"]
  expect_match(storage, "the runner does not require an input index", fixed = TRUE)
  expect_s3_class(result, "galick_result")
  expect_s3_class(d, "galick_vcf_diagnostics")
  expect_equal(nrow(result$errors), 3)
  expect_equal(d$impact$n_records_with_errors, 3)
  expect_equal(d$impact$n_records_with_findings, 3)
  expect_equal(d$impact$n_genotype_cells_using_unsupported_alleles, 4)
  expect_equal(d$impact$n_samples_using_unsupported_alleles, 3)
  expect_equal(d$impact$n_genotype_cells_flagged, 4)
  expect_equal(d$impact$n_genotype_cells_inspected, 20)
  expect_equal(d$unsupported_alleles$n_records, 3)
  expect_equal(d$unsupported_alleles$n_genotype_cells, 4)
  expect_equal(d$unsupported_alleles$n_copies, 6)
  calls <- d$genotype_summary[d$genotype_summary$category == "unsupported_allele", ]
  expect_equal(sum(calls$n_genotype_cells), 4)
  expect_equal(calls$n_genotype_cells[calls$genotype == "2/2"], 2)
  duplicate <- d$issue_impact[d$issue_impact$code == "vcf_duplicated_chrom_pos", ]
  expect_equal(duplicate$n_entries, 1)
  expect_equal(duplicate$n_records, 2)
  expect_setequal(unique(unlist(d$issue_groups$records)), 1:3)
  expect_equal(sum(d$issue_groups$n_records), 3)
  usage <- d$allele_usage[d$allele_usage$record == 3 & d$allele_usage$allele == "-", ]
  expect_equal(usage$n_genotype_cells, 0)
  expect_equal(d$allele_usage$allele[d$allele_usage$record == 1], c("A", "G", "-"))
  expect_length(d$quality_summary$field, 0)
})

test_that("extended profiling follows declarations and reordered FORMAT fields", {
  path <- diagnostic_vcf()
  before <- digest::digest(file = path, algo = "sha256")
  standard <- validate_vcf(path, target = "beagle", progress = FALSE)
  result <- validate_vcf(path, target = "beagle", progress = FALSE, diagnostics = "extended")
  d <- result$summary$diagnostics
  expect_identical(result$ok, standard$ok)
  expect_identical(result$errors, standard$errors)
  expect_identical(result$warnings, standard$warnings)
  expect_identical(d$allele_usage, standard$summary$diagnostics$allele_usage)
  expect_identical(digest::digest(file = path, algo = "sha256"), before)
  expect_setequal(d$quality_summary$field, c("Score", "Depth"))
  expect_equal(d$quality_fields$description[d$quality_fields$field == "Score"], "Custom score, not a probability")
  selected <- d$quality_summary$category == "unsupported_allele"
  score <- d$quality_summary[d$quality_summary$field == "Score" & selected, ]
  expect_equal(score$n, 4)
  expect_equal(score$n_finite, 2)
  expect_equal(score$median, 0.75)
  expect_equal(score$n_nonfinite, 2)
  depth <- d$quality_summary[d$quality_summary$field == "Depth" & selected, ]
  expect_equal(depth$n_nonnumeric, 1)
  expect_equal(depth$median, 2)
  missing <- d$quality_summary[d$quality_summary$field == "Score" & d$quality_summary$category == "missing", ]
  expect_equal(missing$n_missing, 1)
  expect_equal(missing$n_nonnumeric, 1)
  expect_true(is.na(missing$median))
  comparison <- d$quality_comparisons[d$quality_comparisons$field == "Score", ]
  expect_equal(comparison$n_comparable_records, 2)
  expect_equal(comparison$n_lower, 1)
  expect_equal(comparison$n_higher, 1)
  expect_equal(comparison$n_not_comparable_records, 0)
  expect_equal(comparison$n_ratio_records, 2)
  expect_equal(comparison$median_difference, 0, tolerance = 1e-10)
  expect_true(all(c("AD", "Label") %in% d$quality_fields$field[!d$quality_fields$profiled]))
  expect_false("Extra" %in% d$quality_summary$field)
  expect_match(d$interpretation, "unresolved")
})

test_that("duplicate unsupported allele indices do not double count genotype cells", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "1\t1\tm1\tA\tG,-,-\t.\tPASS\t.\tGT\t2|3"
  ), path)
  d <- validate_vcf(path, progress = FALSE)$summary$diagnostics
  expect_equal(d$unsupported_alleles$n_genotype_cells, 1)
  expect_equal(d$unsupported_alleles$n_copies, 2)
  expect_equal(d$impact$n_genotype_cells_using_unsupported_alleles, 1)
})

test_that("partial reports never imply full clearance or inspect beyond the prefix", {
  x <- validate_vcf(
    diagnostic_vcf(), target = "beagle", max_records_to_scan = 1,
    progress = FALSE, diagnostics = "extended"
  )
  d <- x$summary$diagnostics
  expect_equal(d$coverage$scope, "scanned_prefix")
  expect_equal(d$coverage$n_records_total, 4)
  expect_equal(d$coverage$n_records_scanned, 1)
  expect_false(d$coverage$sortedness_checked)
  expect_false(d$coverage$duplicate_coordinates_checked)
  expect_equal(d$records$record, 1)
  expect_true(all(d$allele_usage$record == 1))
  expect_equal(d$unsupported_alleles$n_genotype_cells, 3)
  expect_match(paste(capture_vcf_report(x), collapse = "\n"), "PARTIAL scan")
})

test_that("missing declarations, GT, samples, headers and malformed rows remain inspectable", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "1\t1\tm1\tA\tG\t.\tPASS\t.\tGT\t0/0"), path)
  x <- validate_vcf(path, progress = FALSE, diagnostics = "extended")
  expect_false(x$ok)
  expect_equal(x$summary$diagnostics$coverage$scope, "header_unavailable")
  writeLines(c("##fileformat=VCFv4.2", "#CHROM\tPOS\tID", "1\t1\tm1"), path)
  expect_false(validate_vcf(path, progress = FALSE)$ok)
  writeLines(c(
    "##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "1\t1\tm1\tA\tG,-\t.\tPASS\t.\tDepth\t10",
    "1\t2\tm2\tA\tG,-\t.\tPASS\t.\tGT:Depth\t2/2",
    "1\t3\tm3\tA\tG\t.\tPASS"
  ), path)
  x <- validate_vcf(path, progress = FALSE, diagnostics = "extended")
  expect_false(x$ok)
  expect_equal(x$summary$diagnostics$impact$n_genotype_cells_not_inspected, 3)
  expect_equal(x$summary$diagnostics$genotype_counts$category, "uninterpretable")
  expect_equal(nrow(x$summary$diagnostics$quality_summary), 0)
  writeLines(c(
    "##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    "1\t1\tm1\tA\tG,-\t.\tPASS\t."
  ), path)
  x <- validate_vcf(path, progress = FALSE, diagnostics = "extended")
  expect_false(x$ok)
  expect_equal(x$summary$diagnostics$impact$n_genotype_cells_inspected, 0)
})

test_that("reports keep biological interpretation separate from observed numeric values", {
  x <- validate_vcf(diagnostic_vcf(), target = "beagle", progress = FALSE, diagnostics = "extended")
  expect_match(paste(capture_vcf_report(x), collapse = "\n"), "Entry counts are not counts")
  rendered <- character()
  messages <- capture.output(rendered <- capture.output(print(x)), type = "message")
  expect_true(any(grepl("0=A, 1=G, 2=-", c(rendered, messages), fixed = TRUE)))
  detailed <- paste(capture_vcf_report(x, full = TRUE), collapse = "\n")
  expect_match(detailed, "Within-record comparisons")
  expect_match(detailed, "meaning is unresolved")
  expect_match(detailed, "Action assessment")
  expect_error(validate_vcf(diagnostic_vcf(), diagnostics = "automatic_repair"), "arg")
})

test_that("multi-digit and malformed genotypes retain existing validation decisions", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    paste0("1\t1\tm1\tA\t", paste(rep("C", 11), collapse = ","), "\t.\tPASS\t.\tGT\t10|11\t12/0"),
    "1\t2\tm2\tA\tG\t.\tPASS\t.\tGT\t0/1/1\ta/b"
  ), path)
  x <- validate_vcf(path, target = "beagle", progress = FALSE)
  expect_equal(x$errors$where[[which(x$errors$code == "vcf_invalid_gt")]]$rows, 1:2)
  expect_equal(x$errors$where[[which(x$errors$code == "vcf_beagle_non_diploid")]]$rows, 2)
  expect_equal(x$summary$diagnostics$impact$n_genotype_cells_flagged, 3)
})

test_that("ambiguous and vector declarations do not acquire invented numeric interpretations", {
  meta <- c(
    '##FORMAT=<ID=Q,Number=1,Type=Float,Description="First">',
    '##FORMAT=<ID=Q,Number=1,Type=Float,Description="Second">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Vector">'
  )
  fields <- vcf_diagnostic_fields(meta)
  expect_false(any(fields$profiled))
})

test_that("unsupported REF and symbolic ALT dependencies are generic and union-counted", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "1\t1\tm1\tZ\t<DEL>,N\t.\tPASS\t.\tGT\t0/1"
  ), path)
  x <- validate_vcf(path, progress = FALSE)
  d <- x$summary$diagnostics
  expect_equal(d$impact$n_genotype_cells_using_unsupported_alleles, 1)
  expect_setequal(d$unsupported_alleles$allele, c("Z", "<DEL>"))
  expect_equal(d$unsupported_alleles$n_genotype_cells, c(1, 1))
  expect_equal(d$issue_impact$n_entries[d$issue_impact$code == "vcf_invalid_ref_alt"], 2)
  expect_equal(d$issue_impact$n_records[d$issue_impact$code == "vcf_invalid_ref_alt"], 1)
  expect_equal(d$impact$n_records_with_findings, 1)
})

test_that("numeric fields on uninterpretable FORMAT shapes are not profiled", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", '##FORMAT=<ID=Q,Number=1,Type=Float,Description="Custom">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "1\t1\tm1\tA\tG,-\t.\tPASS\t.\tGT:Q\t2/2:0.9:unexpected"
  ), path)
  x <- validate_vcf(path, progress = FALSE, diagnostics = "extended")
  expect_equal(nrow(x$summary$diagnostics$quality_summary), 0)
  expect_equal(x$summary$diagnostics$genotype_counts$category, "uninterpretable")
})

test_that("diagnostic reports have stable compact and detailed renderings", {
  testthat::local_reproducible_output()
  x <- validate_vcf(diagnostic_vcf(), target = "beagle", progress = FALSE, diagnostics = "extended")
  expect_snapshot(print(x))
  expect_snapshot(summary(x))
})

test_that("quality comparisons identify absent and non-positive baselines", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2", '##FORMAT=<ID=Q,Number=1,Type=Float,Description="Custom">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "1\t1\tm1\tA\tG,-\t.\tPASS\t.\tGT:Q\t2/2:0.8\t2/2:0.9",
    "1\t2\tm2\tA\tG,-\t.\tPASS\t.\tGT:Q\t2/2:0.8\t0/0:0",
    "1\t3\tm3\tA\tG,-\t.\tPASS\t.\tGT:Q\t2/2:0.8\t0/0:-1"
  ), path)
  d <- validate_vcf(path, progress = FALSE, diagnostics = "extended")$summary$diagnostics
  expect_equal(d$quality_comparisons$n_comparable_records, 2)
  expect_equal(d$quality_comparisons$n_not_comparable_records, 1)
  expect_equal(d$quality_comparisons$n_ratio_records, 0)
  expect_true(is.na(d$quality_comparisons$median_positive_baseline_ratio))
  expect_equal(d$quality_comparisons$median_difference, 1.3)
})
