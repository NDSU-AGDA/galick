make_map <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}

test_that("validate_marker_map happy path", {
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3"),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 50)
  )
  result <- validate_marker_map(map)
  expect_s3_class(result, "galick_result")
  expect_true(result$ok)
  expect_equal(result$summary$n_markers_in_map, 3L)
  expect_equal(result$summary$inferred_chrom_style, "chr_prefixed")
  expect_equal(nrow(result$summary$chromosome_summary), 2L)
})

test_that("validate_marker_map aborts on NULL", {
  expect_error(
    validate_marker_map(NULL),
    class = "galick_error_marker_map_invalid_input"
  )
})

test_that("validate_marker_map aborts on non-data.frame", {
  expect_error(
    validate_marker_map(list(a = 1)),
    class = "galick_error_marker_map_invalid_input"
  )
})

test_that("validate_marker_map aborts on empty data.frame", {
  expect_error(
    validate_marker_map(data.frame(marker_id = character(), chrom = character(), pos = integer())),
    class = "galick_error_marker_map_invalid_input"
  )
})

test_that("validate_marker_map aborts on missing required columns", {
  expect_error(
    validate_marker_map(make_map(marker_id = "rs1", chrom = "chr1")),
    class = "galick_error_marker_map_missing_required_columns"
  )
  expect_error(
    validate_marker_map(make_map(marker_id = "rs1", pos = 100)),
    class = "galick_error_marker_map_missing_required_columns"
  )
  expect_error(
    validate_marker_map(make_map(chrom = "chr1", pos = 100)),
    class = "galick_error_marker_map_missing_required_columns"
  )
})

test_that("validate_marker_map aborts when both genotype_table and marker_ids passed", {
  map <- make_map(marker_id = "rs1", chrom = "chr1", pos = 100)
  expect_error(
    validate_marker_map(
      map,
      genotype_table = make_map(sample_id = "s1", rs1 = "A"),
      marker_ids = c("rs1")
    ),
    class = "galick_error_marker_map_invalid_input"
  )
})

test_that("validate_marker_map flags NA marker_id as error", {
  map <- make_map(
    marker_id = c("rs1", NA, "rs3"),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 50)
  )
  result <- validate_marker_map(map)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_invalid_marker_id"))
})

test_that("validate_marker_map flags empty / whitespace marker_id", {
  map <- make_map(
    marker_id = c("rs1", "", "   "),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 50)
  )
  result <- validate_marker_map(map)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_invalid_marker_id"))
})

test_that("validate_marker_map flags invalid pos values", {
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3", "rs4", "rs5"),
    chrom = c("chr1", "chr1", "chr1", "chr1", "chr1"),
    pos = c(100, NA, -1, 0, 1.5)
  )
  result <- validate_marker_map(map)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_invalid_pos"))
})

test_that("validate_marker_map flags invalid chrom values", {
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3"),
    chrom = c("chr1", NA, ""),
    pos = c(100, 200, 300)
  )
  result <- validate_marker_map(map)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_invalid_chrom"))
})

test_that("validate_marker_map flags duplicate marker_ids", {
  map <- make_map(
    marker_id = c("rs1", "rs1", "rs2"),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 50)
  )
  result <- validate_marker_map(map)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_duplicated_marker_ids"))
})

test_that("validate_marker_map flags duplicate (chrom, pos) under as_is", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100, 100)
  )
  result <- validate_marker_map(map, chrom_style = "as_is")
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_duplicated_chrom_pos"))
})

test_that("validate_marker_map under numeric normalization catches chr1+1 duplicates", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "1"),
    pos = c(100, 100)
  )
  result <- validate_marker_map(map, chrom_style = "numeric")
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_duplicated_chrom_pos"))
})

test_that("validate_marker_map under as_is keeps chr1 and 1 separate", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "1"),
    pos = c(100, 100)
  )
  result <- validate_marker_map(map, chrom_style = "as_is")
  expect_false(any(result$errors$code == "marker_map_duplicated_chrom_pos"))
})

test_that("validate_marker_map flags unsorted-within-chrom as warning", {
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3"),
    chrom = c("chr1", "chr1", "chr1"),
    pos = c(200, 100, 300)
  )
  result <- validate_marker_map(map)
  expect_true(result$ok)
  expect_true(any(result$warnings$code == "marker_map_unsorted_within_chrom"))
  expect_equal(result$summary$n_unsorted_within_chrom, 1L)
})

test_that("validate_marker_map chromosome_summary contains per-chrom counts", {
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3", "rs4"),
    chrom = c("chr1", "chr1", "chr2", "chr2"),
    pos = c(100, 200, 50, 25)
  )
  result <- validate_marker_map(map)
  cs <- result$summary$chromosome_summary
  expect_equal(nrow(cs), 2L)
  expect_setequal(cs$chrom, c("chr1", "chr2"))
  expect_equal(cs$n_unsorted[cs$chrom == "chr1"], 0L)
  expect_equal(cs$n_unsorted[cs$chrom == "chr2"], 1L)
})

test_that("validate_marker_map warns on mixed chrom style under auto", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "2"),
    pos = c(100, 200)
  )
  result <- validate_marker_map(map, chrom_style = "auto")
  expect_true(any(result$warnings$code == "marker_map_mixed_chrom_style"))
})

test_that("validate_marker_map warns on unknown chrom style under auto", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("1A", "2B"),
    pos = c(100, 200)
  )
  result <- validate_marker_map(map, chrom_style = "auto")
  expect_true(any(result$warnings$code == "marker_map_unknown_chrom_style"))
})

test_that("validate_marker_map require_alleles = TRUE accepts valid ref/alt", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100, 200),
    ref = c("A", "C"),
    alt = c("T", "G")
  )
  result <- validate_marker_map(map, require_alleles = TRUE)
  expect_true(result$ok)
})

test_that("validate_marker_map require_alleles = TRUE rejects missing columns", {
  map <- make_map(
    marker_id = c("rs1"),
    chrom = c("chr1"),
    pos = c(100)
  )
  result <- validate_marker_map(map, require_alleles = TRUE)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_missing_alleles"))
})

test_that("validate_marker_map require_alleles = TRUE rejects empty ref/alt", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100, 200),
    ref = c("A", ""),
    alt = c("T", "G")
  )
  result <- validate_marker_map(map, require_alleles = TRUE)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_missing_alleles"))
})

test_that("validate_marker_map require_alleles = TRUE rejects ref == alt", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100, 200),
    ref = c("A", "C"),
    alt = c("T", "C")
  )
  result <- validate_marker_map(map, require_alleles = TRUE)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "marker_map_missing_alleles"))
  expect_match(
    paste(result$errors$message, collapse = " "),
    "ref == alt"
  )
})

test_that("validate_marker_map require_alleles = FALSE skips allele checks", {
  map <- make_map(
    marker_id = c("rs1"),
    chrom = c("chr1"),
    pos = c(100)
  )
  result <- validate_marker_map(map, require_alleles = FALSE)
  expect_true(result$ok)
})

test_that("validate_marker_map cross-checks markers_ids", {
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3"),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 50)
  )
  result <- validate_marker_map(map, marker_ids = c("rs1", "rs2", "rs99"))
  expect_true(any(result$warnings$code == "marker_map_markers_in_table_not_in_map"))
  expect_true(any(result$warnings$code == "marker_map_markers_in_map_not_in_table"))
  expect_setequal(result$summary$markers_in_table_not_in_map, "rs99")
  expect_setequal(result$summary$markers_in_map_not_in_table, "rs3")
})

test_that("validate_marker_map cross-checks via genotype_table", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100, 200)
  )
  gt <- data.frame(
    sample_id = c("s1", "s2"),
    rs1 = c("A", "B"),
    rs2 = c("H", "A"),
    rs3 = c("A", "A"),
    stringsAsFactors = FALSE
  )
  result <- validate_marker_map(map, genotype_table = gt)
  expect_true(any(result$warnings$code == "marker_map_markers_in_table_not_in_map"))
  expect_setequal(result$summary$markers_in_table_not_in_map, "rs3")
})

test_that("validate_marker_map cross-check exact match has no warnings", {
  map <- make_map(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100, 200)
  )
  result <- validate_marker_map(map, marker_ids = c("rs1", "rs2"))
  expect_false(any(result$warnings$code == "marker_map_markers_in_table_not_in_map"))
  expect_false(any(result$warnings$code == "marker_map_markers_in_map_not_in_table"))
})

test_that("validate_marker_map snapshot for ok case", {
  testthat::local_reproducible_output()
  map <- make_map(
    marker_id = c("rs1", "rs2", "rs3"),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 200, 50)
  )
  result <- validate_marker_map(map)
  expect_snapshot(print(result))
})

test_that("validate_marker_map snapshot for failing case", {
  testthat::local_reproducible_output()
  map <- make_map(
    marker_id = c("rs1", "rs1", NA),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100, 100, -5)
  )
  result <- validate_marker_map(map)
  expect_snapshot(print(result))
})
