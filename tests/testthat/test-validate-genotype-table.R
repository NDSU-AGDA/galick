make_data <- function(...) {
  data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
}

test_that("validate_genotype_table happy path with AHB encoding", {
  data <- make_data(
    sample_id = c("s1", "s2", "s3"),
    rs1 = c("A", "H", "B"),
    rs2 = c("A", "A", "H"),
    rs3 = c("B", "B", "H")
  )
  result <- validate_genotype_table(data)
  expect_s3_class(result, "beagler_result")
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "AHB")
  expect_equal(result$summary$n_samples, 3L)
  expect_equal(result$summary$n_markers, 3L)
})

test_that("validate_genotype_table happy path with 012 encoding", {
  data <- make_data(
    sample_id = c("s1", "s2", "s3"),
    rs1 = c("0", "1", "2"),
    rs2 = c("0", "1", "1")
  )
  result <- validate_genotype_table(data)
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "012")
})

test_that("validate_genotype_table happy path with VCF encoding", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("0/0", "1/1"),
    rs2 = c("0/1", "0/0")
  )
  result <- validate_genotype_table(data)
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "VCF")
})

test_that("validate_genotype_table happy path with AABB encoding", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("AA", "BB"),
    rs2 = c("AB", "AA")
  )
  result <- validate_genotype_table(data)
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "AABB")
})

test_that("validate_genotype_table happy path with REFHET encoding", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("REF", "HET"),
    rs2 = c("ALT", "REF")
  )
  result <- validate_genotype_table(data)
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "REFHET")
})

test_that("validate_genotype_table happy path with IUPAC encoding", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("A", "T"),
    rs2 = c("C", "G")
  )
  result <- validate_genotype_table(data)
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "IUPAC")
})

test_that("validate_genotype_table auto-detects sample column", {
  data <- make_data(
    sample_id = c("plant_1", "plant_2", "plant_3"),
    rs1 = c("A", "H", "B")
  )
  result <- validate_genotype_table(data)
  expect_equal(result$summary$n_samples, 3L)
  expect_equal(result$summary$n_markers, 1L)
})

test_that("validate_genotype_table auto-detect skips genotype-shaped columns", {
  # 'maybe_sample' looks like genotype tokens; 'pi_id' is the real sample col
  data <- make_data(
    maybe_sample = c("A", "H", "B"),
    pi_id = c("plant_1", "plant_2", "plant_3"),
    rs1 = c("A", "B", "H")
  )
  result <- validate_genotype_table(data)
  expect_equal(result$summary$n_markers, 2L)
})

test_that("validate_genotype_table accepts explicit sample_col", {
  data <- make_data(
    geno_a = c("A", "B"),
    pi_id = c("plant1", "plant2"),
    geno_b = c("H", "A")
  )
  result <- validate_genotype_table(data, sample_col = "pi_id")
  expect_equal(result$summary$n_samples, 2L)
  expect_equal(result$summary$n_markers, 2L)
})

test_that("validate_genotype_table aborts on NULL data", {
  expect_error(
    validate_genotype_table(NULL),
    class = "beagler_error_invalid_input"
  )
})

test_that("validate_genotype_table aborts on non-data.frame data", {
  expect_error(
    validate_genotype_table(list(a = 1, b = 2)),
    class = "beagler_error_invalid_input"
  )
})

test_that("validate_genotype_table aborts on empty / single-column data", {
  expect_error(
    validate_genotype_table(data.frame()),
    class = "beagler_error_invalid_input"
  )
  expect_error(
    validate_genotype_table(data.frame(sample_id = "s1")),
    class = "beagler_error_invalid_input"
  )
})

test_that("validate_genotype_table aborts on missing sample_col arg", {
  data <- make_data(sample_id = c("s1", "s2"), rs1 = c("A", "B"))
  expect_error(
    validate_genotype_table(data, sample_col = "not_there"),
    class = "beagler_error_invalid_sample_col"
  )
})

test_that("validate_genotype_table aborts when auto-detect cannot find sample_col", {
  # All character columns look like genotype tokens
  data <- make_data(
    col_a = c("A", "H", "B"),
    col_b = c("A", "A", "H")
  )
  expect_error(
    validate_genotype_table(data),
    class = "beagler_error_invalid_sample_col"
  )
})

test_that("validate_genotype_table aborts on unsupported ploidy", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("A", "B")
  )
  expect_error(
    validate_genotype_table(data, ploidy = 4L),
    class = "beagler_error_unsupported_ploidy"
  )
})

test_that("validate_genotype_table flags duplicate samples as errors", {
  data <- make_data(
    sample_id = c("s1", "s1", "s2"),
    rs1 = c("A", "B", "H"),
    rs2 = c("A", "A", "H")
  )
  result <- validate_genotype_table(data, sample_col = "sample_id")
  expect_false(result$ok)
  expect_true(any(result$errors$code == "genotype_table_duplicated_samples"))
})

test_that("validate_genotype_table flags duplicate markers as errors", {
  data <- structure(
    list(
      sample_id = c("s1", "s2"),
      rs1 = c("A", "B"),
      rs1 = c("H", "A")
    ),
    class = "data.frame",
    names = c("sample_id", "rs1", "rs1"),
    row.names = c(NA, -2L)
  )
  result <- validate_genotype_table(data, sample_col = "sample_id")
  expect_false(result$ok)
  expect_true(any(result$errors$code == "genotype_table_duplicated_markers"))
})

test_that("validate_genotype_table flags monomorphic markers as warnings", {
  data <- make_data(
    sample_id = c("s1", "s2", "s3"),
    rs_mono = c("A", "A", "A"),
    rs_var = c("A", "B", "H")
  )
  result <- validate_genotype_table(data)
  expect_true(result$ok)
  expect_true(any(result$warnings$code == "genotype_table_monomorphic_markers"))
  expect_equal(result$summary$n_monomorphic_markers, 1L)
})

test_that("validate_genotype_table flags high sample missingness as warning", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("A", NA),
    rs2 = c("B", NA),
    rs3 = c("H", NA),
    rs4 = c("A", NA)
  )
  result <- validate_genotype_table(data)
  expect_true(any(result$warnings$code == "genotype_table_high_missingness_sample"))
})

test_that("validate_genotype_table flags high marker missingness as warning", {
  data <- make_data(
    sample_id = c("s1", "s2", "s3", "s4"),
    rs_bad = c(NA, NA, "A", "B"),
    rs_good = c("A", "B", "A", "H")
  )
  result <- validate_genotype_table(data)
  expect_true(any(result$warnings$code == "genotype_table_high_missingness_marker"))
})

test_that("validate_genotype_table returns mixed-encoding error for cross-encoding tokens", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("A", "0"),
    rs2 = c("H", "1")
  )
  result <- validate_genotype_table(data)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "genotype_table_mixed_encoding"))
})

test_that("validate_genotype_table returns unsupported-encoding error for unknown tokens", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("X", "Q"),
    rs2 = c("Z", "F")
  )
  result <- validate_genotype_table(data)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "genotype_table_unsupported_encoding"))
})

test_that("validate_genotype_table accepts explicit encoding argument", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("A", "B")
  )
  result <- validate_genotype_table(data, encoding = "AHB")
  expect_true(result$ok)
  expect_equal(result$summary$inferred_encoding, "AHB")
})

test_that("validate_genotype_table errors when explicit encoding mismatches data", {
  data <- make_data(
    sample_id = c("s1", "s2"),
    rs1 = c("A", "B")
  )
  result <- validate_genotype_table(data, encoding = "012")
  expect_false(result$ok)
  expect_true(any(result$errors$code == "genotype_table_mixed_encoding"))
})

test_that("validate_genotype_table snapshot for ok case", {
  testthat::local_reproducible_output()
  data <- make_data(
    sample_id = c("s1", "s2", "s3"),
    rs1 = c("A", "H", "B"),
    rs2 = c("H", "A", "B")
  )
  result <- validate_genotype_table(data)
  expect_snapshot(print(result))
})

test_that("validate_genotype_table snapshot for failing case", {
  testthat::local_reproducible_output()
  data <- make_data(
    sample_id = c("s1", "s1", "s2"),
    rs1 = c("A", "B", "H"),
    rs2 = c("A", "A", "A")
  )
  result <- validate_genotype_table(data, sample_col = "sample_id")
  expect_snapshot(print(result))
})
