test_that("new_galick_result with no issues is ok", {
  r <- galick:::new_galick_result()
  expect_s3_class(r, "galick_result")
  expect_true(r$ok)
  expect_equal(nrow(r$errors), 0L)
  expect_equal(nrow(r$warnings), 0L)
})

test_that("new_galick_result auto-derives ok from error count", {
  errs <- tibble::tibble(
    code = "test_code",
    message = "test message",
    where = list(list(line = 1L))
  )
  r <- galick:::new_galick_result(errors = errs)
  expect_false(r$ok)
})

test_that("new_galick_result accepts explicit ok override", {
  errs <- tibble::tibble(code = "x", message = "x", where = list(list()))
  r <- galick:::new_galick_result(errors = errs, ok = TRUE)
  expect_true(r$ok)
})

test_that("new_galick_result captures call into metadata", {
  call <- quote(test_function(x = 1))
  r <- galick:::new_galick_result(call = call)
  expect_identical(r$metadata$call, call)
})

test_that("new_galick_result metadata has package_version and time", {
  r <- galick:::new_galick_result()
  expect_s3_class(r$metadata$package_version, "package_version")
  expect_match(r$metadata$time, "^\\d{4}-\\d{2}-\\d{2}T")
})

test_that("normalize_issues fills missing where column", {
  df <- data.frame(code = "x", message = "y", stringsAsFactors = FALSE)
  out <- galick:::normalize_issues(df)
  expect_true("where" %in% names(out))
  expect_true(is.list(out$where))
})

test_that("normalize_issues errors on missing required columns", {
  expect_error(
    galick:::normalize_issues(data.frame(only = "this")),
    class = "galick_error_invalid_issues"
  )
})

test_that("print.galick_result snapshot for ok result", {
  testthat::local_reproducible_output()
  r <- galick:::new_galick_result(
    summary = list(n_samples = 100L, n_markers = 50L)
  )
  expect_snapshot(print(r))
})

test_that("print.galick_result snapshot for failing result", {
  testthat::local_reproducible_output()
  errs <- tibble::tibble(
    code = c("vcf_unsorted", "vcf_dup_id"),
    message = c("Records out of order", "Marker ID rs1 appears twice"),
    where = list(
      list(chrom = "chr1", pos = 100L),
      list(marker = "rs1")
    )
  )
  warns <- tibble::tibble(
    code = "vcf_low_maf",
    message = "MAF below 0.01 on 5 markers",
    where = list(list())
  )
  r <- galick:::new_galick_result(
    errors = errs,
    warnings = warns,
    summary = list(n_records = 1000L)
  )
  expect_snapshot(print(r))
})

test_that("as.data.frame.galick_result combines errors and warnings", {
  errs <- tibble::tibble(code = "e1", message = "msg1", where = list(list()))
  warns <- tibble::tibble(code = "w1", message = "msg2", where = list(list()))
  r <- galick:::new_galick_result(errors = errs, warnings = warns)
  df <- as.data.frame(r)
  expect_equal(nrow(df), 2L)
  expect_equal(df$severity, c("error", "warning"))
})

test_that("as.data.frame.galick_result flattens uniform where keys to where_*", {
  errs <- tibble::tibble(
    code = c("e1", "e2"),
    message = c("m1", "m2"),
    where = list(
      list(chrom = "chr1", pos = 100L),
      list(chrom = "chr2", pos = 200L)
    )
  )
  r <- galick:::new_galick_result(errors = errs)
  df <- as.data.frame(r)
  expect_true("where_chrom" %in% names(df))
  expect_true("where_pos" %in% names(df))
  expect_equal(df$where_chrom, c("chr1", "chr2"))
})

test_that("as.data.frame.galick_result falls back to JSON when where keys differ", {
  errs <- tibble::tibble(
    code = c("e1", "e2"),
    message = c("m1", "m2"),
    where = list(
      list(chrom = "chr1"),
      list(marker = "rs1")
    )
  )
  r <- galick:::new_galick_result(errors = errs)
  df <- as.data.frame(r)
  expect_true("where" %in% names(df))
  expect_match(df$where[1], "chr1")
  expect_match(df$where[2], "rs1")
})

test_that("format.galick_result returns a one-line summary", {
  r <- galick:::new_galick_result()
  expect_match(format(r), "^<galick_result")
})
