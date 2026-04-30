make_vcf <- function(lines, ext = ".vcf") {
  path <- withr::local_tempfile(fileext = ext, .local_envir = parent.frame())
  if (ext == ".vcf.gz") {
    con <- gzfile(path, "wt")
    writeLines(lines, con)
    close(con)
  } else {
    writeLines(lines, path)
  }
  path
}

simple_vcf <- function(records) {
  c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    records
  )
}

test_that("validate_reference_match returns ok for identical compatible VCFs", {
  records <- c(
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1",
    "chr2\t50\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1"
  )
  t <- make_vcf(simple_vcf(records))
  r <- make_vcf(simple_vcf(records))
  result <- validate_reference_match(t, r)
  expect_s3_class(result, "beagler_result")
  expect_true(result$ok)
  expect_equal(result$summary$n_markers_shared, 3L)
  expect_equal(result$summary$n_chroms_shared, 2L)
})

test_that("validate_reference_match warns on chrom-style mismatch under auto", {
  records_t <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  records_r <- c("1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records_t))
  r <- make_vcf(simple_vcf(records_r))
  result <- validate_reference_match(t, r, chrom_style = "auto")
  expect_true(any(result$warnings$code == "reference_match_chrom_style_mismatch"))
})

test_that("validate_reference_match under chrom_style='numeric' resolves chr/numeric difference", {
  records_t <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  records_r <- c("1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records_t))
  r <- make_vcf(simple_vcf(records_r))
  result <- validate_reference_match(t, r, chrom_style = "numeric")
  expect_true(result$ok)
  expect_equal(result$summary$n_markers_shared, 1L)
})

test_that("validate_reference_match flags REF mismatches as errors", {
  records_t <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  records_r <- c("chr1\t100\trs1\tC\tT\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records_t))
  r <- make_vcf(simple_vcf(records_r))
  result <- validate_reference_match(t, r)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "reference_match_ref_mismatch"))
})

test_that("validate_reference_match flags ALT mismatches as errors", {
  records_t <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  records_r <- c("chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records_t))
  r <- make_vcf(simple_vcf(records_r))
  result <- validate_reference_match(t, r)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "reference_match_alt_mismatch"))
})

test_that("validate_reference_match errors when no chromosomes are shared", {
  records_t <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  records_r <- c("chr2\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records_t))
  r <- make_vcf(simple_vcf(records_r))
  result <- validate_reference_match(t, r)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "reference_match_no_shared_chromosomes"))
})

test_that("validate_reference_match warns when chromosomes overlap but markers do not", {
  records_t <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  records_r <- c("chr1\t999\trs9\tG\tC\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records_t))
  r <- make_vcf(simple_vcf(records_r))
  result <- validate_reference_match(t, r)
  expect_true(any(result$warnings$code == "reference_match_no_shared_markers"))
})

test_that("validate_reference_match accepts .vcf.gz inputs", {
  records <- c("chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0")
  t <- make_vcf(simple_vcf(records))
  r <- make_vcf(simple_vcf(records), ext = ".vcf.gz")
  result <- validate_reference_match(t, r)
  expect_true(result$ok)
})

test_that("validate_reference_match aborts on missing file", {
  expect_error(
    validate_reference_match("/tmp/no-such-target.vcf", "/tmp/no-such-ref.vcf"),
    class = "beagler_error_reference_match_invalid_input"
  )
})

test_that("validate_reference_match aborts on bad extension", {
  bad <- withr::local_tempfile(fileext = ".txt")
  file.create(bad)
  expect_error(
    validate_reference_match(bad, bad),
    class = "beagler_error_reference_match_invalid_input"
  )
})
