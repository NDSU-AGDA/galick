write_vcf_lines <- function(lines) {
  path <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  writeLines(lines, path)
  path
}

minimal_vcf_lines <- function() {
  c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  )
}

test_that("parse_vcf_raw splits meta, header, and body for a valid file", {
  path <- write_vcf_lines(minimal_vcf_lines())
  parsed <- beagler:::parse_vcf_raw(path)

  expect_equal(parsed$meta, c("##fileformat=VCFv4.2", "##contig=<ID=chr1>"))
  expect_equal(parsed$fileformat_line, "##fileformat=VCFv4.2")
  expect_match(parsed$chrom_header_line, "^#CHROM")
  expect_equal(parsed$samples, c("s1", "s2"))
  expect_length(parsed$body_lines, 1L)
  expect_equal(parsed$body_col_counts, 11L)
})

test_that("parse_vcf_raw tolerates a file with no #CHROM header", {
  path <- write_vcf_lines(c("##fileformat=VCFv4.2", "chr1\t100\trs1\tA\tT"))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_true(is.na(parsed$chrom_header_line))
  expect_equal(parsed$header_cols, character(0))
  expect_equal(parsed$samples, character(0))
  expect_length(parsed$body_lines, 1L)
})

test_that("parse_vcf_raw tolerates a file with no ##fileformat line", {
  path <- write_vcf_lines(c(
    "##contig=<ID=chr1>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  ))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_true(is.na(parsed$fileformat_line))
  expect_false(is.na(parsed$chrom_header_line))
})

test_that("parse_vcf_raw tolerates an empty body (header-only)", {
  path <- write_vcf_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1"
  ))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_length(parsed$body_lines, 0L)
  expect_equal(parsed$samples, "s1")
})

test_that("parse_vcf_raw tolerates header with fewer than 9 fixed columns (no samples)", {
  path <- write_vcf_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT"
  ))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_equal(parsed$samples, character(0))
  expect_equal(length(parsed$header_cols), 5L)
})

test_that("parse_vcf_raw preserves per-row column counts when rows are short or wide", {
  path <- write_vcf_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\textra"
  ))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_equal(parsed$body_col_counts, c(11L, 8L, 12L))
})

test_that("parse_vcf_raw preserves raw POS strings (no integer coercion)", {
  path <- write_vcf_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\tnotnumeric\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0",
    "chr1\t-5\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1"
  ))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_equal(parsed$body_split[[1]][2], "notnumeric")
  expect_equal(parsed$body_split[[2]][2], "-5")
})

test_that("parse_vcf_raw aborts only on missing path or file", {
  expect_error(
    beagler:::parse_vcf_raw(NA_character_),
    class = "beagler_error_invalid_argument"
  )
  expect_error(
    beagler:::parse_vcf_raw(file.path(tempdir(), "absent_xyzzy.vcf")),
    class = "beagler_error_invalid_argument"
  )
})

test_that("parse_vcf_raw skips blank lines in body", {
  path <- write_vcf_lines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0",
    ""
  ))
  parsed <- beagler:::parse_vcf_raw(path)
  expect_length(parsed$body_lines, 1L)
})
