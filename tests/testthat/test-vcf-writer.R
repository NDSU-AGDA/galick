make_tiny_vcf_args <- function() {
  list(
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100L, 200L, 50L),
    id = c("rs1", "rs2", "rs3"),
    ref = c("A", "C", "G"),
    alt = c("T", "G", "A"),
    gt_matrix = matrix(
      c(
        "0/0", "0/1", "1/1",
        "0/1", "1/1", "0/0",
        "1/1", "0/0", "0/1"
      ),
      nrow = 3, byrow = TRUE,
      dimnames = list(NULL, c("s1", "s2", "s3"))
    )
  )
}

test_that("write_vcf_v4_2 writes a valid VCF 4.2 file", {
  out <- withr::local_tempfile(fileext = ".vcf")
  args <- make_tiny_vcf_args()
  do.call(galick:::write_vcf_v4_2, c(list(output_path = out), args))

  lines <- readLines(out)
  expect_true("##fileformat=VCFv4.2" %in% lines)
  expect_true(any(grepl("^##fileDate=\\d{8}$", lines)))
  expect_true(any(grepl("^##source=galick/", lines)))
  expect_true("##contig=<ID=chr1>" %in% lines)
  expect_true("##contig=<ID=chr2>" %in% lines)
  expect_true(any(grepl("^##FORMAT=<ID=GT,", lines)))

  header_line <- lines[startsWith(lines, "#CHROM")]
  expect_length(header_line, 1L)
  expect_match(header_line, "s1\\ts2\\ts3$")

  body <- lines[!startsWith(lines, "#")]
  expect_length(body, 3L)
  expect_match(body[1], "^chr1\\t100\\trs1\\tA\\tT\\t\\.\\tPASS\\t\\.\\tGT\\t0/0\\t0/1\\t1/1$")
})

test_that("write_vcf_v4_2 always writes decimal integer positions", {
  args <- make_tiny_vcf_args()
  args$pos <- c(1e5, 2e6, 1e8)
  out <- withr::local_tempfile(fileext = ".vcf")
  do.call(write_vcf_v4_2, c(list(output_path = out), args))
  body <- readLines(out)
  body <- body[!startsWith(body, "#")]
  positions <- vapply(strsplit(body, "\t"), `[`, character(1), 2)
  expect_identical(positions, c("100000", "2000000", "100000000"))
})

test_that("write_vcf_v4_2 omits length from contig lines (per M4 spec)", {
  out <- withr::local_tempfile(fileext = ".vcf")
  args <- make_tiny_vcf_args()
  do.call(galick:::write_vcf_v4_2, c(list(output_path = out), args))

  lines <- readLines(out)
  contig_lines <- lines[grepl("^##contig=", lines)]
  expect_false(any(grepl("length=", contig_lines)))
})

test_that("write_vcf_v4_2 refuses non-.vcf extensions", {
  out <- withr::local_tempfile(fileext = ".vcf.gz")
  args <- make_tiny_vcf_args()
  expect_error(
    do.call(galick:::write_vcf_v4_2, c(list(output_path = out), args)),
    class = "galick_error_unsupported_extension"
  )
})

test_that("write_vcf_v4_2 refuses to overwrite by default", {
  out <- withr::local_tempfile(fileext = ".vcf")
  file.create(out)
  args <- make_tiny_vcf_args()
  expect_error(
    do.call(galick:::write_vcf_v4_2, c(list(output_path = out), args)),
    class = "galick_error_output_exists"
  )
})

test_that("write_vcf_v4_2 overwrites when overwrite = TRUE", {
  out <- withr::local_tempfile(fileext = ".vcf")
  file.create(out)
  args <- make_tiny_vcf_args()
  expect_no_error(
    do.call(galick:::write_vcf_v4_2, c(list(output_path = out, overwrite = TRUE), args))
  )
})

test_that("write_vcf_v4_2 errors on length-mismatched marker vectors", {
  args <- make_tiny_vcf_args()
  args$pos <- c(100L, 200L)
  out <- withr::local_tempfile(fileext = ".vcf")
  expect_error(
    do.call(galick:::write_vcf_v4_2, c(list(output_path = out), args)),
    class = "galick_error_invalid_argument"
  )
})

test_that("write_vcf_v4_2 errors on samples mismatch with gt_matrix columns", {
  args <- make_tiny_vcf_args()
  args$samples <- c("s1", "s2")
  out <- withr::local_tempfile(fileext = ".vcf")
  expect_error(
    do.call(galick:::write_vcf_v4_2, c(list(output_path = out), args)),
    class = "galick_error_invalid_argument"
  )
})
