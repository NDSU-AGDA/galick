make_tiny_vcf_args <- function() {
  list(
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100L, 200L, 50L),
    id = c("rs1", "rs2", "rs3"),
    ref = c("A", "C", "G"),
    alt = c("T", "G", "A"),
    gt_matrix = matrix(
      c("0/0", "0/1", "1/1",
        "0/1", "1/1", "0/0",
        "1/1", "0/0", "0/1"),
      nrow = 3, byrow = TRUE,
      dimnames = list(NULL, c("s1", "s2", "s3"))
    )
  )
}

test_that("read_vcf_minimal round-trips data written by write_vcf_v4_2", {
  out <- withr::local_tempfile(fileext = ".vcf")
  args <- make_tiny_vcf_args()
  do.call(beagler:::write_vcf_v4_2, c(list(output_path = out), args))

  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$chrom, args$chrom)
  expect_equal(parsed$pos, args$pos)
  expect_equal(parsed$id, args$id)
  expect_equal(parsed$ref, args$ref)
  expect_equal(parsed$alt, args$alt)
  expect_equal(parsed$samples, c("s1", "s2", "s3"))
  expect_equal(unname(parsed$gt_matrix), unname(args$gt_matrix))
  expect_equal(colnames(parsed$gt_matrix), colnames(args$gt_matrix))
})

test_that("read_vcf_minimal captures meta lines", {
  out <- withr::local_tempfile(fileext = ".vcf")
  args <- make_tiny_vcf_args()
  do.call(beagler:::write_vcf_v4_2, c(list(output_path = out), args))

  parsed <- beagler:::read_vcf_minimal(out)
  expect_true("##fileformat=VCFv4.2" %in% parsed$meta)
  expect_true(any(grepl("^##contig=", parsed$meta)))
})

test_that("read_vcf_minimal errors on missing file", {
  expect_error(
    beagler:::read_vcf_minimal(file.path(tempdir(), "no-such-vcf-xyzzy.vcf")),
    class = "beagler_error_invalid_argument"
  )
})

test_that("read_vcf_minimal errors on missing #CHROM header", {
  out <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "data\trow\tno\theader"), out)
  expect_error(
    beagler:::read_vcf_minimal(out),
    class = "beagler_error_invalid_argument"
  )
})

test_that("read_vcf_minimal handles header-only file (no records)", {
  out <- withr::local_tempfile(fileext = ".vcf")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##contig=<ID=chr1>",
      paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "s1"),
            collapse = "\t")
    ),
    out
  )
  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$samples, "s1")
  expect_equal(nrow(parsed$gt_matrix), 0L)
})
