make_dose_vcf <- function() {
  path <- withr::local_tempfile(
    fileext = ".vcf",
    .local_envir = parent.frame()
  )
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t./.",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT:DP\t0/0:30\t0/1:25\t1/1:18"
  ), path)
  path
}

test_that("vcf_to_dosage maps GT to ALT-count with samples as rows", {
  vcf <- make_dose_vcf()
  d <- vcf_to_dosage(vcf)
  expect_true(is.matrix(d))
  expect_equal(rownames(d), c("s1", "s2", "s3"))
  expect_equal(colnames(d), c("rs1", "rs2", "rs3"))
  expect_equal(unname(d["s1", "rs1"]), 0)
  expect_equal(unname(d["s2", "rs1"]), 1)
  expect_equal(unname(d["s3", "rs1"]), 2)
})

test_that("vcf_to_dosage produces NA for missing GT", {
  vcf <- make_dose_vcf()
  d <- vcf_to_dosage(vcf)
  expect_true(is.na(d["s3", "rs2"]))
})

test_that("vcf_to_dosage extracts GT from FORMAT=GT:DP cells", {
  vcf <- make_dose_vcf()
  d <- vcf_to_dosage(vcf)
  expect_equal(unname(d["s1", "rs3"]), 0)
  expect_equal(unname(d["s3", "rs3"]), 2)
})

test_that("vcf_to_dosage subsets by samples and markers", {
  vcf <- make_dose_vcf()
  d <- vcf_to_dosage(vcf, samples = c("s1", "s2"), markers = c("rs1", "rs3"))
  expect_equal(rownames(d), c("s1", "s2"))
  expect_equal(colnames(d), c("rs1", "rs3"))
})

test_that("vcf_to_dosage falls back to chrom_pos_ref_alt when ID is missing", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/1"
  ), path)
  d <- vcf_to_dosage(path)
  expect_equal(colnames(d), "chr1_100_A_T")
})

test_that("vcf_to_dosage reads .vcf.gz transparently", {
  vcf <- withr::local_tempfile(fileext = ".vcf.gz")
  con <- gzfile(vcf, "wt")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t1/1"
  ), con)
  close(con)
  d <- vcf_to_dosage(vcf)
  expect_equal(unname(d["s1", "rs1"]), 2)
})

test_that("vcf_to_dosage aborts on bad input", {
  expect_error(vcf_to_dosage(NULL), class = "beagler_error_vcf_to_dosage_invalid_input")
  expect_error(vcf_to_dosage("/no-such-file.vcf"),
    class = "beagler_error_vcf_to_dosage_invalid_input"
  )
})
