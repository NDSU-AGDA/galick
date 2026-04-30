make_hc_vcf <- function() {
  path <- withr::local_tempfile(fileext = ".vcf",
                                 .local_envir = parent.frame())
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t1/1\t./.\t0/0"
  ), path)
  path
}

test_that("vcf_to_hardcalls 012 encoding (default)", {
  vcf <- make_hc_vcf()
  hc <- vcf_to_hardcalls(vcf)
  expect_equal(unname(hc["s1", "rs1"]), "0")
  expect_equal(unname(hc["s2", "rs1"]), "1")
  expect_equal(unname(hc["s3", "rs1"]), "2")
})

test_that("vcf_to_hardcalls AHB encoding", {
  vcf <- make_hc_vcf()
  hc <- vcf_to_hardcalls(vcf, encoding = "AHB")
  expect_equal(unname(hc["s1", "rs1"]), "A")
  expect_equal(unname(hc["s2", "rs1"]), "H")
  expect_equal(unname(hc["s3", "rs1"]), "B")
})

test_that("vcf_to_hardcalls AABB encoding", {
  vcf <- make_hc_vcf()
  hc <- vcf_to_hardcalls(vcf, encoding = "AABB")
  expect_equal(unname(hc["s1", "rs1"]), "AA")
  expect_equal(unname(hc["s2", "rs1"]), "AB")
  expect_equal(unname(hc["s3", "rs1"]), "BB")
})

test_that("vcf_to_hardcalls produces NA for missing GT", {
  vcf <- make_hc_vcf()
  hc <- vcf_to_hardcalls(vcf, encoding = "AHB")
  expect_true(is.na(hc["s2", "rs2"]))
})

test_that("vcf_to_hardcalls subsets by samples and markers", {
  vcf <- make_hc_vcf()
  hc <- vcf_to_hardcalls(vcf, samples = c("s1"), markers = c("rs2"),
                         encoding = "AHB")
  expect_equal(dim(hc), c(1L, 1L))
  expect_equal(unname(hc[1, 1]), "B")
})
