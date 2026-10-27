test_that("parse_vcf_raw transparently reads .vcf.gz", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  )
  path <- withr::local_tempfile(fileext = ".vcf.gz")
  con <- gzfile(path, "wt")
  writeLines(lines, con)
  close(con)

  parsed <- galick:::parse_vcf_raw(path)
  expect_equal(parsed$samples, "s1")
  expect_length(parsed$body_lines, 1L)
})

test_that("validate_vcf summary$has_vcf_gz_extension reflects input extension", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  )

  plain <- withr::local_tempfile(fileext = ".vcf")
  writeLines(lines, plain)
  expect_false(validate_vcf(plain)$summary$has_vcf_gz_extension)

  gz <- withr::local_tempfile(fileext = ".vcf.gz")
  con <- gzfile(gz, "wt")
  writeLines(lines, con)
  close(con)
  expect_true(validate_vcf(gz)$summary$has_vcf_gz_extension)
})
