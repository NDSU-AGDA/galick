make_bridge_vcf <- function() {
  path <- withr::local_tempfile(fileext = ".vcf",
                                 .local_envir = parent.frame())
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    "1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t0/1"
  ), path)
  path
}

test_that("vcf_to_plink aborts on missing input file", {
  expect_error(
    vcf_to_plink("/nonexistent.vcf", "/tmp/out"),
    class = "beagler_error_vcf_to_plink_invalid_input"
  )
})

test_that("vcf_to_plink aborts on bad extension", {
  bad <- withr::local_tempfile(fileext = ".txt")
  file.create(bad)
  expect_error(
    vcf_to_plink(bad, "/tmp/out"),
    class = "beagler_error_vcf_to_plink_invalid_input"
  )
})

test_that("vcf_to_plink aborts on bad plink_prefix", {
  vcf <- make_bridge_vcf()
  expect_error(
    vcf_to_plink(vcf, NA_character_),
    class = "beagler_error_vcf_to_plink_invalid_input"
  )
})

test_that("vcf_to_plink refuses to overwrite existing PLINK files", {
  vcf <- make_bridge_vcf()
  out_dir <- withr::local_tempdir()
  prefix <- file.path(out_dir, "out")
  file.create(paste0(prefix, ".bed"))
  expect_error(
    vcf_to_plink(vcf, prefix),
    class = "beagler_error_vcf_to_plink_output_exists"
  )
})

test_that("plink_to_vcf aborts on missing PLINK input files", {
  vcf <- withr::local_tempfile(fileext = ".vcf")
  expect_error(
    plink_to_vcf("/no/such/prefix", vcf),
    class = "beagler_error_plink_to_vcf_invalid_input"
  )
})

test_that("plink_to_vcf aborts on bad output extension", {
  out_dir <- withr::local_tempdir()
  prefix <- file.path(out_dir, "in")
  file.create(paste0(prefix, c(".bed", ".bim", ".fam")))
  expect_error(
    plink_to_vcf(prefix, "out.txt"),
    class = "beagler_error_plink_to_vcf_invalid_input"
  )
})

test_that("plink_to_vcf refuses to overwrite", {
  out_dir <- withr::local_tempdir()
  prefix <- file.path(out_dir, "in")
  file.create(paste0(prefix, c(".bed", ".bim", ".fam")))
  vcf <- withr::local_tempfile(fileext = ".vcf")
  file.create(vcf)
  expect_error(
    plink_to_vcf(prefix, vcf),
    class = "beagler_error_plink_to_vcf_output_exists"
  )
})

test_that("vcf_to_plink integration: writes BED/BIM/FAM", {
  skip_if_no_plink()
  vcf <- make_bridge_vcf()
  out_dir <- withr::local_tempdir()
  prefix <- file.path(out_dir, "out")
  result <- vcf_to_plink(vcf, prefix)
  expect_s3_class(result, "beagler_run")
  expect_equal(result$status, 0L)
  expect_true(all(file.exists(paste0(prefix, c(".bed", ".bim", ".fam")))))
})

test_that("plink_to_vcf integration: round-trips back to a parseable VCF", {
  skip_if_no_plink()
  vcf <- make_bridge_vcf()
  out_dir <- withr::local_tempdir()
  prefix <- file.path(out_dir, "round")
  vcf_to_plink(vcf, prefix)
  out_vcf <- file.path(out_dir, "round.vcf")
  result <- plink_to_vcf(prefix, out_vcf)
  expect_s3_class(result, "beagler_run")
  expect_equal(result$status, 0L)
  expect_true(file.exists(out_vcf))
  validation <- validate_vcf(out_vcf, target = "general")
  expect_true(validation$ok || nrow(validation$errors) == 0L)
})
