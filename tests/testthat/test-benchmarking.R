make_bench_vcf <- function() {
  path <- withr::local_tempfile(
    fileext = ".vcf",
    .local_envir = parent.frame()
  )
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1\t1/1",
    "chr1\t400\trs4\tT\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1"
  ), path)
  path
}

test_that("mask_genotypes masks roughly the requested fraction", {
  vcf <- make_bench_vcf()
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- mask_genotypes(vcf, out, fraction = 0.5, seed = 42L)
  expect_true(result$ok)
  expect_equal(result$summary$n_eligible, 16L)
  expect_equal(result$summary$n_masked, 8L)
  expect_equal(result$summary$output_path, out)
})

test_that("mask_genotypes is reproducible with seed", {
  vcf <- make_bench_vcf()
  out1 <- withr::local_tempfile(fileext = ".vcf")
  out2 <- withr::local_tempfile(fileext = ".vcf")
  mask_genotypes(vcf, out1, fraction = 0.25, seed = 7L)
  mask_genotypes(vcf, out2, fraction = 0.25, seed = 7L)
  d1 <- vcf_to_dosage(out1)
  d2 <- vcf_to_dosage(out2)
  expect_equal(d1, d2)
})

test_that("mask_genotypes produces .vcf with masked cells", {
  vcf <- make_bench_vcf()
  out <- withr::local_tempfile(fileext = ".vcf")
  mask_genotypes(vcf, out, fraction = 0.5, seed = 1L)
  d_orig <- vcf_to_dosage(vcf)
  d_mask <- vcf_to_dosage(out)
  n_newly_missing <- sum(!is.na(d_orig) & is.na(d_mask))
  expect_gt(n_newly_missing, 0L)
})

test_that("mask_genotypes refuses to overwrite by default", {
  vcf <- make_bench_vcf()
  out <- withr::local_tempfile(fileext = ".vcf")
  file.create(out)
  expect_error(
    mask_genotypes(vcf, out, fraction = 0.1),
    class = "beagler_error_mask_genotypes_output_exists"
  )
})

test_that("mask_genotypes rejects out-of-range fraction", {
  vcf <- make_bench_vcf()
  out <- withr::local_tempfile(fileext = ".vcf")
  expect_error(
    mask_genotypes(vcf, out, fraction = 0),
    class = "beagler_error_mask_genotypes_invalid_input"
  )
  expect_error(
    mask_genotypes(vcf, out, fraction = 1.5),
    class = "beagler_error_mask_genotypes_invalid_input"
  )
})

test_that("compare_imputed_vs_true reports perfect concordance for identical VCFs", {
  vcf <- make_bench_vcf()
  result <- compare_imputed_vs_true(vcf, vcf)
  expect_true(result$ok)
  expect_equal(result$summary$concordance_rate, 1)
  expect_equal(result$summary$rmse, 0)
})

test_that("compare_imputed_vs_true with masked_vcf restricts to masked cells", {
  vcf <- make_bench_vcf()
  masked <- withr::local_tempfile(fileext = ".vcf")
  mask_genotypes(vcf, masked, fraction = 0.5, seed = 9L)

  # Compare 'imputed' = original (perfect) to truth, restricted to the masked cells.
  result <- compare_imputed_vs_true(
    imputed_vcf = vcf,
    true_vcf = vcf,
    masked_vcf = masked
  )
  expect_true(result$ok)
  expect_true(result$summary$used_mask_filter)
  expect_equal(result$summary$concordance_rate, 1)
  expect_gt(result$summary$n_compared, 0L)
})

test_that("compare_imputed_vs_true penalises mismatches", {
  vcf <- make_bench_vcf()
  # Build an "imputed" file that flips half the genotypes
  parsed <- beagler:::parse_vcf_raw(vcf)
  rec <- beagler:::extract_vcf_records(parsed)
  flipped <- rec$gt_matrix
  flipped[1, ] <- c("1/1", "1/1", "0/0", "1/1") # flip rs1 cells
  bad <- withr::local_tempfile(fileext = ".vcf")
  beagler:::write_vcf_v4_2(
    output_path = bad,
    chrom = rec$metadata$chrom,
    pos = rec$metadata$pos,
    id = rec$metadata$id,
    ref = rec$metadata$ref,
    alt = rec$metadata$alt,
    gt_matrix = flipped,
    samples = rec$samples,
    overwrite = TRUE
  )
  result <- compare_imputed_vs_true(bad, vcf)
  expect_lt(result$summary$concordance_rate, 1)
  expect_gt(result$summary$rmse, 0)
})

test_that("compare_imputed_vs_true aborts on missing input", {
  expect_error(
    compare_imputed_vs_true("/no.vcf", "/also-no.vcf"),
    class = "beagler_error_compare_imputed_invalid_input"
  )
})

test_that("benchmark_imputation: integration via Beagle", {
  skip_if_no_java()
  skip_if_no_beagle()
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0\t0/1",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1\t1/1",
    "chr1\t400\trs4\tT\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1",
    "chr1\t500\trs5\tA\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1\t1/1"
  ), vcf)
  out_dir <- withr::local_tempdir()
  result <- benchmark_imputation(vcf, out_dir,
    fraction = 0.2, seed = 42L,
    mode = "phase"
  )
  expect_named(result, c("mask", "run", "comparison"))
  expect_s3_class(result$run, "beagler_run")
  expect_s3_class(result$mask, "beagler_result")
  expect_s3_class(result$comparison, "beagler_result")
})
