test_that("build_beagle_argv produces the expected argv vector for impute mode", {
  args <- beagler:::build_beagle_argv(
    beagle_jar = "/p/beagle.jar",
    gt = "/p/in.vcf",
    out = "/p/out_prefix",
    nthreads = 4L,
    memory_gb = 8L
  )
  expect_equal(args[1], "-Xmx8g")
  expect_equal(args[2], "-jar")
  expect_equal(args[3], "/p/beagle.jar")
  expect_true("gt=/p/in.vcf" %in% args)
  expect_true("out=/p/out_prefix" %in% args)
  expect_true("nthreads=4" %in% args)
  expect_false(any(grepl("^impute=", args)))
})

test_that("build_beagle_argv adds impute=false in phase mode", {
  args <- beagler:::build_beagle_argv(
    beagle_jar = "/p/b.jar",
    gt = "/p/in.vcf",
    out = "/p/out",
    mode = "phase"
  )
  expect_true("impute=false" %in% args)
})

test_that("build_beagle_argv passes ref, map, seed only when given", {
  args <- beagler:::build_beagle_argv(
    beagle_jar = "/p/b.jar",
    gt = "/p/in.vcf",
    out = "/p/out",
    ref = "/p/ref.vcf.gz",
    map = "/p/genetic.map",
    seed = 42L
  )
  expect_true("ref=/p/ref.vcf.gz" %in% args)
  expect_true("map=/p/genetic.map" %in% args)
  expect_true("seed=42" %in% args)
})

test_that("build_beagle_argv appends extra_args verbatim after defaults", {
  args <- beagler:::build_beagle_argv(
    beagle_jar = "/p/b.jar",
    gt = "/p/in.vcf",
    out = "/p/out",
    extra_args = c("burnin=8", "iterations=20")
  )
  expect_equal(tail(args, 2), c("burnin=8", "iterations=20"))
})

test_that("run_beagle rejects out endings of .vcf or .vcf.gz", {
  gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  ), gt)
  expect_error(
    run_beagle(gt = gt, out = "/tmp/foo.vcf"),
    class = "beagler_error_run_beagle_invalid_input"
  )
  expect_error(
    run_beagle(gt = gt, out = "/tmp/foo.vcf.gz"),
    class = "beagler_error_run_beagle_invalid_input"
  )
})

test_that("run_beagle accepts gt of .vcf.gz", {
  skip_if_no_java()
  skip_if_no_beagle()
  gt <- withr::local_tempfile(fileext = ".vcf.gz")
  con <- gzfile(gt, "wt")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
  ), con)
  close(con)
  out_dir <- withr::local_tempdir()
  out_prefix <- file.path(out_dir, "result")
  result <- run_beagle(gt = gt, out = out_prefix, mode = "phase")
  expect_s3_class(result, "beagler_run")
})

test_that("run_beagle rejects an empty reference VCF", {
  skip_if_no_java()
  skip_if_no_beagle()
  gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  ), gt)
  bad_ref <- withr::local_tempfile(fileext = ".vcf")
  file.create(bad_ref)
  out <- withr::local_tempfile()
  expect_error(
    run_beagle(gt = gt, out = out, ref = bad_ref),
    class = "beagler_error_run_beagle_invalid_input"
  )
})

test_that("run_beagle preflight raises when gt fails validate_vcf(target='beagle')", {
  skip_if_no_java()
  skip_if_no_beagle()
  bad_gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\tnotnumeric\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  ), bad_gt)
  out <- withr::local_tempfile()
  expect_error(
    run_beagle(gt = bad_gt, out = out),
    class = "beagler_error_run_beagle_invalid_input"
  )
})

test_that("run_beagle refuses to overwrite <out>.vcf.gz unless overwrite = TRUE", {
  skip_if_no_java()
  skip_if_no_beagle()
  gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  ), gt)
  out_dir <- withr::local_tempdir()
  out_prefix <- file.path(out_dir, "result")
  file.create(paste0(out_prefix, ".vcf.gz"))
  expect_error(
    run_beagle(gt = gt, out = out_prefix),
    class = "beagler_error_run_beagle_output_exists"
  )
})

test_that("run_beagle integration: phase a tiny VCF and return ok beagler_run", {
  skip_if_no_java()
  skip_if_no_beagle()
  gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t400\trs4\tT\tA\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1",
    "chr1\t500\trs5\tA\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1"
  ), gt)
  out_dir <- withr::local_tempdir()
  out_prefix <- file.path(out_dir, "phased")
  result <- run_beagle(gt = gt, out = out_prefix, mode = "phase")

  expect_s3_class(result, "beagler_run")
  expect_equal(result$status, 0L)
  expect_true(file.exists(paste0(out_prefix, ".vcf.gz")))
  expect_true(file.exists(result$run_record))
  expect_true(file.exists(result$stdout_log))
  expect_true(file.exists(result$stderr_log))

  rec <- jsonlite::fromJSON(result$run_record, simplifyVector = FALSE)
  expect_equal(rec$schema_version, "1.0.0")
  expect_equal(rec$status, 0L)
  expect_match(rec$inputs[[1]]$sha256, "^[0-9a-f]{64}$")
})

test_that("run_beagle integration: bad extra_args produce a non-zero beagler_run", {
  skip_if_no_java()
  skip_if_no_beagle()
  gt <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  ), gt)
  out_dir <- withr::local_tempdir()
  out_prefix <- file.path(out_dir, "bad_run")

  result <- run_beagle(
    gt = gt,
    out = out_prefix,
    extra_args = c("nonsense_flag=true")
  )

  expect_s3_class(result, "beagler_run")
  expect_true(result$status != 0L)
})
