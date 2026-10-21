make_multi_chrom_vcf <- function() {
  path <- withr::local_tempfile(
    fileext = ".vcf",
    .local_envir = parent.frame()
  )
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##contig=<ID=chr2>",
    "##contig=<ID=chr3>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1",
    "chr2\t50\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/0",
    "chr2\t150\trs4\tT\tC\t.\tPASS\t.\tGT\t0/1\t0/0",
    "chr3\t300\trs5\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
  ), path)
  path
}

test_that("split_vcf_by_chrom writes one file per CHROM", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  res <- split_vcf_by_chrom(vcf, out_dir)
  expect_length(res$chrom_paths, 3L)
  expect_setequal(names(res$chrom_paths), c("chr1", "chr2", "chr3"))
  for (p in res$chrom_paths) expect_true(file.exists(p))
})

test_that("split_vcf_by_chrom records survive validate_vcf", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  res <- split_vcf_by_chrom(vcf, out_dir)
  for (p in res$chrom_paths) {
    expect_true(validate_vcf(p)$ok)
  }
})

test_that("split_vcf_by_chrom respects chrom_filter", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  res <- split_vcf_by_chrom(vcf, out_dir,
    chrom_filter = c("chr1", "chr3")
  )
  expect_setequal(names(res$chrom_paths), c("chr1", "chr3"))
})

test_that("split_vcf_by_chrom counts records correctly per chrom", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  res <- split_vcf_by_chrom(vcf, out_dir)
  expect_equal(validate_vcf(res$chrom_paths[["chr1"]])$summary$n_records, 2L)
  expect_equal(validate_vcf(res$chrom_paths[["chr2"]])$summary$n_records, 2L)
  expect_equal(validate_vcf(res$chrom_paths[["chr3"]])$summary$n_records, 1L)
})

test_that("split_vcf_by_chrom refuses to overwrite by default", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  split_vcf_by_chrom(vcf, out_dir)
  expect_error(
    split_vcf_by_chrom(vcf, out_dir),
    class = "beagler_error_split_vcf_output_exists"
  )
})

test_that("split_vcf_by_chrom overwrites when overwrite = TRUE", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  split_vcf_by_chrom(vcf, out_dir)
  expect_no_error(split_vcf_by_chrom(vcf, out_dir, overwrite = TRUE))
})

test_that("split_vcf_by_chrom aborts on missing input", {
  expect_error(
    split_vcf_by_chrom("/no-such-file.vcf", tempdir()),
    class = "beagler_error_split_vcf_invalid_input"
  )
})

test_that("split_vcf_by_chrom aborts on missing #CHROM header", {
  bad <- withr::local_tempfile(fileext = ".vcf")
  writeLines("##fileformat=VCFv4.2", bad)
  expect_error(
    split_vcf_by_chrom(bad, tempdir()),
    class = "beagler_error_split_vcf_invalid_input"
  )
})

test_that("split_vcf_by_chrom safe_chrom_filename replaces unsafe chars", {
  expect_equal(beagler:::safe_chrom_filename("chr1"), "chr1")
  expect_equal(beagler:::safe_chrom_filename("chr/1"), "chr_1")
  expect_equal(beagler:::safe_chrom_filename("1A"), "1A")
})

test_that("merge_beagle_outputs concatenates per-chrom outputs", {
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  res <- split_vcf_by_chrom(vcf, out_dir)

  fake_runs <- lapply(res$chrom_paths, function(p) {
    beagler:::new_beagler_run(
      status = 0L,
      runtime_sec = 0.1,
      command = "fake",
      output_paths = p
    )
  })

  merged <- withr::local_tempfile(fileext = ".vcf")
  merge_beagle_outputs(fake_runs, merged)

  v <- validate_vcf(merged)
  expect_true(v$ok)
  expect_equal(v$summary$n_records, 5L)
  expect_setequal(
    unique(v$summary$chromosome_summary$chrom),
    c("chr1", "chr2", "chr3")
  )
})

test_that("merge_beagle_outputs writes .vcf.gz when output extension matches", {
  skip_if(Sys.which("bgzip") == "", "bgzip required for compressed output")
  vcf <- make_multi_chrom_vcf()
  out_dir <- withr::local_tempdir()
  res <- split_vcf_by_chrom(vcf, out_dir)
  fake_runs <- lapply(res$chrom_paths, function(p) {
    beagler:::new_beagler_run(
      status = 0L, runtime_sec = 0.1,
      command = "fake", output_paths = p
    )
  })
  merged <- withr::local_tempfile(fileext = ".vcf.gz")
  merge_beagle_outputs(fake_runs, merged)
  expect_true(validate_vcf(merged)$ok)
})

test_that("merge_beagle_outputs aborts when no output paths exist", {
  fake_run <- beagler:::new_beagler_run(
    status = 0L, runtime_sec = 0.1,
    command = "fake", output_paths = "/no-such-output.vcf.gz"
  )
  expect_error(
    merge_beagle_outputs(list(fake_run), withr::local_tempfile(fileext = ".vcf")),
    class = "beagler_error_merge_beagle_invalid_input"
  )
})

test_that("merge_beagle_outputs aborts on bad output extension", {
  fake <- beagler:::new_beagler_run(status = 0L, output_paths = "/tmp/x.vcf")
  expect_error(
    merge_beagle_outputs(list(fake), "/tmp/out.txt"),
    class = "beagler_error_merge_beagle_invalid_input"
  )
})

test_that("run_beagle_by_chrom integration: one beagler_run per chrom", {
  skip_if_no_java()
  skip_if_no_beagle()
  vcf <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##contig=<ID=chr2>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
    "chr1\t300\trs3\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "chr2\t50\trs4\tT\tA\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1",
    "chr2\t100\trs5\tA\tC\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1"
  ), vcf)
  out_dir <- withr::local_tempdir()
  results <- run_beagle_by_chrom(vcf, out_dir, mode = "phase")
  expect_setequal(names(results), c("chr1", "chr2"))
  for (r in results) {
    expect_s3_class(r, "beagler_run")
    expect_equal(r$status, 0L)
  }
})
