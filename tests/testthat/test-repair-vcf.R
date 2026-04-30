write_vcf_lines <- function(lines) {
  path <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  writeLines(lines, path)
  path
}

clean_vcf_lines <- function() {
  c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##contig=<ID=chr2>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1",
    "chr2\t50\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/1"
  )
}

test_that("repair_vcf is idempotent on a clean valid VCF", {
  input <- write_vcf_lines(clean_vcf_lines())
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)

  expect_s3_class(result, "beagler_result")
  expect_true(result$ok)
  expect_equal(result$summary$n_records_in, 3L)
  expect_equal(result$summary$n_records_out, 3L)
  expect_equal(result$summary$output_path, output)

  parsed <- beagler:::read_vcf_minimal(output)
  expect_equal(parsed$id, c("rs1", "rs2", "rs3"))
})

test_that("repair_vcf step drop_invalid_records removes rows with invalid POS", {
  lines <- clean_vcf_lines()
  lines <- c(lines, "chr2\t-5\trs_bad\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/1")
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = "drop_invalid_records")

  expect_true(result$ok)
  expect_equal(result$summary$n_records_in, 4L)
  expect_equal(result$summary$n_records_out, 3L)
  parsed <- beagler:::read_vcf_minimal(output)
  expect_false("rs_bad" %in% parsed$id)
})

test_that("repair_vcf step drop_invalid_records removes REF == ALT rows", {
  lines <- clean_vcf_lines()
  lines <- c(lines, "chr2\t60\trs_bad\tG\tG\t.\tPASS\t.\tGT\t0/0\t0/1")
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = "drop_invalid_records")
  expect_true(result$ok)
  expect_equal(result$summary$n_records_out, 3L)
})

test_that("repair_vcf step drop_invalid_records removes non-diploid GT", {
  lines <- clean_vcf_lines()
  lines <- c(lines, "chr2\t60\trs_bad\tG\tT\t.\tPASS\t.\tGT\t0\t1")
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = "drop_invalid_records")
  expect_true(result$ok)
  expect_equal(result$summary$n_records_out, 3L)
})

test_that("repair_vcf step filter_biallelic removes multiallelic rows", {
  lines <- clean_vcf_lines()
  lines <- c(lines, "chr2\t60\trs_multi\tA\tT,C\t.\tPASS\t.\tGT\t0/0\t0/1")
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = c("drop_invalid_records", "filter_biallelic"))
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(output)
  expect_false("rs_multi" %in% parsed$id)
})

test_that("repair_vcf step fix_marker_ids fills '.' with chrom_pos_ref_alt", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = c("drop_invalid_records", "fix_marker_ids"))
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(output)
  expect_true("chr1_100_A_T" %in% parsed$id)
})

test_that("repair_vcf step remove_duplicates drops by ID and (chrom, pos)", {
  lines <- clean_vcf_lines()
  lines <- c(
    lines,
    "chr2\t60\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1",
    "chr1\t100\trs9\tA\tT\t.\tPASS\t.\tGT\t0/1\t1/1"
  )
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = c("drop_invalid_records", "remove_duplicates"))
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(output)
  expect_false(any(duplicated(parsed$id)))
})

test_that("repair_vcf step sort orders by (chrom, pos)", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr2\t50\trs3\tG\tA\t.\tPASS\t.\tGT\t1/1\t0/1"
  lines[7] <- "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  lines[8] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = "sort")
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(output)
  expect_equal(parsed$chrom, c("chr1", "chr1", "chr2"))
  expect_equal(parsed$pos, c(100L, 200L, 50L))
})

test_that("repair_vcf rewrites GT:DP cells to GT-only", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT:DP\t0/0:30\t1/1:25"
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)
  expect_true(result$ok)
  out_lines <- readLines(output)
  body <- out_lines[!startsWith(out_lines, "#")]
  expect_match(body[1], "GT\t0/0\t1/1")
  expect_false(any(grepl("DP", out_lines)))
})

test_that("repair_vcf chrom_style = 'numeric' rewrites CHROM in output", {
  lines <- clean_vcf_lines()
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, chrom_style = "numeric")
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(output)
  expect_equal(parsed$chrom, c("1", "1", "2"))
})

test_that("repair_vcf chrom_style = 'numeric' catches chr1+1 duplicates", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0",
    "1\t100\trs2\tA\tT\t.\tPASS\t.\tGT\t0/1"
  )
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, chrom_style = "numeric")
  expect_true(result$ok)
  expect_equal(result$summary$n_records_out, 1L)
})

test_that("repair_vcf returns repair_vcf_no_records_left when everything is dropped", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1>",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t-5\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  )
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "repair_vcf_no_records_left"))
  expect_false(file.exists(output))
})

test_that("repair_vcf aborts on missing #CHROM (file-level blocker)", {
  lines <- c("##fileformat=VCFv4.2", "chr1\t100\trs1\tA\tT")
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "repair_vcf_missing_chrom_header"))
  expect_false(file.exists(output))
})

test_that("repair_vcf aborts on no samples (file-level blocker)", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT"
  )
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "repair_vcf_no_samples"))
  expect_false(file.exists(output))
})

test_that("repair_vcf aborts on duplicate sample names (file-level blocker)", {
  lines <- clean_vcf_lines()
  lines[5] <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts1"
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "repair_vcf_duplicated_samples"))
  expect_false(file.exists(output))
})

test_that("repair_vcf raises on output_vcf existing without overwrite", {
  input <- write_vcf_lines(clean_vcf_lines())
  output <- withr::local_tempfile(fileext = ".vcf")
  file.create(output)
  expect_error(
    repair_vcf(input, output),
    class = "beagler_error_repair_vcf_output_exists"
  )
})

test_that("repair_vcf overwrites when overwrite = TRUE", {
  input <- write_vcf_lines(clean_vcf_lines())
  output <- withr::local_tempfile(fileext = ".vcf")
  file.create(output)
  result <- repair_vcf(input, output, overwrite = TRUE)
  expect_true(result$ok)
  expect_gt(file.size(output), 0)
})

test_that("repair_vcf raises on bad output extension", {
  input <- write_vcf_lines(clean_vcf_lines())
  bad_out <- withr::local_tempfile(fileext = ".vcf.gz")
  expect_error(
    repair_vcf(input, bad_out),
    class = "beagler_error_repair_vcf_unsupported_extension"
  )
})

test_that("repair_vcf raises on unknown step name", {
  input <- write_vcf_lines(clean_vcf_lines())
  output <- withr::local_tempfile(fileext = ".vcf")
  expect_error(
    repair_vcf(input, output, steps = c("not_a_step")),
    class = "beagler_error_repair_vcf_invalid_input"
  )
})

test_that("repair_vcf surfaces repair_vcf_remaining_* when user excludes a needed step", {
  # Multiallelic record present; user excludes filter_biallelic, so the
  # post-validation should catch it and report remaining errors.
  lines <- clean_vcf_lines()
  lines <- c(lines, "chr2\t60\trs_multi\tA\tT,C\t.\tPASS\t.\tGT\t0/0\t0/1")
  input <- write_vcf_lines(lines)
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output, steps = c("drop_invalid_records"))
  expect_false(result$ok)
  expect_true(any(grepl("^repair_vcf_remaining_", result$errors$code)))
  expect_false(file.exists(output))
})

test_that("repair_vcf end-to-end on inst/extdata/tiny_dirty.vcf passes validate_vcf(target = 'beagle')", {
  dirty <- system.file("extdata", "tiny_dirty.vcf", package = "beagler")
  skip_if(dirty == "", "tiny_dirty.vcf fixture not installed")
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(dirty, output)

  expect_true(result$ok, info = paste("repair_vcf errors:",
                                      paste(result$errors$code, collapse = ", ")))
  expect_true(file.exists(output))

  v <- validate_vcf(output, target = "beagle")
  expect_equal(nrow(v$errors), 0L,
               info = paste(v$errors$code, collapse = ", "))
})

test_that("repair_vcf step_summary tibble has one row per step", {
  input <- write_vcf_lines(clean_vcf_lines())
  output <- withr::local_tempfile(fileext = ".vcf")
  result <- repair_vcf(input, output)
  expect_equal(nrow(result$summary$step_summary), 5L)
  expect_setequal(
    result$summary$step_summary$step,
    c("drop_invalid_records", "filter_biallelic",
      "fix_marker_ids", "remove_duplicates", "sort")
  )
})
