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

test_that("validate_vcf happy path returns ok", {
  path <- write_vcf_lines(clean_vcf_lines())
  result <- validate_vcf(path)
  expect_s3_class(result, "galick_result")
  expect_true(result$ok)
  expect_equal(result$summary$n_samples, 2L)
  expect_equal(result$summary$n_records, 3L)
  expect_equal(result$summary$inferred_chrom_style, "chr_prefixed")
})

test_that("validate_vcf accepts .vcf.gz transparently", {
  path <- withr::local_tempfile(fileext = ".vcf.gz")
  con <- gzfile(path, "wt")
  writeLines(clean_vcf_lines(), con)
  close(con)
  result <- validate_vcf(path)
  expect_true(result$ok)
  expect_true(result$summary$has_vcf_gz_extension)
})

test_that("validate_vcf rejects unsupported extensions", {
  path <- withr::local_tempfile(fileext = ".txt")
  file.create(path)
  expect_error(
    validate_vcf(path),
    class = "galick_error_validate_vcf_unsupported_extension"
  )
})

test_that("validate_vcf reports missing fileformat as error", {
  path <- write_vcf_lines(c(
    "##contig=<ID=chr1>",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0"
  ))
  result <- validate_vcf(path)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "vcf_missing_fileformat"))
})

test_that("validate_vcf reports missing #CHROM as error and stops walking", {
  path <- write_vcf_lines(c("##fileformat=VCFv4.2", "chr1\t100\trs1\tA\tT"))
  result <- validate_vcf(path)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "vcf_missing_chrom_header"))
})

test_that("validate_vcf flags rows with wrong column count", {
  lines <- clean_vcf_lines()
  lines[7] <- "chr1\t200\trs2\tC\tG\t.\tPASS\t."
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "vcf_invalid_column_count"))
})

test_that("validate_vcf preserves invalid raw POS as a returned error", {
  lines <- clean_vcf_lines()
  lines[7] <- "chr1\tnotnumeric\trs2\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "vcf_invalid_pos"))
})

test_that("validate_vcf flags negative or zero POS", {
  lines <- clean_vcf_lines()
  lines[7] <- "chr1\t0\trs2\tC\tG\t.\tPASS\t.\tGT\t0/0\t0/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_invalid_pos"))
})

test_that("validate_vcf treats '.' marker IDs as missing (warning), not duplicate", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  lines[7] <- "chr1\t200\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$warnings$code == "vcf_missing_marker_ids"))
  expect_false(any(result$errors$code == "vcf_duplicated_marker_ids"))
})

test_that("validate_vcf flags duplicated non-missing marker IDs", {
  lines <- clean_vcf_lines()
  lines[7] <- "chr1\t200\trs1\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_duplicated_marker_ids"))
})

test_that("validate_vcf flags duplicated (chrom, pos)", {
  lines <- clean_vcf_lines()
  lines[7] <- "chr1\t100\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_duplicated_chrom_pos"))
})

test_that("validate_vcf flags multiallelic ALT as warning under general", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT,C\t.\tPASS\t.\tGT\t0/0\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path, target = "general")
  expect_true(result$ok)
  expect_true(any(result$warnings$code == "vcf_multiallelic"))
})

test_that("validate_vcf flags multiallelic ALT as error under target=beagle", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT,C\t.\tPASS\t.\tGT\t0/0\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path, target = "beagle")
  expect_false(result$ok)
  expect_true(any(result$errors$code == "vcf_multiallelic"))
})

test_that("validate_vcf flags REF == ALT under target=beagle", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tA\t.\tPASS\t.\tGT\t0/0\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path, target = "beagle")
  expect_true(any(result$errors$code == "vcf_beagle_ref_eq_alt"))
})

test_that("validate_vcf flags non-diploid GT under target=beagle", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path, target = "beagle")
  expect_true(any(result$errors$code == "vcf_beagle_non_diploid"))
})

test_that("validate_vcf flags GT in FORMAT=GT:DP combos correctly", {
  lines <- clean_vcf_lines()
  lines[5] <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2"
  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT:DP\t0/0:30\t1/1:25"
  lines[7] <- "chr1\t200\trs2\tC\tG\t.\tPASS\t.\tGT:DP\t0/1:15\t1/1:42"
  lines[8] <- "chr2\t50\trs3\tG\tA\t.\tPASS\t.\tGT:DP\t1/1:50\t0/1:33"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(result$ok)
})

test_that("validate_vcf rejects malformed GT tokens via strict regex", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\ta/b\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_invalid_gt"))
})

test_that("validate_vcf rejects ambiguous GT tokens like '..' and '0..'", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t..\t1/1"
  path <- write_vcf_lines(lines)
  expect_true(any(validate_vcf(path)$errors$code == "vcf_invalid_gt"))

  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT\t0..\t1/1"
  path2 <- write_vcf_lines(lines)
  expect_true(any(validate_vcf(path2)$errors$code == "vcf_invalid_gt"))
})

test_that("validate_vcf reports missing GT in FORMAT", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tDP\t30\t25"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_invalid_format_field"))
})

test_that("validate_vcf flags unsorted records with one warning per chrom", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t200\trs1\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  lines[7] <- "chr1\t100\trs2\tC\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$warnings$code == "vcf_unsorted_records"))
})

test_that("validate_vcf no_samples is warning under general, error under beagle", {
  lines <- c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    "chr1\t100\trs1\tA\tT\t.\tPASS\t.\tGT"
  )
  path <- write_vcf_lines(lines)
  general <- validate_vcf(path, target = "general")
  expect_true(any(general$warnings$code == "vcf_no_samples"))
  expect_false(any(general$errors$code == "vcf_no_samples"))

  beagle <- validate_vcf(path, target = "beagle")
  expect_true(any(beagle$errors$code == "vcf_no_samples"))
})

test_that("validate_vcf flags duplicated sample names", {
  lines <- clean_vcf_lines()
  lines[5] <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_duplicated_samples"))
})

test_that("validate_vcf target=plink flags missing IDs as error", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path, target = "plink")
  expect_true(any(result$errors$code == "vcf_plink_missing_ids"))
})

test_that("validate_vcf target=beagle warns when extension is plain .vcf", {
  path <- write_vcf_lines(clean_vcf_lines())
  result <- validate_vcf(path, target = "beagle")
  expect_true(any(result$warnings$code == "vcf_beagle_not_bgzf"))
})

test_that("validate_vcf max_records_to_scan limits scan and skips full-only checks", {
  lines <- c(
    clean_vcf_lines(),
    "chr2\t100\trs4\tT\tA\t.\tPASS\t.\tGT\t0/0\t1/1"
  )
  # Make rs1 + rs4 conflict on (chr1, 100) if all scanned
  lines[length(lines)] <- "chr1\t100\trs4\tT\tA\t.\tPASS\t.\tGT\t0/0\t1/1"
  path <- write_vcf_lines(lines)

  partial <- validate_vcf(path, max_records_to_scan = 2)
  expect_true(partial$summary$scanned_partial)
  expect_equal(partial$summary$n_records_scanned, 2L)
  expect_true(is.na(partial$summary$n_records))
  expect_false(any(partial$errors$code == "vcf_duplicated_chrom_pos"))

  full <- validate_vcf(path)
  expect_true(any(full$errors$code == "vcf_duplicated_chrom_pos"))
})

test_that("validate_vcf summary has_vcf_gz_extension is FALSE for .vcf input", {
  path <- write_vcf_lines(clean_vcf_lines())
  result <- validate_vcf(path)
  expect_false(result$summary$has_vcf_gz_extension)
})

test_that("validate_vcf reports invalid REF (non-DNA chars)", {
  lines <- clean_vcf_lines()
  lines[6] <- "chr1\t100\trs1\tQQ\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  path <- write_vcf_lines(lines)
  result <- validate_vcf(path)
  expect_true(any(result$errors$code == "vcf_invalid_ref_alt"))
})
