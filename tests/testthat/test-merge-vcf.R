write_merge_fixture <- function(path, sample, records, alt_header = NULL) {
  meta <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=1000>",
    if (!is.null(alt_header)) alt_header,
    paste(
      c(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "INFO", "FORMAT", sample
      ),
      collapse = "\t"
    )
  )
  writeLines(c(meta, records), path)
  path
}

test_that("by_sample reindexes ALT calls and retains duplicate records", {
  a <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"),
    "a",
    c(
      "chr1\t10\tdup1\tA\tG\t.\tPASS\t.\tGT\t0/1",
      "chr1\t10\tdup2\tA\tG\t.\tPASS\t.\tGT\t1/1",
      "chr1\t20\tbadref\tA\tC\t.\tPASS\t.\tGT\t0/1",
      "chr1\t30\treindex\tC\tG\t.\tPASS\t.\tGT\t1/1"
    ),
    "##first-only=retained"
  )
  b <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"),
    "b",
    c(
      "chr1\t10\tdup1\tA\tG\t.\tPASS\t.\tGT\t1/1",
      "chr1\t10\tdup2\tA\tG\t.\tPASS\t.\tGT\t0/0",
      "chr1\t20\tbadref\tG\tC\t.\tPASS\t.\tGT\t0/1",
      "chr1\t30\treindex\tC\tT,G\t.\tPASS\t.\tGT\t2|1"
    )
  )
  output <- withr::local_tempfile(fileext = ".vcf")

  report <- merge_vcf(c(a, b), output, mode = "by_sample")
  parsed <- galick:::read_genotype_vcf(output)

  expect_equal(nrow(parsed$fixed), 3L)
  expect_equal(parsed$fixed$ID, c("dup1", "dup2", "reindex"))
  expect_equal(parsed$gt[1:2, "a"], c("0/1", "1/1"))
  expect_equal(parsed$gt[1:2, "b"], c("1/1", "0/0"))
  expect_equal(parsed$fixed$ALT[3], "G,T")
  expect_equal(unname(parsed$gt[3, "b"]), "1|2")
  expect_equal(unname(report$summary$sites_in), c(4L, 4L))
  expect_equal(report$summary$sites_out, 3L)
  expect_true(any(report$summary$sites_dropped$reason == "REF mismatch"))
  expect_true("##first-only=retained" %in% readLines(output))
})

test_that("by_sample handles sample collisions only by named policy", {
  a <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "s",
    "chr1\t1\tx\tA\tG\t.\tPASS\t.\tGT\t1/1"
  )
  b <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "s",
    "chr1\t1\tx\tA\tT,G\t.\tPASS\t.\tGT\t2/2"
  )

  expect_error(
    merge_vcf(c(a, b), withr::local_tempfile(fileext = ".vcf"),
      mode = "by_sample"
    ),
    class = "galick_error_merge_sample_collision"
  )

  suffix_output <- withr::local_tempfile(fileext = ".vcf")
  suffix_report <- merge_vcf(
    c(a, b), suffix_output,
    mode = "by_sample", sample_collision = "suffix"
  )
  expect_equal(
    galick:::read_genotype_vcf(suffix_output)$samples,
    c("s", "s.1")
  )
  expect_equal(suffix_report$summary$collisions$resolution, "suffix")

  coalesce_output <- withr::local_tempfile(fileext = ".vcf")
  coalesce_report <- merge_vcf(
    c(a, b), coalesce_output,
    mode = "by_sample",
    sample_collision = "coalesce"
  )
  coalesced <- galick:::read_genotype_vcf(coalesce_output)
  expect_equal(coalesced$samples, "s")
  expect_equal(unname(coalesced$gt[, "s"]), "1/1")
  expect_equal(coalesce_report$summary$samples_out, 1L)

  conflicting <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "s",
    "chr1\t1\tx\tA\tG\t.\tPASS\t.\tGT\t0/0"
  )
  expect_error(
    merge_vcf(
      c(a, conflicting), withr::local_tempfile(fileext = ".vcf"),
      mode = "by_sample", sample_collision = "coalesce"
    ),
    class = "galick_error_merge_coalesce_conflict"
  )
})

test_that("by_marker row-binds different markers and aligns sample columns", {
  a <- withr::local_tempfile(fileext = ".vcf")
  b <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2",
    "chr1\t20\tb\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1"
  ), a)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts2\ts1",
    "chr1\t10\ta\tC\tT\t.\tPASS\t.\tGT\t0/0\t1/1"
  ), b)
  output <- withr::local_tempfile(fileext = ".vcf")

  report <- merge_vcf(c(a, b), output, mode = "by_marker")
  parsed <- galick:::read_genotype_vcf(output)
  expect_equal(parsed$fixed$ID, c("a", "b"))
  expect_equal(parsed$samples, c("s1", "s2"))
  expect_equal(
    unname(parsed$gt[parsed$fixed$ID == "a", ]),
    c("1/1", "0/0")
  )
  expect_equal(report$summary$sites_out, 2L)
  expect_equal(report$summary$samples_out, 2L)
})

test_that("merge drops unresolvable calls and reports the reason", {
  a <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "a",
    c(
      "chr1\t1\tbad\tA\tG\t.\tPASS\t.\tGT\t2/2",
      "chr1\t2\tgood\tA\tG\t.\tPASS\t.\tGT\t0/1"
    )
  )
  b <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "b",
    c(
      "chr1\t1\tbad\tA\tG\t.\tPASS\t.\tGT\t0/1",
      "chr1\t2\tgood\tA\tG\t.\tPASS\t.\tGT\t1/1"
    )
  )
  output <- withr::local_tempfile(fileext = ".vcf")

  report <- merge_vcf(c(a, b), output, mode = "by_sample")
  expect_equal(galick:::read_genotype_vcf(output)$fixed$ID, "good")
  expect_true(any(
    report$summary$sites_dropped$reason == "unresolvable alleles"
  ))
})

test_that("by_marker refuses the same canonical record across files", {
  a <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "s",
    "chr1\t1\tx\tA\tG,T\t.\tPASS\t.\tGT\t1/2"
  )
  b <- write_merge_fixture(
    withr::local_tempfile(fileext = ".vcf"), "s",
    "chr1\t1\tx\tA\tT,G\t.\tPASS\t.\tGT\t2/1"
  )

  expect_error(
    merge_vcf(
      c(a, b),
      withr::local_tempfile(fileext = ".vcf"),
      mode = "by_marker"
    ),
    class = "galick_error_merge_marker_contract"
  )
})
