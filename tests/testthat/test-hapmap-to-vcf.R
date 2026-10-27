tiny_hapmap_text <- function() {
  paste(
    c(
      "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\ts1\ts2\ts3",
      "rs1\tA/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\tW\tT",
      "rs2\tC/G\tchr1\t200\t+\tNA\tNA\tNA\tNA\tNA\tNA\tS\tG\tC",
      "rs3\tG/A\tchr2\t50\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\tG\tR"
    ),
    collapse = "\n"
  )
}

write_tiny_hapmap <- function() {
  path <- withr::local_tempfile(fileext = ".hmp.txt", .local_envir = parent.frame())
  writeLines(tiny_hapmap_text(), path)
  path
}

test_that("read_hapmap parses the metadata, samples, and gt_matrix", {
  path <- write_tiny_hapmap()
  parsed <- galick:::read_hapmap(path)
  expect_equal(parsed$samples, c("s1", "s2", "s3"))
  expect_equal(nrow(parsed$gt_matrix), 3L)
  expect_equal(ncol(parsed$gt_matrix), 3L)
  expect_equal(parsed$metadata$marker_id, c("rs1", "rs2", "rs3"))
  expect_equal(parsed$metadata$chrom, c("chr1", "chr1", "chr2"))
  expect_equal(parsed$metadata$pos, c(100L, 200L, 50L))
  expect_equal(parsed$metadata$ref, c("A", "C", "G"))
  expect_equal(parsed$metadata$alt, c("T", "G", "A"))
})

test_that("read_hapmap aborts on missing required columns", {
  path <- withr::local_tempfile(fileext = ".hmp.txt")
  writeLines("foo\tbar\tbaz\n1\t2\t3", path)
  expect_error(
    galick:::read_hapmap(path),
    class = "galick_error_hapmap_invalid_input"
  )
})

test_that("read_hapmap aborts when no sample columns remain", {
  path <- withr::local_tempfile(fileext = ".hmp.txt")
  writeLines(
    paste(
      c(
        "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode",
        "rs1\tA/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA"
      ),
      collapse = "\n"
    ),
    path
  )
  expect_error(
    galick:::read_hapmap(path),
    class = "galick_error_hapmap_invalid_input"
  )
})

test_that("hapmap_to_vcf happy path writes a valid VCF with decoded IUPAC", {
  path <- write_tiny_hapmap()
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- hapmap_to_vcf(path, out)

  expect_s3_class(result, "galick_result")
  expect_true(result$ok)
  expect_equal(result$summary$output_path, out)

  parsed <- galick:::read_vcf_minimal(out)
  expect_equal(parsed$samples, c("s1", "s2", "s3"))
  # rs1 ref=A alt=T: A->0/0, W->0/1, T->1/1
  expect_equal(parsed$gt_matrix[1, ], c(s1 = "0/0", s2 = "0/1", s3 = "1/1"))
  # rs2 ref=C alt=G: S->0/1, G->1/1, C->0/0
  expect_equal(parsed$gt_matrix[2, ], c(s1 = "0/1", s2 = "1/1", s3 = "0/0"))
  # rs3 ref=G alt=A: A->1/1, G->0/0, R->0/1
  expect_equal(parsed$gt_matrix[3, ], c(s1 = "1/1", s2 = "0/0", s3 = "0/1"))
})

test_that("hapmap_to_vcf writes ##reference= when build is provided", {
  path <- write_tiny_hapmap()
  out <- withr::local_tempfile(fileext = ".vcf")
  hapmap_to_vcf(path, out, build = "GRCh38")
  parsed <- galick:::read_vcf_minimal(out)
  expect_true(any(grepl("^##reference=GRCh38$", parsed$meta)))
})

test_that("hapmap_to_vcf rejects non-.vcf output extension and writes nothing", {
  path <- write_tiny_hapmap()
  out <- withr::local_tempfile(fileext = ".vcf.gz")
  result <- hapmap_to_vcf(path, out)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "hapmap_to_vcf_unsupported_extension"))
  expect_false(file.exists(out))
})

test_that("hapmap_to_vcf refuses to overwrite by default", {
  path <- write_tiny_hapmap()
  out <- withr::local_tempfile(fileext = ".vcf")
  file.create(out)
  result <- hapmap_to_vcf(path, out)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "hapmap_to_vcf_output_exists"))
})

test_that("hapmap_to_vcf flags non-A/C/G/T alleles", {
  bad_text <- paste(
    c(
      "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\ts1",
      "rs1\tAT/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA"
    ),
    collapse = "\n"
  )
  path <- withr::local_tempfile(fileext = ".hmp.txt")
  writeLines(bad_text, path)
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- hapmap_to_vcf(path, out)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "hapmap_to_vcf_invalid_alleles"))
})

test_that("hapmap_to_vcf flags IUPAC tokens inconsistent with REF/ALT", {
  # rs1 ref=A alt=T but cell s1 uses 'Y' (C/T het) which is inconsistent
  bad_text <- paste(
    c(
      "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\ts1\ts2",
      "rs1\tA/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tY\tA"
    ),
    collapse = "\n"
  )
  path <- withr::local_tempfile(fileext = ".hmp.txt")
  writeLines(bad_text, path)
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- hapmap_to_vcf(path, out)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "hapmap_to_vcf_invalid_iupac_token"))
  expect_false(file.exists(out))
})

test_that("hapmap_to_vcf accepts N as missing", {
  text_with_n <- paste(
    c(
      "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\ts1\ts2",
      "rs1\tA/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tN\tA"
    ),
    collapse = "\n"
  )
  path <- withr::local_tempfile(fileext = ".hmp.txt")
  writeLines(text_with_n, path)
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- hapmap_to_vcf(path, out)
  expect_true(result$ok)
  parsed <- galick:::read_vcf_minimal(out)
  expect_equal(parsed$gt_matrix[1, ], c(s1 = "./.", s2 = "0/0"))
})
