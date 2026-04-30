make_csv_data <- function() {
  data.frame(
    sample_id = c("s1", "s2", "s3"),
    rs1 = c("A", "H", "B"),
    rs2 = c("H", "B", "A"),
    rs3 = c("B", "A", "H"),
    stringsAsFactors = FALSE
  )
}

make_map_data <- function() {
  data.frame(
    marker_id = c("rs1", "rs2", "rs3"),
    chrom = c("chr1", "chr1", "chr2"),
    pos = c(100L, 200L, 50L),
    ref = c("A", "C", "G"),
    alt = c("T", "G", "A"),
    stringsAsFactors = FALSE
  )
}

test_that("csv_to_vcf happy path writes a valid VCF", {
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(make_csv_data(), make_map_data(), out)

  expect_s3_class(result, "beagler_result")
  expect_true(result$ok)
  expect_equal(result$summary$output_path, out)
  expect_true(file.exists(out))

  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$samples, c("s1", "s2", "s3"))
  expect_equal(parsed$chrom, c("chr1", "chr1", "chr2"))
  expect_equal(parsed$pos, c(100L, 200L, 50L))
  expect_equal(parsed$id, c("rs1", "rs2", "rs3"))
  expect_equal(parsed$ref, c("A", "C", "G"))
  expect_equal(parsed$alt, c("T", "G", "A"))
})

test_that("csv_to_vcf round-trips AHB cells to VCF GT", {
  out <- withr::local_tempfile(fileext = ".vcf")
  csv_to_vcf(make_csv_data(), make_map_data(), out)
  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$gt_matrix[1, ], c(s1 = "0/0", s2 = "0/1", s3 = "1/1"))
  expect_equal(parsed$gt_matrix[2, ], c(s1 = "0/1", s2 = "1/1", s3 = "0/0"))
  expect_equal(parsed$gt_matrix[3, ], c(s1 = "1/1", s2 = "0/0", s3 = "0/1"))
})

test_that("csv_to_vcf accepts paths for csv and marker_map", {
  csv_path <- withr::local_tempfile(fileext = ".csv")
  map_path <- withr::local_tempfile(fileext = ".tsv")
  vroom::vroom_write(make_csv_data(), csv_path, delim = ",")
  vroom::vroom_write(make_map_data(), map_path, delim = "\t")

  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(csv_path, map_path, out)
  expect_true(result$ok)
  expect_true(file.exists(out))
})

test_that("csv_to_vcf reports unmappable markers as a $errors and writes nothing", {
  csv <- make_csv_data()
  csv$rs_extra <- c("A", "H", "B")
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(csv, make_map_data(), out)

  expect_false(result$ok)
  expect_true(any(result$errors$code == "csv_to_vcf_unmappable_markers"))
  expect_false(file.exists(out))
  expect_true(is.na(result$summary$output_path))
})

test_that("csv_to_vcf reports invalid alleles and writes nothing", {
  bad_map <- make_map_data()
  bad_map$ref[1] <- "AT"  # not A/C/G/T
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(make_csv_data(), bad_map, out)

  expect_false(result$ok)
  expect_true(any(result$errors$code == "csv_to_vcf_invalid_alleles"))
  expect_false(file.exists(out))
})

test_that("csv_to_vcf rejects non-.vcf output extension and writes nothing", {
  out <- withr::local_tempfile(fileext = ".vcf.gz")
  result <- csv_to_vcf(make_csv_data(), make_map_data(), out)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "csv_to_vcf_unsupported_extension"))
  expect_false(file.exists(out))
})

test_that("csv_to_vcf refuses to overwrite by default", {
  out <- withr::local_tempfile(fileext = ".vcf")
  file.create(out)
  result <- csv_to_vcf(make_csv_data(), make_map_data(), out)
  expect_false(result$ok)
  expect_true(any(result$errors$code == "csv_to_vcf_output_exists"))
})

test_that("csv_to_vcf overwrites when overwrite = TRUE", {
  out <- withr::local_tempfile(fileext = ".vcf")
  file.create(out)
  result <- csv_to_vcf(make_csv_data(), make_map_data(), out, overwrite = TRUE)
  expect_true(result$ok)
  expect_true(file.exists(out))
  expect_gt(file.size(out), 0)
})

test_that("csv_to_vcf aggregates upstream errors from M2 and M3 (no write)", {
  bad_csv <- make_csv_data()
  bad_csv$sample_id <- c("s1", "s1", "s2")  # duplicate sample
  bad_map <- make_map_data()
  bad_map$pos[1] <- -1  # invalid pos

  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(bad_csv, bad_map, out, sample_col = "sample_id")

  expect_false(result$ok)
  codes <- result$errors$code
  expect_true(any(codes == "genotype_table_duplicated_samples"))
  expect_true(any(codes == "marker_map_invalid_pos"))
  expect_false(file.exists(out))
})

test_that("csv_to_vcf supports 012 encoding round-trip", {
  csv <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    rs1 = c("0", "1", "2"),
    rs2 = c("1", "2", "0"),
    stringsAsFactors = FALSE
  )
  map <- data.frame(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100L, 200L),
    ref = c("A", "C"),
    alt = c("T", "G"),
    stringsAsFactors = FALSE
  )
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(csv, map, out)
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$gt_matrix[1, ], c(s1 = "0/0", s2 = "0/1", s3 = "1/1"))
  expect_equal(parsed$gt_matrix[2, ], c(s1 = "0/1", s2 = "1/1", s3 = "0/0"))
})

test_that("csv_to_vcf supports IUPAC encoding round-trip", {
  csv <- data.frame(
    sample_id = c("s1", "s2", "s3"),
    rs1 = c("A", "W", "T"),  # ref=A alt=T: A->0/0, W->0/1, T->1/1
    rs2 = c("S", "G", "C"),  # ref=C alt=G: S->0/1, G->1/1, C->0/0
    stringsAsFactors = FALSE
  )
  map <- data.frame(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100L, 200L),
    ref = c("A", "C"),
    alt = c("T", "G"),
    stringsAsFactors = FALSE
  )
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(csv, map, out, encoding = "IUPAC")
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$gt_matrix[1, ], c(s1 = "0/0", s2 = "0/1", s3 = "1/1"))
  expect_equal(parsed$gt_matrix[2, ], c(s1 = "0/1", s2 = "1/1", s3 = "0/0"))
})

test_that("csv_to_vcf preserves missing values as ./.", {
  csv <- data.frame(
    sample_id = c("s1", "s2"),
    rs1 = c("A", NA),
    rs2 = c(".", "B"),
    stringsAsFactors = FALSE
  )
  map <- data.frame(
    marker_id = c("rs1", "rs2"),
    chrom = c("chr1", "chr1"),
    pos = c(100L, 200L),
    ref = c("A", "C"),
    alt = c("T", "G"),
    stringsAsFactors = FALSE
  )
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- csv_to_vcf(csv, map, out)
  expect_true(result$ok)
  parsed <- beagler:::read_vcf_minimal(out)
  expect_equal(parsed$gt_matrix[1, ], c(s1 = "0/0", s2 = "./."))
  expect_equal(parsed$gt_matrix[2, ], c(s1 = "./.", s2 = "1/1"))
})
