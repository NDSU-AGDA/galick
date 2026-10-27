test_that("Galick exposes the complete API under its new namespace", {
  expect_equal(unname(getNamespaceName(asNamespace("galick"))), "galick")
  expect_length(getNamespaceExports("galick")[!grepl("[.]", getNamespaceExports("galick"))], 28L)
  result <- validate_genotype_table(
    data.frame(sample_id = c("plant_1", "plant_2"), m1 = c(0, 2)),
    encoding = "012"
  )
  expect_s3_class(result, "galick_result")
  expect_match(format(result), "galick_result", fixed = TRUE)
  expect_equal(result$metadata$package_version, utils::packageVersion("galick"))
  expect_error(validate_genotype_table(NULL), class = "galick_error_invalid_input")
})

test_that("bundled conversion writes the Galick VCF source", {
  csv <- system.file("extdata", "tiny.csv", package = "galick")
  map <- system.file("extdata", "tiny_map.tsv", package = "galick")
  output <- withr::local_tempfile(fileext = ".vcf")
  expect_true(csv_to_vcf(csv, map, output)$ok)
  expect_true(any(startsWith(readLines(output), "##source=galick/")))
  expect_identical(dim(vcf_to_dosage(output)), c(5L, 20L))
})

test_that("Galick run-record serialization preserves JSON array and object types", {
  path <- write_run_record_json(
    withr::local_tempdir(),
    command = "java",
    versions = list(galick = "0.1.0", r = "4.5.3", os = "test"),
    runtime_sec = 0, status = 0L,
    started_at = "2026-10-27T12:00:00Z", ended_at = "2026-10-27T12:00:00Z"
  )
  text <- paste(readLines(path), collapse = "\n")
  expect_match(text, '"command": \\["java"\\]')
  expect_match(text, '"args": \\{\\}')
  record <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expect_equal(record$schema_version, "2.0.0")
  expect_named(record$versions, c("galick", "r", "os"))
})
