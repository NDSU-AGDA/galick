test_that("every exported function has a help example", {
  source_dir <- testthat::test_path("..", "..")
  help <- if (dir.exists(file.path(source_dir, "man"))) {
    tools::Rd_db(dir = source_dir)
  } else {
    tools::Rd_db("galick")
  }
  exports <- getNamespaceExports("galick")
  exports <- exports[!grepl(".", exports, fixed = TRUE)]
  example_tag <- paste0(intToUtf8(92L), "examples")
  for (name in exports) {
    page <- help[[paste0(name, ".Rd")]]
    expect_false(is.null(page), info = name)
    tags <- vapply(page, function(node) attr(node, "Rd_tag"), character(1))
    expect_true(example_tag %in% tags, info = name)
  }
})

test_that("instructional founders support validation and small simulations", {
  founders <- system.file("extdata", "tiny_founders.vcf", package = "galick")
  before <- digest::digest(file = founders, algo = "sha256")
  checked <- validate_vcf(founders, target = "beagle", progress = FALSE)
  expect_true(checked$ok)
  expect_identical(checked$summary$n_records, 12L)
  expect_identical(checked$summary$n_samples, 4L)
  map <- data.frame(chrom = "1", pos = seq(100, 1200, 100), cm = seq(0, 55, 5))
  output <- file.path(withr::local_tempdir(), "simulation")
  simulation <- simulate_vcf(founders, map, output, scenario = "f2",
    parents = c("f1", "f2"), n_target = 10L, n_reference = 0L, seed = 42L
  )
  expect_true(simulation$ok)
  expect_true(validate_vcf(simulation$paths[["truth"]], progress = FALSE)$ok)
  expect_identical(digest::digest(file = founders, algo = "sha256"), before)
})
