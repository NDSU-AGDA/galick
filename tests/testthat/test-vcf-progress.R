test_that("progress does not change validation results", {
  path <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
  withr::local_options(cli.progress_handlers = character())
  for (target in c("general", "beagle", "plink")) {
    for (limit in c(2, Inf)) {
      quiet <- validate_vcf(path, target, limit, progress = FALSE)
      shown <- validate_vcf(path, target, limit, progress = TRUE)
      expect_identical(shown$errors, quiet$errors)
      expect_identical(shown$warnings, quiet$warnings)
      expect_identical(shown$summary, quiet$summary)
    }
  }
})

test_that("progress stages render and disabled progress stays silent", {
  path <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
  withr::local_options(list(
    cli.progress_handlers_only = "cli", cli.progress_show_after = 0,
    cli.progress_clear = FALSE, cli.num_colors = 1, width = 140
  ))
  output <- capture.output(invisible(validate_vcf(path, target = "beagle", progress = TRUE)), type = "message")
  expect_true(any(grepl("Reading and parsing VCF", output, fixed = TRUE)))
  expect_true(any(grepl("Checking VCF structure", output, fixed = TRUE)))
  expect_true(any(grepl("Genotypes", output, fixed = TRUE)))
  expect_true(any(grepl("Beagle checks", output, fixed = TRUE)))
  expect_true(any(grepl("100%", output, fixed = TRUE)))
  expect_true(any(grepl("ETA", output, fixed = TRUE)))
  expect_length(capture.output(invisible(validate_vcf(path, progress = FALSE)), type = "message"), 0)
})

test_that("progress is opt-in outside interactive sessions and validates its argument", {
  expect_identical(formals(validate_vcf)$progress, quote(interactive()))
  path <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
  for (value in list(NULL, NA, 1, "yes", c(TRUE, FALSE))) {
    expect_error(validate_vcf(path, progress = value), class = "galick_error_invalid_argument")
  }
})

test_that("progress bars close after parser errors and early returns", {
  withr::local_options(cli.progress_handlers = character())
  state <- get("clienv", envir = asNamespace("cli"))
  before <- names(state$progress)
  expect_error(validate_vcf(tempfile(fileext = ".vcf"), progress = TRUE), "File not found")
  expect_identical(names(state$progress), before)
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines("##fileformat=VCFv4.2", path)
  expect_false(validate_vcf(path, progress = TRUE)$ok)
  expect_identical(names(state$progress), before)
})
