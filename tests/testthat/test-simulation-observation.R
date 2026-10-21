test_that("missingness calibration, heterogeneity and errors have explicit denominators", {
  withr::local_seed(1001)
  fixed <- data.table::data.table(
    CHROM = "1", POS = 1:200, ID = paste0("m", 1:200),
    REF = "A", ALT = "G", QUAL = ".", FILTER = "PASS", INFO = "."
  )
  x <- new_genotype_vcf(fixed, matrix("0|1", 200, 200), paste0("s", 1:200), character())
  settings <- list(
    keep = rep(TRUE, 200), missingness = list(rate = 0.3, sample_sd = 1, marker_sd = 1),
    error_rate = 0.2
  )
  observed <- simulation_observe(x, settings, 2L)
  expect_equal(observed$summary$expected_missing_rate, 0.3, tolerance = 1e-8)
  expect_lt(abs(observed$summary$realized_missing_rate - 0.3), 0.015)
  expect_lt(abs(observed$summary$realized_error_rate - 0.2), 0.015)
  expect_gt(stats::sd(observed$by_sample$realized_missing_rate), 0.05)
  expect_gt(stats::sd(observed$by_marker$realized_missing_rate), 0.05)
  expect_equal(nrow(observed$events), observed$summary$n_errors + observed$summary$n_missing)
  expect_true(all(observed$events$observed[observed$events$event == "error"] %in% c("0/0", "1/1")))
  expect_true(all(x$gt == "0|1"))
})

test_that("publication errors do not leave a partial destination", {
  f <- simulation_fixture()
  destination <- file.path(withr::local_tempdir(), "bundle")
  local_mocked_bindings(write_genotype_vcf = function(...) stop("Write failed"))
  expect_error(
    simulate_vcf(f$path, f$map, destination, n_target = 2, n_reference = 0, population_size = 4, generations = 1),
    "Write failed"
  )
  expect_false(dir.exists(destination))
  expect_length(list.files(dirname(destination)), 0)
})
