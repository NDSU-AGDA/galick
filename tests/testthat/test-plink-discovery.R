test_that("PLINK discovery rejects the PuTTY namesake", {
  executable <- withr::local_tempfile()
  file.create(executable)
  testthat::local_mocked_bindings(run_cli = function(...) {
    list(status = 0L, stdout = "plink: Release 0.83\nPuTTY Link", stderr = "")
  })
  expect_false(is_genetics_plink(executable))
  expect_error(find_plink(executable), class = "galick_error_plink_incompatible")
  expect_error(find_plink(executable), "PuTTY")
})

test_that("PLINK version probes handle supported banners and failed commands", {
  for (banner in c("PLINK v1.90b7.7 64-bit", "PLINK v2.00a5.12LM 64-bit")) {
    testthat::local_mocked_bindings(run_cli = function(...) {
      list(status = 0L, stdout = banner, stderr = "")
    })
    expect_true(is_genetics_plink("plink"))
  }
  testthat::local_mocked_bindings(run_cli = function(...) {
    list(status = 1L, stdout = "PLINK v1.90", stderr = "failed")
  })
  expect_false(is_genetics_plink("plink"))
  testthat::local_mocked_bindings(run_cli = function(...) stop("unavailable"))
  expect_false(is_genetics_plink("plink"))
})

test_that("automatic PLINK discovery refuses a PuTTY-only PATH", {
  testthat::local_mocked_bindings(is_genetics_plink = function(...) FALSE)
  expect_error(find_plink(), class = "galick_error_plink_not_found")
  expect_condition(skip_if_no_plink(), class = "skip")
})
