test_that("install_beagle rejects unsupported version with classed error", {
  expect_error(
    install_beagle(version = "5.5"),
    class = "beagler_error_beagle_unsupported_version"
  )
  expect_error(
    install_beagle(version = "4.2"),
    class = "beagler_error_beagle_unsupported_version"
  )
})

test_that("install_beagle aborts on invalid dest", {
  expect_error(
    install_beagle(dest = c("a", "b")),
    class = "beagler_error_invalid_argument"
  )
  expect_error(
    install_beagle(dest = NA_character_),
    class = "beagler_error_invalid_argument"
  )
})

skip_if_no_network <- function() {
  if (!identical(Sys.getenv("BEAGLER_TEST_NETWORK"), "true")) {
    testthat::skip("Set BEAGLER_TEST_NETWORK=true to run network-dependent tests.")
  }
}

test_that("install_beagle downloads to a tempdir and beagle_version returns 5.4", {
  skip_if_no_network()
  skip_if_no_java()
  dest_dir <- withr::local_tempdir()
  dest <- file.path(dest_dir, "beagle.jar")
  out <- install_beagle(dest = dest)
  expect_true(file.exists(dest))
  expect_match(beagle_version(beagle_jar = dest), "^5\\.")
})
