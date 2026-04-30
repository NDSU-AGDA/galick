test_that("with_run_dir creates a temp dir when NULL and sets the option", {
  observed <- beagler:::with_run_dir(NULL, getOption("beagler.run_dir"))
  expect_match(observed, "beagler_run_")
  expect_true(dir.exists(observed))
})

test_that("with_run_dir restores the previous option after evaluation", {
  prev <- getOption("beagler.run_dir")
  beagler:::with_run_dir(NULL, NULL)
  expect_identical(getOption("beagler.run_dir"), prev)
})

test_that("with_run_dir restores the previous option even on error", {
  prev <- getOption("beagler.run_dir")
  expect_error(
    beagler:::with_run_dir(NULL, stop("boom"))
  )
  expect_identical(getOption("beagler.run_dir"), prev)
})

test_that("with_run_dir uses the given dir when provided", {
  d <- withr::local_tempdir()
  observed <- beagler:::with_run_dir(d, getOption("beagler.run_dir"))
  expect_identical(observed, normalizePath(d, mustWork = TRUE))
})

test_that("with_run_dir creates the dir if it does not exist", {
  base <- withr::local_tempdir()
  newdir <- file.path(base, "subdir1", "subdir2")
  expect_false(dir.exists(newdir))
  observed <- beagler:::with_run_dir(newdir, getOption("beagler.run_dir"))
  expect_true(dir.exists(newdir))
})
