test_that("with_run_dir creates a temp dir when NULL and sets the option", {
  observed <- galick:::with_run_dir(NULL, getOption("galick.run_dir"))
  expect_match(observed, "galick_run_")
  expect_true(dir.exists(observed))
})

test_that("with_run_dir restores the previous option after evaluation", {
  prev <- getOption("galick.run_dir")
  galick:::with_run_dir(NULL, NULL)
  expect_identical(getOption("galick.run_dir"), prev)
})

test_that("with_run_dir restores the previous option even on error", {
  prev <- getOption("galick.run_dir")
  expect_error(
    galick:::with_run_dir(NULL, stop("boom"))
  )
  expect_identical(getOption("galick.run_dir"), prev)
})

test_that("with_run_dir uses the given dir when provided", {
  d <- withr::local_tempdir()
  observed <- galick:::with_run_dir(d, getOption("galick.run_dir"))
  expect_identical(observed, normalizePath(d, mustWork = TRUE))
})

test_that("with_run_dir creates the dir if it does not exist", {
  base <- withr::local_tempdir()
  newdir <- file.path(base, "subdir1", "subdir2")
  expect_false(dir.exists(newdir))
  observed <- galick:::with_run_dir(newdir, getOption("galick.run_dir"))
  expect_true(dir.exists(newdir))
})
