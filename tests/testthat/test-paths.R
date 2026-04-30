test_that("normalize_path_safe expands ~ to an absolute path", {
  home <- withr::local_tempdir()
  withr::local_envvar(HOME = home)
  p <- beagler:::normalize_path_safe("~/foo")
  expect_false(startsWith(p, "~"))
  expect_match(p, "foo$")
})

test_that("normalize_path_safe rejects non-character or NA input", {
  expect_error(
    beagler:::normalize_path_safe(NA_character_),
    class = "beagler_error_invalid_path"
  )
  expect_error(
    beagler:::normalize_path_safe(c("a", "b")),
    class = "beagler_error_invalid_path"
  )
  expect_error(
    beagler:::normalize_path_safe(123),
    class = "beagler_error_invalid_path"
  )
})

test_that("normalize_path_safe with must_exist = TRUE errors on missing file", {
  missing <- file.path(tempdir(), "definitely-not-here-xyzzy")
  expect_error(beagler:::normalize_path_safe(missing, must_exist = TRUE))
})

test_that("assert_readable returns the absolute path invisibly", {
  f <- withr::local_tempfile()
  writeLines("hi", f)
  result <- beagler:::assert_readable(f)
  expect_match(result, "[Tt]mp")
})

test_that("assert_readable errors on missing file", {
  expect_error(
    beagler:::assert_readable(file.path(tempdir(), "missing-xyzzy"))
  )
})

test_that("assert_writable accepts a writable parent dir for non-existing file", {
  d <- withr::local_tempdir()
  expect_no_error(beagler:::assert_writable(file.path(d, "newfile.txt")))
})

test_that("assert_writable accepts an existing writable file", {
  f <- withr::local_tempfile()
  writeLines("hi", f)
  expect_no_error(beagler:::assert_writable(f))
})

test_that("assert_writable errors when parent dir does not exist", {
  expect_error(
    beagler:::assert_writable("/nonexistent_xyz_dir/file.txt"),
    class = "beagler_error_path_not_writable"
  )
})
