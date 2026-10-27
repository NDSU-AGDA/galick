rscript_bin <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

test_that("run_cli runs a successful command and captures stdout", {
  result <- galick:::run_cli(rscript_bin, c("-e", "cat('hello')"))
  expect_equal(result$status, 0L)
  expect_match(result$stdout, "hello")
  expect_true(result$runtime_sec >= 0)
  expect_length(result$command, 3L)
  expect_equal(result$command[2:3], c("-e", "cat('hello')"))
})

test_that("run_cli reports non-zero exit status without raising", {
  result <- galick:::run_cli(rscript_bin, c("-e", "quit(status = 7)"))
  expect_equal(result$status, 7L)
})

test_that("run_cli aborts with class when executable not found", {
  expect_error(
    galick:::run_cli("definitely_not_an_executable_xyzzy"),
    class = "galick_error_executable_not_found"
  )
})

test_that("run_cli respects timeout", {
  result <- galick:::run_cli(
    rscript_bin,
    c("-e", "Sys.sleep(5)"),
    timeout = 1
  )
  expect_equal(result$status, -1L)
  expect_match(result$stderr, "timeout")
})

test_that("run_cli writes captured stdout to a log file", {
  tmp <- withr::local_tempfile()
  galick:::run_cli(
    rscript_bin,
    c("-e", "cat('logged-stdout')"),
    stdout_log = tmp
  )
  expect_true(file.exists(tmp))
  content <- paste(readLines(tmp, warn = FALSE), collapse = "")
  expect_match(content, "logged-stdout")
})

test_that("run_cli writes captured stderr to a log file", {
  tmp <- withr::local_tempfile()
  galick:::run_cli(
    rscript_bin,
    c("-e", "writeLines('hi-from-stderr', stderr())"),
    stderr_log = tmp
  )
  expect_true(file.exists(tmp))
  content <- paste(readLines(tmp, warn = FALSE), collapse = "")
  expect_match(content, "hi-from-stderr")
})
