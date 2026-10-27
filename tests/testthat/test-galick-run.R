test_that("new_galick_run defaults are coerced correctly", {
  r <- galick:::new_galick_run()
  expect_s3_class(r, "galick_run")
  expect_identical(r$status, NA_integer_)
  expect_identical(r$runtime_sec, NA_real_)
})

test_that("new_galick_run accepts a successful run", {
  r <- galick:::new_galick_run(
    status = 0L,
    runtime_sec = 1.23,
    command = c("Rscript", "-e", "1+1"),
    output_paths = "/tmp/output.vcf"
  )
  expect_equal(r$status, 0L)
  expect_equal(r$runtime_sec, 1.23)
  expect_length(r$command, 3L)
})

test_that("new_galick_run metadata has package_version and time", {
  r <- galick:::new_galick_run()
  expect_s3_class(r$metadata$package_version, "package_version")
  expect_match(r$metadata$time, "^\\d{4}-\\d{2}-\\d{2}T")
})

test_that("new_galick_run captures call into metadata", {
  call <- quote(run_beagle(gt = "x"))
  r <- galick:::new_galick_run(call = call)
  expect_identical(r$metadata$call, call)
})

test_that("print.galick_run snapshot for successful run", {
  testthat::local_reproducible_output()
  r <- galick:::new_galick_run(
    status = 0L,
    runtime_sec = 3.45,
    command = c("java", "-jar", "beagle.jar", "gt=in.vcf", "out=out"),
    output_paths = "/tmp/out.vcf.gz",
    stdout_log = "/tmp/stdout.log",
    stderr_log = "/tmp/stderr.log"
  )
  expect_snapshot(print(r))
})

test_that("print.galick_run snapshot for failed run", {
  testthat::local_reproducible_output()
  r <- galick:::new_galick_run(
    status = 7L,
    runtime_sec = 0.5,
    command = c("java", "-jar", "beagle.jar"),
    stderr_log = "/tmp/stderr.log"
  )
  expect_snapshot(print(r))
})

test_that("as.data.frame.galick_run produces a single row", {
  r <- galick:::new_galick_run(status = 0L, runtime_sec = 1.0)
  df <- as.data.frame(r)
  expect_equal(nrow(df), 1L)
  expect_equal(df$status, 0L)
})

test_that("format.galick_run returns a one-line summary", {
  r <- galick:::new_galick_run(status = 0L, runtime_sec = 1.5)
  expect_match(format(r), "status=0")
  expect_match(format(r), "runtime=1.50s")
})
