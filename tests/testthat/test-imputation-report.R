make_fake_run <- function(run_dir = NULL) {
  if (is.null(run_dir)) {
    run_dir <- withr::local_tempdir(.local_envir = parent.frame())
  }
  stdout_log <- file.path(run_dir, "stdout.log")
  stderr_log <- file.path(run_dir, "stderr.log")
  writeLines("ok", stdout_log)
  writeLines("", stderr_log)

  rec_path <- beagler:::write_run_record_json(
    run_dir = run_dir,
    command = c("/usr/bin/java", "-jar", "beagle.jar", "gt=in.vcf", "out=out"),
    versions = list(beagler = "0.1.0", r = "4.5.3", os = "darwin",
                    beagle = "5.4", java = "17"),
    args = list(gt = "in.vcf", out = "out", mode = "phase"),
    inputs = list(list(path = "in.vcf",
                       sha256 = paste0(rep("a", 64), collapse = ""))),
    outputs = list(list(path = "out.vcf.gz",
                        sha256 = paste0(rep("b", 64), collapse = ""))),
    runtime_sec = 12.3,
    status = 0L,
    started_at = "2026-04-28T15:00:00Z",
    ended_at = "2026-04-28T15:00:12Z"
  )

  beagler:::new_beagler_run(
    status = 0L,
    runtime_sec = 12.3,
    command = c("/usr/bin/java", "-jar", "beagle.jar", "gt=in.vcf", "out=out"),
    version = list(beagler = "0.1.0", r = "4.5.3", os = "darwin",
                   beagle = "5.4", java = "17"),
    output_paths = "out.vcf.gz",
    stdout_log = stdout_log,
    stderr_log = stderr_log,
    run_record = rec_path
  )
}

test_that("imputation_report writes a Markdown file with status, command, inputs", {
  run <- make_fake_run()
  out <- withr::local_tempfile(fileext = ".md")
  imputation_report(run, out)
  text <- readLines(out)
  expect_match(paste(text, collapse = "\n"), "# beagler run report")
  expect_match(paste(text, collapse = "\n"), "OK \\(0\\)")
  expect_match(paste(text, collapse = "\n"), "12\\.30 sec")
  expect_match(paste(text, collapse = "\n"), "beagle\\.jar gt=in\\.vcf")
  expect_match(paste(text, collapse = "\n"), "in\\.vcf")
  expect_match(paste(text, collapse = "\n"), "out\\.vcf\\.gz")
})

test_that("imputation_report includes the version block", {
  run <- make_fake_run()
  out <- withr::local_tempfile(fileext = ".md")
  imputation_report(run, out)
  text <- paste(readLines(out), collapse = "\n")
  expect_match(text, "beagler\\*\\*: 0\\.1\\.0")
  expect_match(text, "beagle\\*\\*: 5\\.4")
})

test_that("imputation_report aborts on non-beagler_run input", {
  expect_error(
    imputation_report(list(), tempfile()),
    class = "beagler_error_imputation_report_invalid_input"
  )
})

test_that("imputation_report aborts on bad file argument", {
  run <- make_fake_run()
  expect_error(
    imputation_report(run, NA_character_),
    class = "beagler_error_imputation_report_invalid_input"
  )
  expect_error(
    imputation_report(run, c("a", "b")),
    class = "beagler_error_imputation_report_invalid_input"
  )
})

test_that("imputation_report tolerates a run with NA status", {
  run <- beagler:::new_beagler_run()
  out <- withr::local_tempfile(fileext = ".md")
  imputation_report(run, out)
  text <- paste(readLines(out), collapse = "\n")
  expect_match(text, "status NA")
})
