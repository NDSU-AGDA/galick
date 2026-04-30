test_that("log_step appends a JSONL entry to run.log", {
  d <- withr::local_tempdir()
  beagler:::log_step(d, "validate_vcf", target = "beagle")
  log_path <- file.path(d, "run.log")
  expect_true(file.exists(log_path))

  lines <- readLines(log_path)
  expect_length(lines, 1L)
  parsed <- jsonlite::fromJSON(lines[1L], simplifyVector = FALSE)
  expect_equal(parsed$step, "validate_vcf")
  expect_equal(parsed$target, "beagle")
  expect_match(parsed$time, "^\\d{4}-\\d{2}-\\d{2}T")
})

test_that("log_step appends multiple entries (one per call)", {
  d <- withr::local_tempdir()
  beagler:::log_step(d, "step1")
  beagler:::log_step(d, "step2", extra = 42)
  lines <- readLines(file.path(d, "run.log"))
  expect_length(lines, 2L)
})

test_that("log_step errors if run_dir does not exist", {
  expect_error(
    beagler:::log_step("/tmp/nonexistent_xyzzy_logging", "anything"),
    class = "beagler_error_invalid_argument"
  )
})

test_that("log_step rejects non-character step", {
  d <- withr::local_tempdir()
  expect_error(
    beagler:::log_step(d, NA_character_),
    class = "beagler_error_invalid_argument"
  )
})

test_that("write_run_record_json writes a v1.0.0 schema record", {
  d <- withr::local_tempdir()
  out <- beagler:::write_run_record_json(
    run_dir = d,
    command = c("Rscript", "-e", "1+1"),
    versions = list(beagler = "0.0.0.9000", r = "4.5.3", os = "darwin"),
    args = list(gt = "input.vcf"),
    inputs = list(),
    outputs = list(),
    runtime_sec = 0.05,
    status = 0L,
    started_at = "2026-04-28T15:00:00Z",
    ended_at = "2026-04-28T15:00:00Z"
  )
  expect_true(file.exists(out))
  parsed <- jsonlite::fromJSON(out, simplifyVector = FALSE)
  expect_equal(parsed$schema_version, "1.0.0")
  expect_equal(unlist(parsed$command), c("Rscript", "-e", "1+1"))
  expect_equal(parsed$status, 0L)
})

test_that("write_run_record_json rejects missing required versions fields", {
  d <- withr::local_tempdir()
  expect_error(
    beagler:::write_run_record_json(
      run_dir = d,
      command = c("a"),
      versions = list(beagler = "x"),
      runtime_sec = 0,
      status = 0L,
      started_at = "2026-04-28T15:00:00Z",
      ended_at = "2026-04-28T15:00:00Z"
    ),
    class = "beagler_error_invalid_run_record"
  )
})

test_that("write_run_record_json rejects negative runtime", {
  d <- withr::local_tempdir()
  expect_error(
    beagler:::write_run_record_json(
      run_dir = d,
      command = c("a"),
      versions = list(beagler = "x", r = "y", os = "z"),
      runtime_sec = -1,
      status = 0L,
      started_at = "2026-04-28T15:00:00Z",
      ended_at = "2026-04-28T15:00:00Z"
    ),
    class = "beagler_error_invalid_run_record"
  )
})

test_that("write_run_record_json rejects empty command", {
  d <- withr::local_tempdir()
  expect_error(
    beagler:::write_run_record_json(
      run_dir = d,
      command = character(),
      versions = list(beagler = "x", r = "y", os = "z"),
      runtime_sec = 0,
      status = 0L,
      started_at = "2026-04-28T15:00:00Z",
      ended_at = "2026-04-28T15:00:00Z"
    ),
    class = "beagler_error_invalid_run_record"
  )
})
