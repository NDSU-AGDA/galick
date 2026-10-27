test_that("log_step appends a JSONL entry to run.log", {
  d <- withr::local_tempdir()
  galick:::log_step(d, "validate_vcf", target = "beagle")
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
  galick:::log_step(d, "step1")
  galick:::log_step(d, "step2", extra = 42)
  lines <- readLines(file.path(d, "run.log"))
  expect_length(lines, 2L)
})

test_that("log_step errors if run_dir does not exist", {
  expect_error(
    galick:::log_step("/tmp/nonexistent_xyzzy_logging", "anything"),
    class = "galick_error_invalid_argument"
  )
})

test_that("log_step rejects non-character step", {
  d <- withr::local_tempdir()
  expect_error(
    galick:::log_step(d, NA_character_),
    class = "galick_error_invalid_argument"
  )
})

test_that("write_run_record_json writes a v2.0.0 schema record", {
  d <- withr::local_tempdir()
  out <- galick:::write_run_record_json(
    run_dir = d,
    command = c("Rscript", "-e", "1+1"),
    versions = list(galick = "0.0.0.9000", r = "4.5.3", os = "darwin"),
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
  expect_equal(parsed$schema_version, "2.0.0")
  expect_equal(unlist(parsed$command), c("Rscript", "-e", "1+1"))
  expect_equal(parsed$status, 0L)
})

test_that("write_run_record_json rejects missing required versions fields", {
  d <- withr::local_tempdir()
  expect_error(
    galick:::write_run_record_json(
      run_dir = d,
      command = c("a"),
      versions = list(galick = "x"),
      runtime_sec = 0,
      status = 0L,
      started_at = "2026-04-28T15:00:00Z",
      ended_at = "2026-04-28T15:00:00Z"
    ),
    class = "galick_error_invalid_run_record"
  )
})

test_that("write_run_record_json rejects negative runtime", {
  d <- withr::local_tempdir()
  expect_error(
    galick:::write_run_record_json(
      run_dir = d,
      command = c("a"),
      versions = list(galick = "x", r = "y", os = "z"),
      runtime_sec = -1,
      status = 0L,
      started_at = "2026-04-28T15:00:00Z",
      ended_at = "2026-04-28T15:00:00Z"
    ),
    class = "galick_error_invalid_run_record"
  )
})

test_that("write_run_record_json rejects empty command", {
  d <- withr::local_tempdir()
  expect_error(
    galick:::write_run_record_json(
      run_dir = d,
      command = character(),
      versions = list(galick = "x", r = "y", os = "z"),
      runtime_sec = 0,
      status = 0L,
      started_at = "2026-04-28T15:00:00Z",
      ended_at = "2026-04-28T15:00:00Z"
    ),
    class = "galick_error_invalid_run_record"
  )
})

test_that("run records reject invalid scalar fields before writing", {
  arguments <- list(
    run_dir = withr::local_tempdir(), command = "java",
    versions = list(galick = "0.1.0", r = "4.5.3", os = "test"),
    runtime_sec = 0, status = 0L,
    started_at = "2026-10-27T12:00:00Z", ended_at = "2026-10-27T12:00:00Z"
  )
  cases <- list(
    list(runtime_sec = NA_real_), list(runtime_sec = Inf),
    list(runtime_sec = c(0, 1)), list(status = 0.5), list(status = NA_integer_),
    list(status = c(0L, 1L)), list(command = NA_character_),
    list(versions = list(galick = NA_character_, r = "4.5.3", os = "test")),
    list(started_at = "yesterday"), list(ended_at = "2026-10-26T12:00:00Z"),
    list(args = list("unnamed")), list(inputs = list(list(path = "input.vcf"))),
    list(outputs = list(list(path = "output.vcf", sha256 = "invalid")))
  )
  for (case in cases) {
    expect_error(do.call(write_run_record_json, utils::modifyList(arguments, case)),
      class = "galick_error_invalid_run_record"
    )
  }
  expect_false(file.exists(file.path(arguments$run_dir, "run_record.json")))
})
