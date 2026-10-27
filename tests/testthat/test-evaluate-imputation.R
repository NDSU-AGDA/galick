test_that("evaluation selects on tuning and confirms using isolated test jobs", {
  f <- evaluation_fixture()
  observed <- list()
  observe <- function(input, seed, directory) {
    observed[[basename(directory)]] <<- list(
      samples = input$samples, seed = seed,
      gt = input$gt, frozen = file.exists(file.path(dirname(dirname(directory)), "selection.json"))
    )
  }
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth, observer = observe))
  result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    fraction = 0.25, repeats = 2L
  )
  expect_s3_class(result, "galick_evaluation")
  expect_s3_class(result, "galick_result")
  expect_true(result$ok)
  expect_equal(result$selection$decision, "confirmed")
  expect_equal(result$selection$recommended_candidate, "good")
  expect_length(result$runs, 8)
  expect_identical(observed$base_tune_1$gt, observed$good_tune_1$gt)
  expect_identical(observed$base_tune_1$seed, observed$good_tune_1$seed)
  expect_false(observed$base_tune_1$frozen)
  expect_true(observed$base_test_1$frozen)
  expect_length(intersect(observed$base_test_1$samples, observed$base_tune_1$samples), 0)
  expect_equal(nrow(as.data.frame(result)), 4)
  expect_output(print(result), "completion_rate")
  expect_equal(jsonlite::fromJSON(result$summary$evaluation_record)$schema_version, "1.0.0")
  for (manifest in result$design$manifests) {
    expect_equal(file_sha256(manifest$manifest_path), manifest$manifest_sha256)
  }
})

test_that("test outcomes cannot rerank or change the frozen tuning decision", {
  f <- evaluation_fixture()
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth, "test_reversal"))
  result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    fraction = 0.25, repeats = 1L
  )
  expect_equal(result$selection$candidate, "good")
  expect_equal(result$selection$decision, "not_confirmed")
  expect_true(is.na(result$selection$recommended_candidate))
  frozen <- jsonlite::fromJSON(result$selection$frozen_path)
  expect_equal(frozen$tuning$candidate, "good")
})

test_that("missing predictions, missing DS, swaps and failures cannot win", {
  f <- evaluation_fixture()
  for (behavior in c("missing", "absent_ds", "swap", "drop", "failed")) {
    local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth, behavior))
    result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
      field = if (behavior == "absent_ds") "DS" else "GT", fraction = 0.25, repeats = 1L
    )
    expect_equal(result$selection$decision, "incomplete", info = behavior)
    expect_true(is.na(result$selection$recommended_candidate))
    expect_true(all(result$metrics$n_predicted + result$metrics$n_missing == result$metrics$n_eligible))
    expect_false(any(result$metrics$partition == "test"))
  }
})

test_that("sample permutations and DS evaluation give correct paired scores", {
  f <- evaluation_fixture()
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth, "permuted"))
  result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    field = "DS", fraction = 0.25, repeats = 1L
  )
  expect_equal(result$selection$decision, "confirmed")
  expect_equal(result$metrics$rmse[result$metrics$candidate == "good"], c(0, 0))
})

test_that("budget and output guards reject before any engine launch", {
  f <- evaluation_fixture()
  local_mocked_bindings(run_beagle = function(...) stop("Must not launch"))
  directory <- withr::local_tempdir()
  expect_error(evaluate_imputation(f$path, f$samples, f$candidates, "base", directory, max_runs = 1),
    class = "galick_error_evaluation_invalid_input"
  )
  expect_length(list.files(directory), 0)
  writeLines("existing", file.path(directory, "file"))
  expect_error(evaluate_imputation(f$path, f$samples, f$candidates, "base", directory),
    class = "galick_error_evaluation_invalid_input"
  )
})

test_that("reference sample leakage is rejected before execution", {
  f <- evaluation_fixture()
  f$truth$gt[, ] <- gsub("/", "|", f$truth$gt, fixed = TRUE)
  reference <- withr::local_tempfile(fileext = ".vcf")
  write_genotype_vcf(f$truth, reference)
  candidates <- lapply(f$candidates, function(x) c(x, list(ref = reference)))
  local_mocked_bindings(run_beagle = function(...) stop("Must not launch"))
  expect_error(evaluate_imputation(f$path, f$samples, candidates, "base", withr::local_tempdir()),
    class = "galick_error_evaluation_invalid_input"
  )
})

test_that("array evaluation removes hidden markers and leaves sources unchanged", {
  f <- evaluation_fixture()
  reference <- evaluation_reference(f)
  before <- file_sha256(c(f$path, reference))
  candidates <- lapply(f$candidates, function(x) c(x, list(ref = reference)))
  observed <- numeric()
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth, observer = function(input, seed, directory) {
    observed <<- c(observed, input$fixed$POS)
  }))
  panel <- data.frame(chrom = "1", pos = c(100, 300), ref = "A", alt = "G")
  result <- evaluate_imputation(f$path, f$samples, candidates, "base", withr::local_tempdir(),
    design = "array", observed_markers = panel, repeats = 1L, field = "DS"
  )
  expect_true(all(observed %in% panel$pos))
  expect_equal(result$selection$decision, "confirmed")
  expect_true(all(result$metrics$n_eligible == 40))
  expect_equal(file_sha256(c(f$path, reference)), before)
})

test_that("changing test truth cannot affect tuning masks or selection", {
  f <- evaluation_fixture()
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth))
  first <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    fraction = 0.25, repeats = 1L
  )
  changed <- f$truth
  changed$gt[, 5:8] <- "1/1"
  path <- withr::local_tempfile(fileext = ".vcf")
  write_genotype_vcf(changed, path)
  local_mocked_bindings(run_beagle = evaluation_mock_runner(changed))
  second <- evaluate_imputation(path, f$samples, f$candidates, "base", withr::local_tempdir(),
    fraction = 0.25, repeats = 1L
  )
  expect_equal(
    first$metrics[first$metrics$partition == "tune", ],
    second$metrics[second$metrics$partition == "tune", ]
  )
  expect_equal(first$selection$frozen_sha256, second$selection$frozen_sha256)
  expect_equal(first$design$manifests$tune_1$input_sha256, second$design$manifests$tune_1$input_sha256)
})

test_that("ties test only the baseline and do not invent a winner", {
  f <- evaluation_fixture()
  f$candidates$base$extra_args <- "iterations=2"
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth))
  result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    fraction = 0.25, repeats = 1L
  )
  expect_equal(result$selection$decision, "no_clear_winner")
  expect_true(is.na(result$selection$candidate))
  expect_length(result$runs, 3)
  expect_equal(result$metrics$candidate[result$metrics$partition == "test"], "base")
})

test_that("real Beagle completes paired sporadic evaluation", {
  skip_if_no_java()
  skip_if_no_beagle()
  f <- evaluation_fixture()
  candidates <- list(base = list(), alternative = list(extra_args = "iterations=4"))
  result <- evaluate_imputation(f$path, f$samples, candidates, "base", withr::local_tempdir(),
    fraction = 0.25, repeats = 1L, memory_gb = 1L
  )
  expect_true(result$ok)
  expect_true(all(result$metrics$status == "ok"))
  expect_true(all(result$metrics$completion_rate == 1))
  expect_lte(length(result$runs), 4)
  expect_true(all(vapply(result$runs, function(x) file.exists(x$run_record), logical(1))))
})

test_that("real Beagle completes reference-based DS array evaluation", {
  skip_if_no_java()
  skip_if_no_beagle()
  f <- evaluation_fixture()
  reference <- evaluation_reference(f)
  candidates <- list(base = list(ref = reference), alternative = list(ref = reference, extra_args = "iterations=4"))
  panel <- data.frame(chrom = "1", pos = seq(100, 1200, 200), ref = "A", alt = "G")
  result <- evaluate_imputation(f$path, f$samples, candidates, "base", withr::local_tempdir(),
    design = "array", observed_markers = panel, field = "DS", repeats = 1L, memory_gb = 1L
  )
  expect_true(result$ok)
  expect_true(all(result$metrics$status == "ok"))
  expect_true(all(result$metrics$completion_rate == 1))
  expect_true(all(result$metrics$n_eligible == 24))
  expect_true(all(vapply(result$runs, function(x) x$status == 0L, logical(1))))
})

test_that("three candidates respect exact job budget and candidate-order-independent masks", {
  f <- evaluation_fixture()
  f$candidates$other <- list(extra_args = "iterations=3")
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth))
  first <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    repeats = 1L, fraction = 0.25, max_runs = 5L
  )
  second <- evaluate_imputation(f$path, f$samples, rev(f$candidates), "base", withr::local_tempdir(),
    repeats = 1L, fraction = 0.25, max_runs = 5L
  )
  expect_length(first$runs, 5)
  expect_length(second$runs, 5)
  expect_equal(first$selection$candidate, second$selection$candidate)
  expect_equal(first$design$manifests$tune_1$manifest_sha256, second$design$manifests$tune_1$manifest_sha256)
})

test_that("unexpected job exceptions retain diagnostics and cannot select", {
  f <- evaluation_fixture()
  local_mocked_bindings(run_beagle = function(...) stop("Launch failed"))
  result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    repeats = 1L, fraction = 0.25
  )
  expect_false(result$ok)
  expect_equal(result$selection$decision, "incomplete")
  expect_equal(nrow(result$errors), 2L)
  expect_true(all(grepl("Launch failed", result$errors$message)))
  expect_true(file.exists(result$summary$evaluation_record))
})

test_that("changed final held-out manifest invalidates a recommendation", {
  f <- evaluation_fixture()
  observer <- function(input, seed, directory) {
    if (grepl("good_test_", directory, fixed = TRUE)) {
      root <- dirname(dirname(directory))
      cat("\n", file = file.path(root, "test_1_heldout.tsv"), append = TRUE)
    }
  }
  local_mocked_bindings(run_beagle = evaluation_mock_runner(f$truth, observer = observer))
  result <- evaluate_imputation(f$path, f$samples, f$candidates, "base", withr::local_tempdir(),
    repeats = 1L, fraction = 0.25
  )
  expect_false(result$ok)
  expect_equal(result$selection$decision, "incomplete")
  expect_true(is.na(result$selection$recommended_candidate))
  expect_true("evaluation_input_changed" %in% result$errors$code)
})
