test_that("panel design confirms an improvement and records fixed scoring evidence", {
  f <- panel_fixture()
  args <- panel_arguments(f)
  args$repeats <- 2L
  args$max_runs <- 16L
  observed <- list()
  local_mocked_bindings(run_beagle = panel_mock_runner(f$truth, observer = function(input, seed, dir) {
    sample.int(100, 1)
    observed[[basename(dir)]] <<- list(
      input = input, seed = seed,
      frozen = file.exists(file.path(dirname(dirname(dir)), "selection.json"))
    )
  }))
  before <- file_sha256(c(f$path, f$reference))
  set.seed(713)
  rng <- .Random.seed
  result <- do.call(design_marker_panel, args)
  expect_identical(.Random.seed, rng)
  expect_true(result$ok)
  expect_equal(result$selection$decision, "confirmed")
  expect_s3_class(result, "galick_panel_design")
  expect_s3_class(result, "galick_result")
  expect_equal(result$recommended_panel$pos, c(100, 100))
  expect_lte(length(result$runs), 16)
  expect_true(all(result$metrics$n_eligible == 320))
  expect_true(all(result$metrics$n_missing == 0))
  expect_equal(file_sha256(c(f$path, f$reference)), before)
  expect_output(print(result), "completion_rate")
  expect_equal(nrow(as.data.frame(result)), nrow(result$metrics))
  for (x in observed) {
    expect_true(all(x$input$fixed$POS <= 400))
    expect_true(all(x$input$formats == "GT"))
    expect_length(x$input$format_fields, 0)
  }
  expect_equal(observed$baseline_design_1$seed, observed$ld_seed_design_1$seed)
  expect_false(observed$baseline_tune_1$frozen)
  expect_true(observed$baseline_test_1$frozen)
  record <- jsonlite::fromJSON(file.path(args$output_dir, "panel_design.json"), simplifyVector = FALSE)
  expect_equal(record$schema_version, "1.0.0")
  expect_true(panel_integrity(record$artifacts))
  expect_equal(nrow(result$search_history), 4)
})

test_that("ties, missing predictions and failed confirmation never recommend a panel", {
  f <- panel_fixture()
  for (behavior in c("tie", "reversal", "missing", "missing_baseline", "absent_ds", "failed")) {
    local_mocked_bindings(run_beagle = panel_mock_runner(f$truth, behavior))
    result <- do.call(design_marker_panel, panel_arguments(f))
    expect_null(result$recommended_panel, info = behavior)
    expect_true("recommended_panel" %in% names(result))
    expected <- if (behavior == "reversal") {
      "not_confirmed"
    } else if (
               behavior %in% c("absent_ds", "failed", "missing_baseline")) {
      "incomplete"
    } else {
      "no_clear_winner"
    }
    expect_equal(result$selection$decision, expected, info = behavior)
  }
})

test_that("test truth cannot change design search or tuning selection", {
  f <- panel_fixture()
  local_mocked_bindings(run_beagle = panel_mock_runner(f$truth))
  first <- do.call(design_marker_panel, panel_arguments(f))
  changed <- f$truth
  cols <- f$samples$sample_id[f$samples$partition == "test"]
  changed$gt[, match(cols, changed$samples)] <- "0/0"
  write_genotype_vcf(changed, f$path, overwrite = TRUE)
  second <- do.call(design_marker_panel, panel_arguments(f))
  expect_equal(first$candidate_panel, second$candidate_panel)
  expect_equal(
    first$search_history[c("proxy", "rmse", "accepted")],
    second$search_history[c("proxy", "rmse", "accepted")]
  )
  expect_equal(first$selection$candidate, second$selection$candidate)
  changed$gt[, match(f$samples$sample_id[f$samples$partition == "tune"], changed$samples)] <- "0/0"
  write_genotype_vcf(changed, f$path, overwrite = TRUE)
  third <- do.call(design_marker_panel, panel_arguments(f))
  expect_equal(first$candidate_panel, third$candidate_panel)
  expect_equal(
    first$search_history[c("proxy", "rmse", "accepted")],
    third$search_history[c("proxy", "rmse", "accepted")]
  )
})

test_that("an overall improvement cannot hide a group regression", {
  f <- panel_fixture()
  local_mocked_bindings(run_beagle = panel_mock_runner(f$truth, "group_harm"))
  result <- do.call(design_marker_panel, panel_arguments(f))
  expect_true(result$ok)
  expect_equal(result$selection$decision, "no_clear_winner")
  expect_null(result$recommended_panel)
  seed <- result$search_history[result$search_history$panel == "ld_seed", ]
  expect_lt(seed$rmse, result$search_history$rmse[1])
  expect_false(seed$accepted)
})

test_that("original missing scoring calls have identical denominators across panels", {
  f <- panel_fixture()
  f$truth$gt[5, 41] <- "./."
  write_genotype_vcf(f$truth, f$path, overwrite = TRUE)
  args <- panel_arguments(f)
  args$max_proposals <- 0L
  args$field <- "GT"
  local_mocked_bindings(run_beagle = panel_mock_runner(f$truth))
  result <- do.call(design_marker_panel, args)
  expect_true(result$ok)
  expect_true(all(result$metrics$n_eligible[result$metrics$partition == "tune"] == 159))
  expect_true(all(result$metrics$n_eligible[result$metrics$partition == "design"] == 160))
})

test_that("input tampering stops execution and invalidates recommendations", {
  f <- panel_fixture()
  local_mocked_bindings(run_beagle = panel_mock_runner(f$truth, observer = function(input, seed, dir) {
    writeLines("changed", file.path(dirname(dirname(dir)), "experiment.json"))
  }))
  result <- do.call(design_marker_panel, panel_arguments(f))
  expect_false(result$ok)
  expect_equal(result$selection$decision, "incomplete")
  expect_null(result$recommended_panel)
  expect_true("panel_design_input_changed" %in% result$errors$code)
})

test_that("input contracts and budgets reject without launching Beagle", {
  f <- panel_fixture()
  args <- panel_arguments(f)
  local_mocked_bindings(run_beagle = function(...) stop("Must not launch"))
  changes <- list(
    list(max_runs = 1), list(n_markers = 20), list(max_genotype_cells = 1),
    list(max_ld_pairs = 1), list(min_spacing_bp = 1000, n_markers = 4), list(repeats = 0), list(seed = -1),
    list(extra_args = "gt=forbidden"), list(extra_args = "ne=100"),
    list(scoring_markers = f$candidates), list(mandatory_markers = f$scoring),
    list(mandatory_markers = f$candidates[1, ], excluded_markers = f$candidates[1, ]),
    list(candidate_markers = f$candidates[FALSE, ])
  )
  for (change in changes) {
    modified <- args
    modified[names(change)] <- change
    expect_error(do.call(design_marker_panel, modified),
      class = "galick_error_panel_design_invalid_input", info = paste(names(change), collapse = ",")
    )
  }
  expect_length(list.files(args$output_dir), 0)
  writeLines("existing", file.path(args$output_dir, "keep"))
  expect_error(do.call(design_marker_panel, args), class = "galick_error_panel_design_invalid_input")
})

test_that("sample groups, blocks and exact marker identities are checked", {
  f <- panel_fixture()
  expect_error(panel_marker_rows(f$candidates[c(1, 1), ], f$truth, "markers"),
    class = "galick_error_panel_design_invalid_input"
  )
  expect_error(panel_marker_rows(data.frame(), f$truth, "markers"), class = "galick_error_panel_design_invalid_input")
  bad <- f$samples
  bad$block <- "same"
  expect_error(panel_samples(bad, f$truth), class = "galick_error_panel_design_invalid_input")
  expect_error(panel_samples(f$samples[-1, ], f$truth), class = "galick_error_panel_design_invalid_input")
  bad <- f$samples
  bad$sample_id[1] <- "unknown"
  expect_error(panel_samples(bad, f$truth), class = "galick_error_panel_design_invalid_input")
  bad$sample_id[1] <- NA_character_
  expect_error(panel_samples(bad, f$truth), class = "galick_error_panel_design_invalid_input")
  expect_error(panel_samples(data.frame(), f$truth), class = "galick_error_panel_design_invalid_input")
  bad <- f$samples
  bad$group <- NULL
  expect_true(all(panel_samples(bad, f$truth)$group == "all"))
})

test_that("missing scoring truth and ineligible mandatory markers are rejected", {
  f <- panel_fixture()
  x <- f$truth
  x$gt[1, 1:40] <- "./."
  write_genotype_vcf(x, f$path, overwrite = TRUE)
  expect_error(panel_inputs(
    f$path, f$reference, f$samples, f$candidates, f$scoring,
    f$candidates[1, ], NULL, 1e6
  ), class = "galick_error_panel_design_invalid_input")
  x$gt[c(5:6, 11:12), 1:20] <- "./."
  write_genotype_vcf(x, f$path, overwrite = TRUE)
  expect_error(panel_inputs(
    f$path, f$reference, f$samples, f$candidates, f$scoring,
    NULL, NULL, 1e6
  ), class = "galick_error_panel_design_invalid_input")
})

test_that("real Beagle completes a bounded panel experiment", {
  skip_if_no_java()
  skip_if_no_beagle()
  f <- panel_fixture()
  args <- panel_arguments(f)
  args$max_proposals <- 0L
  args$n_markers <- 4L
  args$extra_args <- c("iterations=2", "burnin=1")
  args$memory_gb <- 1L
  result <- do.call(design_marker_panel, args)
  expect_true(result$ok, info = paste(result$errors$message, collapse = "; "))
  expect_true(all(result$metrics$completion_rate == 1))
  expect_lte(result$summary$n_runs, 6)
  expect_true(file.exists(file.path(args$output_dir, "panel_design.json")))
  for (i in seq_len(nrow(result$metrics))) {
    metric <- result$metrics[i, ]
    job <- result$runs[[paste(metric$candidate, metric$partition, 1, sep = "_")]]
    con <- gzfile(job$output_paths, "rt")
    lines <- readLines(con)
    close(con)
    header <- strsplit(lines[startsWith(lines, "#CHROM")], "\t")[[1]]
    body <- do.call(rbind, strsplit(lines[!startsWith(lines, "#")], "\t"))
    rows <- match(paste(f$scoring$chrom, f$scoring$pos), paste(body[, 1], body[, 2]))
    expected <- genotype_counts(f$truth)[c(5:6, 11:12), match(header[-(1:9)], f$truth$samples), drop = FALSE]
    prediction <- t(vapply(rows, function(row) {
      index <- match("DS", strsplit(body[row, 9], ":", fixed = TRUE)[[1]])
      vapply(strsplit(body[row, -(1:9)], ":", fixed = TRUE), function(x) as.numeric(x[index]), numeric(1))
    }, numeric(ncol(expected))))
    expect_equal(sqrt(mean((expected - prediction)^2)), metric$rmse)
  }
})

test_that("reference leakage, absent markers and non-SNP input cannot reach execution", {
  f <- panel_fixture()
  expect_error(panel_inputs(f$path, f$path, f$samples, f$candidates, f$scoring, NULL, NULL, 1e6),
    class = "galick_error_panel_design_invalid_input"
  )
  ref <- read_genotype_vcf(f$reference)
  ref$fixed <- ref$fixed[-6]
  ref$gt <- ref$gt[-6, , drop = FALSE]
  ref$formats <- ref$formats[-6]
  write_genotype_vcf(ref, f$reference, overwrite = TRUE)
  expect_error(panel_inputs(f$path, f$reference, f$samples, f$candidates, f$scoring, NULL, NULL, 1e6),
    class = "galick_error_panel_design_invalid_input"
  )
  f$truth$fixed$REF[1] <- "AA"
  write_genotype_vcf(f$truth, f$path, overwrite = TRUE)
  expect_error(panel_inputs(f$path, f$reference, f$samples, f$candidates, f$scoring, NULL, NULL, 1e6),
    class = "galick_error_panel_design_invalid_input"
  )
})
