test_that("evaluation plans are deterministic, paired and restore caller RNG", {
  f <- evaluation_fixture()
  set.seed(123)
  before <- .Random.seed
  plans <- evaluation_plan(f$truth, f$samples, "sporadic", NULL, 0.25, 2L, 42L)
  expect_identical(.Random.seed, before)
  again <- evaluation_plan(f$truth, f$samples, "sporadic", NULL, 0.25, 2L, 42L)
  expect_identical(plans, again)
  expect_false(identical(plans$tune_1$eligible, plans$tune_2$eligible))
  expect_equal(sum(plans$tune_1$eligible), 12L)
  expect_setequal(unique(plans$tune_1$samples$partition), "tune")
  expect_setequal(unique(plans$test_1$samples$partition), "test")
})

test_that("engine inputs strip truth-derived fields and phase", {
  f <- evaluation_fixture()
  x <- f$truth
  x$gt[, ] <- gsub("/", "|", x$gt, fixed = TRUE)
  x$fixed$INFO <- "AF=0.5;AC=8;AN=16"
  x$format_fields <- list(DS = matrix("1", 12, 8), PS = matrix("7", 12, 8))
  x$meta <- c(
    x$meta, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">",
    "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"DS\">"
  )
  plan <- evaluation_plan(x, f$samples, "sporadic", NULL, 0.25, 1L, 1L)[[1]]
  input <- evaluation_input(x, plan, "sporadic")
  expect_length(input$format_fields, 0L)
  expect_false(any(grepl("|", input$gt, fixed = TRUE)))
  expect_true(all(input$fixed$INFO == "."))
  expect_false(any(grepl("^##INFO=|^##FORMAT=", input$meta)))
  expect_true(all(input$gt[plan$eligible] == "./."))
  expect_identical(x$fixed$INFO, rep("AF=0.5;AC=8;AN=16", 12))
})

test_that("array design physically removes untyped records and retains heldout denominator", {
  f <- evaluation_fixture()
  panel <- data.frame(chrom = "1", pos = c(100, 300), ref = "A", alt = "G")
  plans <- evaluation_plan(f$truth, f$samples, "array", panel, 0.1, 2L, 1L)
  input <- evaluation_input(f$truth, plans$tune_1, "array")
  expect_equal(input$fixed$POS, c(100, 300))
  expect_equal(sum(plans$tune_1$eligible), 40)
  expect_identical(plans$tune_1$eligible, plans$tune_2$eligible)
  expect_false(identical(plans$tune_1$seed, plans$tune_2$seed))
  expect_equal(nrow(evaluation_manifest(f$truth, plans$tune_1)), 40)
})

test_that("invalid partitions, blocks, candidates and zero eligibility fail before execution", {
  f <- evaluation_fixture()
  expect_error(evaluation_samples(transform(f$samples, block = "one"), f$truth$samples),
    class = "galick_error_evaluation_invalid_input"
  )
  expect_error(evaluation_samples(transform(f$samples, partition = "tune"), f$truth$samples))
  expect_error(evaluation_samples(f$samples[c(1, 1, 5), ], f$truth$samples))
  expect_error(evaluation_candidates(list(a = list(), b = list(extra_args = "excludesamples=x")), "a"))
  expect_error(evaluation_candidates(list(a = list(), b = list(extra_args = "ne=100")), "a"))
  expect_error(evaluation_plan(f$truth, f$samples, "sporadic", NULL, 0.00001, 1L, 1L))
  f$truth$gt[, ] <- "./."
  expect_error(evaluation_plan(f$truth, f$samples, "sporadic", NULL, 0.25, 1L, 1L))
})

test_that("originally missing and partially called truth is never eligible", {
  f <- evaluation_fixture()
  f$truth$gt[1, ] <- "./."
  f$truth$gt[2, ] <- "0/."
  plan <- evaluation_plan(f$truth, f$samples, "sporadic", NULL, 0.5, 1L, 42L)[[1]]
  expect_false(any(plan$eligible[1:2, ]))
  expect_equal(sum(plan$eligible), 20)
})

test_that("model settings and array identities are validated before execution", {
  f <- evaluation_fixture()
  for (argument in c("iterations=bogus", "iterations=0", "em=maybe", "window=1", "err=Inf")) {
    expect_error(evaluation_candidates(list(a = list(), b = list(extra_args = argument)), "a"),
      class = "galick_error_evaluation_invalid_input"
    )
  }
  expect_error(evaluation_candidates(stats::setNames(rep(list(list()), 5), letters[1:5]), "a"))
  panel <- data.frame(chrom = "1", pos = 100, ref = "G", alt = "A")
  expect_error(evaluation_plan(f$truth, f$samples, "array", panel, 0.1, 1L, 1L))
  panel <- data.frame(chrom = "1", pos = 1.5, ref = "A", alt = "G")
  expect_error(evaluation_plan(f$truth, f$samples, "array", panel, 0.1, 1L, 1L))
})
