test_that("packing oracle preserves anchors and detects infeasible spacing", {
  expect_equal(panel_capacity(1:5, integer(), seq(0, 8, 2), 3), 3)
  expect_equal(panel_capacity(1:5, 2L, seq(0, 8, 2), 3), 2)
  expect_equal(panel_capacity(1:5, c(2L, 3L), seq(0, 8, 2), 3), -Inf)
  for (anchors in list(integer(), 1L, 2L, c(1L, 5L))) {
    sets <- unlist(lapply(0:5, function(k) combn(1:5, k, simplify = FALSE)), recursive = FALSE)
    feasible <- vapply(sets, function(s) {
      all(anchors %in% s) &&
        (length(s) < 2 || all(diff(seq(0, 8, 2)[s]) >= 3))
    }, logical(1))
    expect_equal(panel_capacity(1:5, anchors, seq(0, 8, 2), 3), max(lengths(sets)[feasible]))
  }
})

test_that("LD seed, quotas and swaps preserve every hard constraint", {
  f <- panel_fixture()
  input <- panel_inputs(f$path, f$reference, f$samples, f$candidates, f$scoring, NULL, NULL, 1e6)
  constraints <- panel_constraints(input, 2, 0)
  ld <- panel_ld(input, 1000, 100)
  baseline <- panel_seed(input, constraints)
  informed <- panel_seed(input, constraints, ld)
  expect_equal(input$truth$fixed$POS[baseline], c(200, 200))
  expect_equal(input$truth$fixed$POS[informed], c(100, 100))
  expect_gt(panel_proxy(informed, ld), panel_proxy(baseline, ld))
  expect_true(panel_feasible(informed, constraints, TRUE))
  expect_false(panel_feasible(input$candidates, constraints))
  expect_false(panel_feasible(integer(), constraints, TRUE))
  proposal <- panel_proposal(input, constraints, ld, baseline, paste(baseline, collapse = ","))
  expect_true(panel_feasible(proposal, constraints, TRUE))
  expect_length(setdiff(proposal, baseline), 1)
  zero <- panel_ld(input, 0, 100)
  expect_equal(zero$n_pairs, 0)
  expect_equal(panel_proxy(informed, zero), 0)
  expect_error(panel_ld(input, 1000, 1), class = "galick_error_panel_design_invalid_input")
  input$counts[, ] <- 0
  expect_equal(panel_proxy(informed, panel_ld(input, 1000, 100)), 0)
  input$counts[, ] <- NA_real_
  expect_equal(panel_proxy(informed, panel_ld(input, 1000, 100)), 0)
})

test_that("mandatory, excluded and distance rules apply to both seed methods", {
  f <- panel_fixture()
  input <- panel_inputs(
    f$path, f$reference, f$samples, f$candidates, f$scoring,
    f$candidates[c(1, 5), ], f$candidates[2, ], 1e6
  )
  constraints <- panel_constraints(input, 4, 150)
  ld <- panel_ld(input, 1000, 100)
  for (seed in list(panel_seed(input, constraints), panel_seed(input, constraints, ld))) {
    expect_length(seed, 4)
    expect_true(all(input$mandatory %in% seed))
    expect_true(panel_feasible(seed, constraints, TRUE))
  }
  input$mandatory <- input$candidates
  expect_error(panel_constraints(input, 2, 0), class = "galick_error_panel_design_invalid_input")
  constraints <- panel_constraints(input, length(input$candidates), 0)
  expect_null(panel_proposal(input, constraints, ld, input$candidates, character()))
})
