test_that("evaluation metrics match independent arithmetic and undefined correlations", {
  cells <- tibble::tibble(
    marker = c(1, 1, 1, 2), group = "A", truth = c(0, 1, 2, 1),
    prediction = c(0, 2, NA, 1), gt = c(0, 2, NA, 1)
  )
  score <- evaluation_metrics(cells)
  expect_equal(score$n_eligible, 4)
  expect_equal(score$n_predicted + score$n_missing, score$n_eligible)
  expect_equal(score$rmse, sqrt(1 / 3))
  expect_equal(score$concordance_rate, 0.5)
  expect_equal(score$completion_rate, 0.75)
  expect_equal(score$dosage_r2, 1)
  expect_equal(score$n_undefined_r2, 1)
})

test_that("selection requires paired completeness, practical margins and group safety", {
  metrics <- data.frame(
    candidate = c("base", "good", "other"), status = "ok",
    completion_rate = 1, rmse = c(0.5, 0.2, 0.4)
  )
  groups <- data.frame(
    candidate = rep(metrics$candidate, each = 2), group = rep(c("A", "B"), 3),
    rmse = rep(metrics$rmse, each = 2)
  )
  choose <- function(m = metrics, g = groups) evaluation_select(m, g, metrics$candidate, "base", 0.01, 0.01)
  expect_equal(choose()$candidate, "good")
  expect_equal(choose(transform(metrics, completion_rate = c(1, 0.99, 1)))$decision, "incomplete")
  expect_equal(choose(transform(metrics, rmse = c(0.5, 0.2, 0.2)))$decision, "no_clear_winner")
  groups$rmse[groups$candidate == "good" & groups$group == "B"] <- 0.7
  expect_equal(choose(g = groups)$decision, "no_clear_winner")
})

test_that("RMSE is pooled by cell counts rather than averaging group RMSE", {
  cells <- tibble::tibble(
    marker = 1, group = c("A", rep("B", 9)), truth = 0,
    prediction = c(2, rep(0, 9)), gt = c(2, rep(0, 9))
  )
  expect_equal(evaluation_metrics(cells)$rmse, sqrt(4 / 10))
  expect_false(identical(evaluation_metrics(cells)$rmse, 1))
})
