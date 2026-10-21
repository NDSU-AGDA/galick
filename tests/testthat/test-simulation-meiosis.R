test_that("Haldane meiosis obeys analytic recombination probabilities", {
  withr::local_seed(519L)
  map <- c(0, 0.01, 0.1, 0.5, 1)
  hap <- cbind(rep(0L, 5), rep(1L, 5))
  simulated <- replicate(10000, simulation_meiosis(hap, map)$haplotypes[, 1])
  empirical <- rowMeans(simulated != rep(simulated[1, ], each = 5))
  expected <- (1 - exp(-2 * map)) / 2
  expect_equal(empirical[1], 0)
  expect_true(all(abs(empirical - expected) < 0.025))
  expect_true(abs(mean(simulated[1, ]) - 0.5) < 0.025)
})

test_that("zero-distance markers co-segregate and events replay parental inheritance", {
  withr::local_seed(3L)
  hap <- cbind(0:1, 1:0)
  for (i in 1:50) {
    result <- simulation_meiosis(hap, c(0, 0))
    expect_equal(sum(result$haplotypes), 1)
    expect_equal(nrow(result$events), 1)
    expect_true(result$events$parent_homolog %in% 1:2)
  }
})

test_that("random-bivalent tetraploid gametes match duplex segregation", {
  withr::local_seed(917L)
  hap <- matrix(c(0L, 0L, 1L, 1L), nrow = 1)
  dosage <- replicate(12000, sum(simulation_meiosis(hap, 0)$haplotypes))
  frequency <- tabulate(dosage + 1L, 3) / length(dosage)
  expect_true(all(abs(frequency - c(1, 4, 1) / 6) < 0.02))
  result <- simulation_meiosis(matrix(rep(1:4 %% 2, each = 3), 3), c(0, 0.1, 0.2))
  starts <- result$events[result$events$event == "start", ]
  expect_length(unique(c(starts$pair_1, starts$pair_2)), 4)
})
