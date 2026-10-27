# Development-only comparison; AlphaSimR is not a galick dependency.
pkgload::load_all(".")
library(AlphaSimR)
n <- 10000L
map_m <- seq(0, 1, length.out = 11)
settings <- list(v = 1, p = 0, quadProb = 0, nThreads = 1)
make_population <- function(haplotypes, ploidy) {
  founder <- AlphaSimR::newMapPop(list(map_m), list(haplotypes), ploidy = ploidy)
  parameters <- AlphaSimR::SimParam$new(founder)
  parameters$v <- settings$v
  parameters$p <- settings$p
  parameters$quadProb <- settings$quadProb
  parameters$nThreads <- 1L
  list(pop = AlphaSimR::newPop(founder, simParam = parameters), parameters = parameters)
}
set.seed(5081)
diploid <- make_population(rbind(rep(0, 11), rep(1, 11), rep(0, 11), rep(1, 11)), 2L)
offspring <- AlphaSimR::makeCross(diploid$pop, matrix(c(1, 2), nrow = 1),
  nProgeny = n,
  simParam = diploid$parameters
)
hap <- AlphaSimR::pullSegSiteHaplo(offspring, simParam = diploid$parameters)[seq(1, 2 * n, 2), , drop = FALSE]
alpha_recombination <- colMeans(hap != hap[, 1])
set.seed(6091)
own_hap <- replicate(n, galick:::simulation_meiosis(cbind(rep(0L, 11), rep(1L, 11)), map_m)$haplotypes[, 1])
own_recombination <- rowMeans(own_hap != rep(own_hap[1, ], each = 11))
expected_recombination <- (1 - exp(-2 * map_m)) / 2
stopifnot(
  max(abs(alpha_recombination - expected_recombination)) < 0.025,
  max(abs(own_recombination - expected_recombination)) < 0.025,
  max(abs(own_recombination - alpha_recombination)) < 0.035
)

set.seed(7011)
tetraploid <- make_population(rbind(
  rep(0, 11), rep(0, 11), rep(1, 11), rep(1, 11),
  rep(0, 11), rep(0, 11), rep(1, 11), rep(1, 11)
), 4L)
offspring <- AlphaSimR::makeCross(tetraploid$pop, matrix(c(1, 2), nrow = 1),
  nProgeny = n,
  simParam = tetraploid$parameters
)
hap <- AlphaSimR::pullSegSiteHaplo(offspring, simParam = tetraploid$parameters)
alpha_gamete <- hap[seq(1, 4 * n, 4), 1] + hap[seq(2, 4 * n, 4), 1]
alpha_tetraploid <- tabulate(alpha_gamete + 1L, 3) / n
set.seed(8011)
own_tetraploid <- tabulate(replicate(n, sum(galick:::simulation_meiosis(
  matrix(c(0L, 0L, 1L, 1L), nrow = 1), 0
)$haplotypes)) + 1L, 3) / n
stopifnot(
  max(abs(alpha_tetraploid - c(1, 4, 1) / 6)) < 0.03,
  max(abs(own_tetraploid - c(1, 4, 1) / 6)) < 0.03,
  max(abs(alpha_tetraploid - own_tetraploid)) < 0.04
)

set.seed(9011)
parents <- make_population(rbind(rep(0, 11), rep(0, 11), rep(1, 11), rep(1, 11)), 2L)
f1 <- AlphaSimR::makeCross(parents$pop, matrix(c(1, 2), nrow = 1), simParam = parents$parameters)
f2 <- AlphaSimR::self(f1, nProgeny = n, simParam = parents$parameters)
dosage <- AlphaSimR::pullSegSiteGeno(f2, simParam = parents$parameters)
alpha_f2 <- tabulate(dosage[, 1] + 1L, 3) / n
set.seed(10011)
own_f2 <- tabulate(replicate(n, {
  h <- matrix(c(0L, 1L), nrow = 1)
  sum(galick:::simulation_meiosis(h, 0)$haplotypes) + sum(galick:::simulation_meiosis(h, 0)$haplotypes)
}) + 1L, 3) / n
stopifnot(
  max(abs(alpha_f2 - c(0.25, 0.5, 0.25))) < 0.03,
  max(abs(own_f2 - c(0.25, 0.5, 0.25))) < 0.03,
  max(abs(alpha_f2 - own_f2)) < 0.04
)
set.seed(11011)
lines <- AlphaSimR::makeCross(parents$pop, matrix(c(1, 2), nrow = 1),
  nProgeny = 1500,
  simParam = parents$parameters
)
alpha_heterozygosity <- 1
for (i in 1:5) {
  lines <- AlphaSimR::self(lines, simParam = parents$parameters)
  dosage <- AlphaSimR::pullSegSiteGeno(lines, simParam = parents$parameters)
  alpha_heterozygosity <- c(alpha_heterozygosity, mean(dosage == 1))
}
set.seed(12011)
own_lines <- matrix(rep(c(0L, 1L), 1500), ncol = 2, byrow = TRUE)
own_heterozygosity <- 1
for (i in 1:5) {
  for (j in seq_len(nrow(own_lines))) {
    h <- matrix(own_lines[j, ], nrow = 1)
    own_lines[j, ] <- c(
      galick:::simulation_meiosis(h, 0)$haplotypes,
      galick:::simulation_meiosis(h, 0)$haplotypes
    )
  }
  own_heterozygosity <- c(own_heterozygosity, mean(own_lines[, 1] != own_lines[, 2]))
}
stopifnot(
  max(abs(alpha_heterozygosity - 2^-(0:5))) < 0.04,
  max(abs(own_heterozygosity - 2^-(0:5))) < 0.04,
  max(abs(own_heterozygosity - alpha_heterozygosity)) < 0.04
)
output <- list(
  schema_version = "1.0.0", backend = "AlphaSimR", version = as.character(packageVersion("AlphaSimR")),
  r_version = as.character(getRversion()), settings = settings, n_gametes = n, n_selfed_lines = 1500L,
  implementation_sha256 = galick:::file_sha256("R/simulation-meiosis.R"),
  validation_script_sha256 = galick:::file_sha256("tools/validate-simulation-alphasimr.R"),
  seeds = c(5081, 6091, 7011, 8011, 9011, 10011, 11011, 12011), map_m = map_m,
  alpha_recombination = alpha_recombination, own_recombination = own_recombination,
  alpha_tetraploid_gametes = alpha_tetraploid, own_tetraploid_gametes = own_tetraploid,
  alpha_f2 = alpha_f2, own_f2 = own_f2, alpha_heterozygosity = alpha_heterozygosity,
  own_heterozygosity = own_heterozygosity,
  tolerances = list(recombination = 0.025, between_engines = 0.04, segregation = 0.03, selfing = 0.04)
)
args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args)) args[1] else "tests/testthat/fixtures/simulation-alphasimr.json"
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite::toJSON(output, pretty = TRUE, auto_unbox = TRUE, digits = NA), path)
cat("AlphaSimR validation passed; results:", path, "\n")
