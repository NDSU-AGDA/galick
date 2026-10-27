test_that("diploid simulations publish valid bundles without changing sources or RNG", {
  f <- simulation_fixture()
  hash <- file_sha256(f$path)
  withr::local_seed(456)
  before <- .Random.seed
  result <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    n_target = 12L,
    n_reference = 8L, population_size = 30L, generations = 2L, seed = 42L
  )
  expect_identical(.Random.seed, before)
  expect_s3_class(result, "galick_simulation")
  expect_s3_class(result, "galick_result")
  expect_true(result$ok)
  expect_equal(file_sha256(f$path), hash)
  expect_true(all(file.exists(result$paths)))
  expect_true(validate_vcf(result$paths[["truth"]], target = "beagle")$ok)
  expect_true(validate_vcf(result$paths[["observed"]], target = "beagle")$ok)
  expect_true(validate_reference_match(result$paths[["truth"]], result$paths[["reference"]], require_phased = TRUE)$ok)
  truth <- read_genotype_vcf(result$paths[["truth"]])
  observed <- read_genotype_vcf(result$paths[["observed"]])
  reference <- read_genotype_vcf(result$paths[["reference"]])
  expect_length(intersect(truth$samples, reference$samples), 0)
  expect_false(any(grepl("|", observed$gt, fixed = TRUE)))
  expect_equal(ncol(truth$gt), 12)
  expect_equal(ncol(reference$gt), 8)
  expect_equal(nrow(as.data.frame(result)), 36)
  record <- jsonlite::fromJSON(result$paths[["simulation_record"]], simplifyVector = FALSE)
  expect_equal(record$schema_version, "1.0.0")
  for (artifact in record$artifacts) {
    expect_equal(file_sha256(file.path(result$summary$output_dir, artifact$path)), artifact$sha256)
  }
})

test_that("observations cannot change complete genetic truth", {
  f <- simulation_fixture()
  first <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    n_target = 8, n_reference = 4,
    population_size = 15, generations = 1, missingness = list(rate = 0), error_rate = 0, seed = 7
  )
  second <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    n_target = 8, n_reference = 4,
    population_size = 15, generations = 1, missingness = list(rate = 0.5, sample_sd = 1, marker_sd = 1),
    error_rate = 1, seed = 7, observed_markers = data.frame(chrom = "1", pos = c(100, 300), ref = "A", alt = "G")
  )
  expect_equal(file_sha256(first$paths[["truth"]]), file_sha256(second$paths[["truth"]]))
  expect_equal(file_sha256(first$paths[["reference"]]), file_sha256(second$paths[["reference"]]))
  expect_equal(nrow(read_genotype_vcf(second$paths[["observed"]])$gt), 2)
  expect_equal(second$observation$expected_missing_rate, 0.5, tolerance = 1e-8)
  expect_equal(second$observation$n_errors, second$observation$n_called)
})

test_that("F2 and selfed lines have the declared pedigree and segregation", {
  f <- simulation_fixture(biparental = TRUE)
  f2 <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    scenario = "f2", parents = c("f1", "f2"),
    n_target = 600, n_reference = 0, missingness = list(rate = 0), error_rate = 0
  )
  dosage <- genotype_counts(read_genotype_vcf(f2$paths[["truth"]]))
  expect_true(all(abs(tabulate(dosage[1, ] + 1L, 3) / 600 - c(0.25, 0.5, 0.25)) < 0.06))
  expect_false("reference" %in% names(f2$paths))
  selfed <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    scenario = "selfed", parents = c("f1", "f2"),
    n_target = 200, n_reference = 0, generations = 5, missingness = list(rate = 0), error_rate = 0
  )
  dosage <- genotype_counts(read_genotype_vcf(selfed$paths[["truth"]]))
  expect_lt(abs(mean(dosage == 1) - 1 / 16), 0.03)
  pedigree <- data.table::fread(selfed$paths[["pedigree"]])
  expect_true(all(pedigree$parent_1[pedigree$generation > 1] == pedigree$parent_2[pedigree$generation > 1]))
})

test_that("tetraploid bundles retain four copies and remain outside Beagle", {
  f <- simulation_fixture(4L)
  result <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    inheritance = "autotetraploid_bivalent",
    n_target = 10, n_reference = 5, population_size = 20, generations = 1,
    missingness = list(rate = 1), error_rate = 0
  )
  expect_equal(detect_genotype_ploidy(result$paths[["truth"]]), 4L)
  expect_true(validate_vcf(result$paths[["truth"]])$ok)
  expect_false(validate_vcf(result$paths[["truth"]], target = "beagle")$ok)
  expect_true(all(read_genotype_vcf(result$paths[["observed"]])$gt == "./././."))
  expect_error(run_beagle(result$paths[["truth"]], tempfile()), class = "galick_error_run_beagle_invalid_input")
})

test_that("population labels separate mating pools", {
  f <- simulation_fixture()
  metadata <- data.frame(sample_id = paste0("f", 1:4), population = rep(c("A", "B"), each = 2))
  result <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    founder_data = metadata,
    n_target = c(A = 3, B = 4), n_reference = c(A = 2, B = 1), population_size = 8, generations = 2
  )
  pedigree <- data.table::fread(result$paths[["pedigree"]])
  parents <- match(pedigree$parent_1, pedigree$sample_id)
  created <- pedigree$generation > 0
  expect_equal(pedigree$population[parents[created]], pedigree$population[created])
  expect_equal(sum(result$samples$role == "target"), 7)
})

test_that("invalid founders, maps, inheritance and resource plans publish nothing", {
  f <- simulation_fixture()
  destination <- file.path(withr::local_tempdir(), "bundle")
  expect_error(simulate_vcf(f$path, f$map, destination, max_allele_cells = 1),
    class = "galick_error_simulation_resource_limit"
  )
  expect_false(dir.exists(destination))
  expect_error(simulate_vcf(f$path, f$map[-1, ], destination))
  expect_error(simulate_vcf(f$path, f$map, destination, inheritance = "allopolyploid"),
    class = "galick_error_simulation_unsupported_inheritance"
  )
  expect_error(simulate_vcf(f$path, f$map, destination, scenario = "f2", parents = c("f1", "f2")))
  map <- f$map
  map$cm[2] <- -1
  expect_error(simulate_vcf(f$path, map, destination))
  expect_false(dir.exists(destination))
  dir.create(destination)
  writeLines("keep", file.path(destination, "existing"))
  expect_error(simulate_vcf(f$path, f$map, destination))
  expect_equal(readLines(file.path(destination, "existing")), "keep")
})

test_that("unresolved phase and incorrect founder ploidy are refused", {
  f <- simulation_fixture()
  x <- read_genotype_vcf(f$path)
  x$format_fields$PS <- matrix(rep(c("1", "2"), 24), 12, 4)
  x$formats <- rep("GT:PS", nrow(x$gt))
  x$meta <- c(x$meta, '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">')
  path <- withr::local_tempfile(fileext = ".vcf")
  write_genotype_vcf(x, path)
  expect_error(simulate_vcf(path, f$map, tempfile()), "phase blocks")
  expect_error(simulate_vcf(f$path, f$map, tempfile(), inheritance = "autotetraploid_bivalent"), "declared ploidy")
})

test_that("independent parsers accept synthetic phase and tetraploid slots", {
  skip_if_no_bcftools()
  f <- simulation_fixture(4L)
  result <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    inheritance = "autotetraploid_bivalent",
    n_target = 4, n_reference = 2, population_size = 8, generations = 1
  )
  query <- run_cli(Sys.which("bcftools"), c("query", "-f", "[%GT\\n]", result$paths[["truth"]]))
  expect_equal(query$status, 0L)
  expect_match(query$stdout, "[01]\\|[01]\\|[01]\\|[01]")
  expect_false(grepl("not defined", query$stderr, fixed = TRUE))
})

test_that("diploid simulated truth completes evaluation without leakage", {
  skip_if_no_java()
  skip_if_no_beagle()
  f <- simulation_fixture()
  result <- simulate_vcf(f$path, f$map, withr::local_tempdir(),
    n_target = 8, n_reference = 8,
    population_size = 20, generations = 2
  )
  samples <- as.data.frame(result$samples[result$samples$role == "target", ])
  samples$partition <- rep(c("tune", "test"), each = 4)
  study <- evaluate_imputation(result$paths[["truth"]], samples,
    list(
      base = list(ref = result$paths[["reference"]]),
      alternative = list(ref = result$paths[["reference"]], extra_args = "iterations=4")
    ),
    baseline = "base", output_dir = withr::local_tempdir(), repeats = 1, fraction = 0.25, memory_gb = 1L
  )
  expect_true(study$ok)
  expect_true(all(study$metrics$status == "ok"))
})
