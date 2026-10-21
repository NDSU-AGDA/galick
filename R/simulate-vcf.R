#' Simulate founder-conditioned breeding genotypes
#'
#' Generates linked descendants using Haldane meiosis (Poisson gametic crossovers,
#' no interference) and a supplied genetic map. This is not de novo demographic
#' simulation. Initial LD, allele ascertainment and population structure come from
#' the founders. No mutation, migration, selection or read simulation is performed.
#' No new runtime dependency is required. Plain VCF bundles are staged and published
#' only after validation; nonempty destinations are refused.
#'
#' @param founder_vcf Complete chromosome-consistent phased biallelic A/C/G/T SNP VCF.
#'   Multiple heterozygous phase blocks per chromosome/sample are refused.
#' @param map Data frame or tab-delimited file with chrom, pos, cm. Exact coverage,
#'   unique positive physical positions and nondecreasing nonnegative cM are required.
#' @param output_dir New or empty directory for the output bundle.
#' @param scenario Random mating, biparental F2, or independent selfed lineages.
#' @param inheritance Diploid or random-bivalent autotetraploid (outbred only).
#'   Tetraploids exclude quadrivalents, double reduction and preferential pairing.
#'   These files do not enable polyploid Beagle or evaluation support.
#' @param founder_data Optional data frame assigning every sample_id to a population.
#' @param parents Two homozygous founder sample IDs for f2/selfed; NULL for outbred.
#' @param n_target,n_reference Counts per population; named vectors for multiple
#'   outbred populations. Targets must be positive; reference counts may be zero.
#' @param population_size Outbred census size per population (not effective size).
#' @param generations Forward generations from founders (generation 0). Defaults:
#'   outbred 3, f2 2, selfed 6. F2 requires 2; F6 has five selfings after F1.
#' @param observed_markers Optional array data frame with chrom, pos, ref, alt.
#' @param missingness List with rate between 0 and 1 and sample_sd/marker_sd
#'   between 0 and 10 (inclusive).
#'   Gaussian logit effects are calibrated to the expected retained-cell missing rate.
#' @param error_rate Probability of flipping one allele copy in each retained called
#'   target genotype. No likelihood, depth or genotype-probability fields are invented.
#' @param seed Nonnegative integer seed; caller RNG state is restored. Observation
#'   settings use a separate seed and cannot change genetic truth.
#' @param max_allele_cells Conservative ceiling on simulation allele-array entries,
#'   checked after reading founders but before population-array allocation. This is
#'   not a bound on total R memory or output size; GT strings and ledgers add overhead.
#' @return A beagler_simulation inheriting from beagler_result, with paths, samples,
#'   observation rates, diagnostics and versioned model/input/output provenance.
#' @export
simulate_vcf <- function(founder_vcf, map, output_dir,
                         scenario = c("outbred", "f2", "selfed"),
                         inheritance = c("diploid", "autotetraploid_bivalent"),
                         founder_data = NULL, parents = NULL, n_target = 100L, n_reference = 100L,
                         population_size = 500L, generations = NULL, observed_markers = NULL,
                         missingness = list(rate = 0.05, sample_sd = 0, marker_sd = 0),
                         error_rate = 0.001, seed = 2026L, max_allele_cells = 2e7) {
  the_call <- match.call()
  scenario <- match.arg(scenario)
  if (length(inheritance) != 1L && !identical(inheritance, c("diploid", "autotetraploid_bivalent"))) {
    simulation_abort("Specify one supported inheritance model.", "beagler_error_simulation_unsupported_inheritance")
  }
  if (identical(inheritance, c("diploid", "autotetraploid_bivalent"))) inheritance <- "diploid"
  if (!is.character(inheritance) || length(inheritance) != 1L || is.na(inheritance) ||
        !inheritance %in% c("diploid", "autotetraploid_bivalent") ||
        (inheritance != "diploid" && scenario != "outbred")) {
    simulation_abort("Unsupported inheritance/scenario combination.",
      "beagler_error_simulation_unsupported_inheritance"
    )
  }
  ploidy <- if (inheritance == "diploid") 2L else 4L
  simulation_number(seed, "seed", 0, .Machine$integer.max, TRUE)
  simulation_number(population_size, "population_size", 2, .Machine$integer.max, TRUE)
  simulation_number(max_allele_cells, "max_allele_cells", 1)
  if (is.null(generations)) {
    generations <- switch(scenario,
      outbred = 3L,
      f2 = 2L,
      selfed = 6L
    )
  }
  simulation_number(generations, "generations", if (scenario == "outbred") 1 else 2, 1000, TRUE)
  if (scenario == "f2" && generations != 2) simulation_abort("F2 requires generations = 2.")
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    simulation_abort("output_dir must be a new or empty directory.")
  }
  simulation_output_available(output_dir)
  founders <- simulation_founders(founder_vcf, map, ploidy, founder_data)
  if (scenario == "outbred") {
    if (!is.null(parents)) simulation_abort("parents applies only to f2/selfed.")
    groups <- unique(founders$metadata$population)
    if (any(table(founders$metadata$population) < 2L)) {
      simulation_abort("Each outbred population needs at least two founders.")
    }
  } else {
    if (!is.character(parents) || length(parents) != 2L || anyNA(parents) ||
          anyDuplicated(parents) || any(!parents %in% founders$vcf$samples)) {
      simulation_abort("Supply two distinct homozygous founder IDs as parents.")
    }
    calls <- founders$vcf$gt[, match(parents, founders$vcf$samples), drop = FALSE]
    if (any(!calls %in% c("0|0", "1|1"))) simulation_abort("Cross parents must be homozygous at every retained locus.")
    groups <- "cross"
  }
  n_target <- simulation_counts(n_target, groups, "n_target", 1)
  n_reference <- simulation_counts(n_reference, groups, "n_reference", 0)
  totals <- n_target + n_reference
  if (any(totals > .Machine$integer.max)) simulation_abort("Requested cohort is too large.")
  if (scenario == "outbred" && any(totals > population_size)) {
    simulation_abort("Census must accommodate disjoint target/reference draws.")
  }
  census <- if (scenario == "outbred") population_size * length(groups) else sum(pmax(2, totals))
  estimated <- nrow(founders$vcf$gt) * ploidy * (ncol(founders$vcf$gt) + 2 * census + 3 * sum(totals))
  if (!is.finite(estimated) || estimated > max_allele_cells) {
    simulation_abort("Simulation exceeds max_allele_cells before population allocation.",
      "beagler_error_simulation_resource_limit"
    )
  }
  settings <- simulation_noise_settings(missingness, error_rate, observed_markers, founders)
  source_hash <- file_sha256(founder_vcf)
  map_file_hash <- if (is.character(map)) file_sha256(map) else NULL
  output_dir <- normalize_path_safe(output_dir, must_exist = FALSE)
  if (!dir.exists(dirname(output_dir))) dir.create(dirname(output_dir), recursive = TRUE)
  stage <- tempfile("beagler_simulation_", tmpdir = dirname(output_dir))
  if (!dir.create(stage)) simulation_abort("Cannot create staging directory.")
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  event_path <- file.path(stage, "recombination.tsv")
  event_count <- 0
  event_sink <- function(events) {
    data.table::fwrite(events, event_path, sep = "\t", append = event_count > 0)
    event_count <<- event_count + nrow(events)
  }
  withr::local_seed(as.integer(seed))
  population <- simulation_populations(
    founders, scenario, parents, n_target, n_reference,
    population_size, generations, event_sink
  )
  truth <- simulation_vcf(founders, population, "target")
  reference <- if (sum(n_reference)) simulation_vcf(founders, population, "reference") else NULL
  observation_seed <- as.integer((as.double(seed) + 104729) %% .Machine$integer.max)
  withr::local_seed(observation_seed)
  observed <- simulation_observe(truth, settings, ploidy)
  if (!identical(file_sha256(founder_vcf), source_hash) ||
        (!is.null(map_file_hash) && !identical(file_sha256(map), map_file_hash))) {
    simulation_abort("An input changed during simulation.")
  }
  new_beagler_simulation(
    founders, population, truth, reference, observed, stage, output_dir,
    list(
      model_version = "1.0.0", scenario = scenario, inheritance = inheritance, ploidy = ploidy,
      seed = seed, observation_seed = observation_seed, rng_kind = RNGkind(), n_target = n_target,
      n_reference = n_reference, population_size = population_size, generations = generations,
      parents = parents, missingness = settings$missingness, error_rate = error_rate,
      estimated_allele_cells = estimated, max_allele_cells = max_allele_cells,
      founder_path = normalize_path_safe(founder_vcf, TRUE), founder_sha256 = source_hash,
      map_file_sha256 = map_file_hash
    ), the_call
  )
}

#' @keywords internal
#' @noRd
simulation_output_available <- function(path) {
  if (file.exists(path) && (!dir.exists(path) || length(list.files(path, all.files = TRUE, no.. = TRUE)))) {
    simulation_abort("Refusing a nonempty output destination.")
  }
}
