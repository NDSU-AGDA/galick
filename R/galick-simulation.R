#' @keywords internal
#' @noRd
new_galick_simulation <- function(founders, population, truth, reference, observed, stage, destination,
                                  parameters, call) {
  names <- c(
    truth = "truth.vcf", observed = "observed.vcf", reference = "reference.vcf",
    samples = "samples.tsv", pedigree = "pedigree.tsv", recombination = "recombination.tsv",
    map = "genetic.map", marker_diagnostics = "marker_diagnostics.tsv", ld = "ld.tsv",
    observation_events = "observation_events.tsv", sample_missingness = "sample_missingness.tsv",
    marker_missingness = "marker_missingness.tsv", simulation_record = "simulation.json"
  )
  if (is.null(reference)) names <- names[names != "reference.vcf"]
  write_genotype_vcf(truth, file.path(stage, "truth.vcf"))
  write_genotype_vcf(observed$vcf, file.path(stage, "observed.vcf"))
  if (!is.null(reference)) write_genotype_vcf(reference, file.path(stage, "reference.vcf"))
  for (name in c("truth", "observed", if (!is.null(reference)) "reference")) {
    path <- file.path(stage, names[[name]])
    if (!validate_vcf(path)$ok) simulation_abort("Generated VCF failed validation; bundle was not published.")
    assert_resolvable_genotypes(read_genotype_vcf(path))
  }
  data.table::fwrite(population$samples, file.path(stage, "samples.tsv"), sep = "\t")
  data.table::fwrite(population$pedigree, file.path(stage, "pedigree.tsv"), sep = "\t", na = "NA")
  map <- data.frame(
    chrom = founders$map$chrom, marker_id = founders$vcf$fixed$ID,
    cm = founders$map$cm, pos = founders$map$pos
  )
  data.table::fwrite(map, file.path(stage, "genetic.map"), sep = "\t", col.names = FALSE)
  data.table::fwrite(observed$events, file.path(stage, "observation_events.tsv"), sep = "\t")
  data.table::fwrite(observed$by_sample, file.path(stage, "sample_missingness.tsv"), sep = "\t")
  data.table::fwrite(observed$by_marker, file.path(stage, "marker_missingness.tsv"), sep = "\t")
  samples <- as.data.frame(population$samples)
  target_metadata <- samples[samples$role == "target", , drop = FALSE]
  diagnostics <- list(
    simulation_diagnostics(founders$vcf, founders$metadata, "founder", parameters$ploidy)$by_marker,
    simulation_diagnostics(
      truth, target_metadata,
      "target", parameters$ploidy
    )$by_marker
  )
  ld <- list(
    simulation_ld(founders$vcf, founders$metadata, founders$map, "founder"),
    simulation_ld(truth, target_metadata, founders$map, "target")
  )
  if (!is.null(reference)) {
    ref_metadata <- samples[samples$role == "reference", , drop = FALSE]
    diagnostics[[3]] <- simulation_diagnostics(reference, ref_metadata, "reference", parameters$ploidy)$by_marker
    ld[[3]] <- simulation_ld(reference, ref_metadata, founders$map, "reference")
  }
  diagnostics <- data.table::rbindlist(diagnostics)
  ld <- data.table::rbindlist(ld)
  data.table::fwrite(diagnostics, file.path(stage, "marker_diagnostics.tsv"), sep = "\t", na = "NA")
  data.table::fwrite(ld, file.path(stage, "ld.tsv"), sep = "\t", na = "NA")
  artifacts <- lapply(names[names != "simulation.json"], function(name) {
    list(path = name, sha256 = file_sha256(file.path(stage, name)))
  })
  parameters$map <- founders$map
  parameters$founder_data <- founders$metadata
  parameters$map_sha256 <- digest::digest(jsonlite::toJSON(founders$map, digits = NA),
    algo = "sha256", serialize = FALSE
  )
  parameters$founder_data_sha256 <- digest::digest(jsonlite::toJSON(founders$metadata, digits = NA),
    algo = "sha256", serialize = FALSE
  )
  parameters$observed_markers <- observed$vcf$fixed[, c("CHROM", "POS", "REF", "ALT"), with = FALSE]
  parameters$ld_pair_selection <- "Up to ten adjacent marker lags per chromosome; dosage correlation within population."
  parameters$transmission_validation <- "Every gamete replayed against its recorded parental homolog switches."
  record_parameters <- parameters
  record_parameters$n_target <- as.list(parameters$n_target)
  record_parameters$n_reference <- as.list(parameters$n_reference)
  evaluation_write_json(
    list(
      schema_version = "1.0.0", parameters = record_parameters,
      observation = observed$summary, artifacts = artifacts,
      package_version = as.character(utils::packageVersion("galick")), r_version = as.character(getRversion())
    ),
    file.path(stage, "simulation.json")
  )
  simulation_output_available(destination)
  existed <- dir.exists(destination)
  if (existed) unlink(destination, recursive = TRUE)
  if (!file.rename(stage, destination)) {
    if (existed) dir.create(destination)
    simulation_abort("Could not publish the completed simulation bundle.")
  }
  paths <- stats::setNames(file.path(destination, names), names(names))
  object <- new_galick_result(summary = list(
    scenario = parameters$scenario, inheritance = parameters$inheritance,
    n_target = sum(parameters$n_target), n_reference = sum(parameters$n_reference), n_markers = nrow(truth$gt),
    generations = parameters$generations, output_dir = destination
  ), call = call)
  object$paths <- paths
  object$samples <- tibble::as_tibble(population$samples)
  object$observation <- observed$summary
  object$diagnostics <- tibble::as_tibble(diagnostics)
  object$ld <- tibble::as_tibble(ld)
  object$model <- parameters
  class(object) <- c("galick_simulation", class(object))
  object
}

#' @export
print.galick_simulation <- function(x, ...) {
  cli::cli_h1("Founder-conditioned simulation: {x$summary$scenario}")
  cli::cli_text("{x$summary$n_target} targets; {x$summary$n_reference} references; {x$summary$n_markers} SNPs")
  cli::cli_text("Inheritance: {x$summary$inheritance}; generation: {x$summary$generations}")
  cli::cli_text("Outputs: {.file {x$summary$output_dir}}")
  invisible(x)
}

#' @export
as.data.frame.galick_simulation <- function(x, ...) {
  as.data.frame(x$diagnostics)
}
