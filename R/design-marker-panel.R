#' Design and independently confirm a fixed-size imputation marker panel
#'
#' Constructs an evenly spaced baseline and a local LD-informed panel using
#' design samples only. A bounded single-swap search measures actual Beagle
#' recovery. Tuning selects one frozen challenger; untouched test samples confirm
#' it against the baseline. No global optimality or assay feasibility is claimed.
#'
#' @param truth_vcf Diploid biallelic A/C/G/T SNP VCF with trusted called genotypes.
#' @param reference_vcf Complete phased reference VCF with no truth sample IDs.
#' @param sample_data Data frame with sample_id and partition (design/tune/test),
#'   optional group and block. Each group needs 20 design and 10 tune/test samples.
#'   Declared blocks cannot cross partitions. Relationships are not inferred.
#' @param candidate_markers,scoring_markers Disjoint data frames with chrom, pos,
#'   ref, alt matching truth and reference exactly. Scoring markers never enter
#'   any panel; missing original calls are not scored. Design-only eligibility
#'   requires candidate call rate at least 95 percent and pooled MAF at least 1 percent.
#' @param n_markers Exact marker count, including obligatory markers.
#' @param output_dir New or empty directory. No automatic resume or overwrite.
#' @param mandatory_markers,excluded_markers Optional exact marker identity tables.
#'   Mandatory markers must be eligible candidates and cannot be excluded.
#' @param min_spacing_bp Minimum within-chromosome separation, including mandatory sites.
#' @param ld_window_bp Maximum physical distance for design-only candidate/scoring LD.
#'   Pairwise squared genotype correlation requires 20 jointly called samples
#'   and nonzero variance. Undefined pairs contribute zero proxy coverage.
#' @param map Optional Beagle genetic-map path. Without it Beagle assumes 1 cM/Mb.
#' @param extra_args Frozen Beagle model settings, with the same restrictions as
#'   [evaluate_imputation()]. No parameter tuning is performed.
#' @param field Score expected DS or GT allele counts; DS never falls back to GT.
#' @param repeats Number of paired engine-seed repeats, not biological replicates.
#' @param max_proposals Maximum measured swap proposals after the two seed panels.
#' @param max_runs Job budget; reserve (max_proposals + 6) * repeats before execution.
#'   Version probes are separate; cached panels can reduce actual job counts.
#' @param min_improvement Minimum absolute tuning and test RMSE improvement.
#' @param max_group_regression Maximum group RMSE regression relative to baseline.
#'   Margins are practical thresholds, not significance tests.
#' @param max_genotype_cells Maximum total input genotype cells, checked by streaming
#'   before genotype allocation. This is not a total RAM limit.
#' @param max_ld_pairs Maximum local candidate/scoring pairs summed over groups.
#' @param seed Nonnegative master seed. Caller RNG state is restored.
#' @param nthreads,memory_gb Shared Beagle and JVM resources. Jobs run sequentially.
#' @param beagle_jar,java Optional Beagle runtime paths.
#' @return A galick_panel_design inheriting from galick_result, with
#'   baseline_panel, candidate_panel, recommended_panel (NULL unless confirmed),
#'   metrics, by_group, search_history, selection, runs and per_marker_paths.
#'   Execution errors are collected; invalid designs raise
#'   galick_error_panel_design_invalid_input. An ok result need not have a winner.
#' @export
design_marker_panel <- function(truth_vcf, reference_vcf, sample_data,
                                candidate_markers, scoring_markers, n_markers, output_dir,
                                mandatory_markers = NULL, excluded_markers = NULL,
                                min_spacing_bp = 0L, ld_window_bp = 1e6,
                                map = NULL, extra_args = character(), field = c("DS", "GT"),
                                repeats = 3L, max_proposals = 12L, max_runs = 60L,
                                min_improvement = 0.01, max_group_regression = 0.01,
                                max_genotype_cells = 2e7, max_ld_pairs = 2e6,
                                seed = 2026L, nthreads = 1L, memory_gb = 4L,
                                beagle_jar = NULL, java = NULL) {
  the_call <- match.call()
  field <- match.arg(field)
  for (name in c("n_markers", "repeats", "max_runs", "nthreads", "memory_gb")) {
    panel_number(get(name), name, 1, TRUE)
  }
  for (name in c("max_proposals", "seed")) panel_number(get(name), name, 0, TRUE)
  for (name in c("min_spacing_bp", "ld_window_bp", "min_improvement", "max_group_regression")) {
    panel_number(get(name), name)
  }
  for (name in c("max_genotype_cells", "max_ld_pairs")) panel_number(get(name), name, 1)
  withr::local_seed(seed)
  if ((as.double(max_proposals) + 6) * repeats > max_runs) panel_abort("max_runs is below the reserved plan.")
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    panel_abort("output_dir must be new or empty.")
  }
  occupied <- file.exists(output_dir) && (
    !dir.exists(output_dir) || length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L
  )
  if (occupied) {
    panel_abort("output_dir must be new or empty.")
  }
  settings <- tryCatch(evaluation_candidates(list(a = list(map = map, extra_args = extra_args), b = list()), "a")[[1]],
    error = function(e) panel_abort(conditionMessage(e))
  )
  tracked <- lapply(unique(c(truth_vcf, reference_vcf, map)), function(path) {
    list(path = normalize_path_safe(path, must_exist = TRUE), sha256 = file_sha256(path))
  })
  input <- panel_inputs(
    truth_vcf, reference_vcf, sample_data, candidate_markers, scoring_markers,
    mandatory_markers, excluded_markers, max_genotype_cells
  )
  constraints <- panel_constraints(input, n_markers, min_spacing_bp)
  ld <- panel_ld(input, ld_window_bp, max_ld_pairs)
  baseline <- panel_seed(input, constraints)
  initial <- panel_seed(input, constraints, ld)
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) panel_abort("Cannot create output_dir.")
  output_dir <- normalize_path_safe(output_dir, must_exist = TRUE)
  config <- list(
    output_dir = output_dir, reference_vcf = normalize_path_safe(reference_vcf, must_exist = TRUE),
    map = settings$map, extra_args = settings$extra_args, field = field, repeats = repeats,
    max_proposals = max_proposals, max_runs = max_runs, min_improvement = min_improvement,
    max_group_regression = max_group_regression, min_spacing_bp = min_spacing_bp, ld_window_bp = ld_window_bp,
    max_genotype_cells = max_genotype_cells, max_ld_pairs = max_ld_pairs, n_ld_pairs = ld$n_pairs,
    nthreads = nthreads, memory_gb = memory_gb, seed = seed, beagle_jar = beagle_jar, java = java,
    quotas = as.list(constraints$quotas), mandatory_markers = panel_markers(input, input$mandatory),
    candidate_markers = panel_markers(input, input$candidates), removed_markers = panel_markers(input, input$removed)
  )
  evaluation_write_json(config, file.path(output_dir, "experiment.json"))
  tracked[[length(tracked) + 1L]] <- list(
    path = file.path(output_dir, "experiment.json"),
    sha256 = file_sha256(file.path(output_dir, "experiment.json"))
  )
  plans <- panel_plans(input, repeats, seed, output_dir)
  state <- panel_execute(input, constraints, ld, baseline, initial, plans, config, tracked)
  new_galick_panel_design(state, input, baseline, config, the_call)
}
