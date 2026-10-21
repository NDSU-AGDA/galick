#' Compare Beagle configurations using separate tuning and test samples
#'
#' Runs paired sporadic-call masking or reference-based array thinning. Selection
#' uses tuning RMSE only and requires complete predictions. The frozen candidate
#' is compared with the baseline on test samples, without reranking. Thresholds
#' are practical margins, not significance tests. Repeated counts are prediction
#' attempts, not independent biological observations. Truth quality and biological
#' independence are the caller's responsibility; renamed duplicates are not found.
#'
#' @param truth_vcf VCF with diploid biallelic A/C/G/T SNP truth calls.
#' @param sample_data Data frame with unique `sample_id`, `partition` (tune/test),
#'   optional reporting `group`, and optional `block` that cannot cross partitions.
#'   Only listed samples are evaluated. All truth sample IDs are excluded from references.
#' @param candidates Two to four named lists with optional `ref`, `map`,
#'   `extra_args`. Allowed settings: ne, em, err, burnin, iterations, phase-states,
#'   imp-states, window, overlap. Fixed ne requires em=false.
#' @param baseline Name of the baseline candidate.
#' @param output_dir New or empty output directory; automatic resume is not supported.
#' @param design Sporadic missing calls or an observed-marker array.
#' @param observed_markers Array data frame with chrom, pos, ref, alt matching truth exactly.
#' @param fraction Sporadic masking fraction strictly between zero and one; ignored for array.
#' @param field Evaluate GT counts or expected DS; no fallback when DS is absent.
#' @param repeats Number of paired repeats, at least one. Array repeats change engine seeds,
#'   not the observed-marker panel. Samples run separately by partition, not one at a time.
#' @param seed Nonnegative integer master seed. Caller RNG state is restored.
#' @param min_improvement Minimum absolute RMSE reduction over baseline and runner-up.
#' @param max_group_regression Maximum allowed absolute group RMSE increase over baseline.
#' @param nthreads,memory_gb Shared positive integer Beagle/JVM resource settings.
#' @param max_runs Maximum Beagle computations (not version probes). Reserve the worst
#'   case `(length(candidates) + 2) * repeats` before launching any jobs.
#' @param beagle_jar,java Optional runtime paths passed to [run_beagle()].
#' @return A `beagler_evaluation` inheriting from `beagler_result`, with `metrics`,
#'   `by_group`, `selection`, `runs`, `design` and `per_marker_paths`. Missing
#'   predictions prevent selection. Execution/scoring failures are collected;
#'   invalid experiment designs raise `beagler_error_evaluation_invalid_input`.
#' @export
evaluate_imputation <- function(truth_vcf, sample_data, candidates, baseline, output_dir,
                                design = c("sporadic", "array"), observed_markers = NULL,
                                fraction = 0.1, field = c("GT", "DS"), repeats = 3L,
                                seed = 2026L, min_improvement = 0.01,
                                max_group_regression = 0.01, nthreads = 1L,
                                memory_gb = 4L, max_runs = 20L, beagle_jar = NULL, java = NULL) {
  the_call <- match.call()
  design <- match.arg(design)
  field <- match.arg(field)
  for (name in c("repeats", "nthreads", "memory_gb", "max_runs")) {
    evaluation_number(get(name), name, 1, .Machine$integer.max, integer = TRUE)
  }
  evaluation_number(seed, "seed", 0, .Machine$integer.max, integer = TRUE)
  evaluation_number(min_improvement, "min_improvement", 0)
  evaluation_number(max_group_regression, "max_group_regression", 0)
  if (design == "sporadic") {
    evaluation_number(fraction, "fraction", 0, 1)
    if (fraction == 0 || fraction == 1) evaluation_abort("fraction must be strictly between zero and one.")
  }
  candidates <- evaluation_candidates(candidates, baseline)
  if ((length(candidates) + 2) * as.double(repeats) > max_runs) {
    evaluation_abort("max_runs cannot accommodate the worst-case tuning and confirmation plan.")
  }
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir) ||
        (file.exists(output_dir) && (!dir.exists(output_dir) ||
                                       length(list.files(output_dir, all.files = TRUE, no.. = TRUE))))) {
    evaluation_abort("output_dir must be new or empty.")
  }
  validate_vcf_path(truth_vcf, "truth_vcf")
  if (!validate_vcf(truth_vcf, target = "beagle")$ok) evaluation_abort("Truth failed VCF validation.")
  truth <- read_genotype_vcf(truth_vcf)
  if (!nrow(truth$gt) || !ncol(truth$gt) || any(!truth$fixed$REF %in% c("A", "C", "G", "T")) ||
        any(!truth$fixed$ALT %in% c("A", "C", "G", "T"))) {
    evaluation_abort("Evaluation requires nonempty diploid biallelic A/C/G/T SNP truth.")
  }
  samples <- evaluation_samples(sample_data, truth$samples)
  plans <- evaluation_plan(truth, samples, design, observed_markers, fraction, repeats, seed)
  preflight_warnings <- list()
  for (candidate in candidates) {
    if (design == "array" && is.null(candidate$ref)) evaluation_abort("Every array candidate requires a reference.")
    if (!is.null(candidate$ref)) {
      checked <- validate_reference_match(truth_vcf, candidate$ref, require_phased = TRUE)
      if (!checked$ok || checked$summary$n_markers_shared == 0L) evaluation_abort("Reference preflight failed.")
      for (i in seq_len(nrow(checked$warnings))) {
        preflight_warnings[[length(preflight_warnings) + 1L]] <- list(
          code = checked$warnings$code[i], message = checked$warnings$message[i],
          where = c(checked$warnings$where[[i]], list(reference = candidate$ref))
        )
      }
      reference <- read_genotype_vcf(candidate$ref)
      observed <- plans[[1]]$observed
      common <- vcf_variant_keys(truth) %in% vcf_variant_keys(reference)
      if (!setequal(truth$fixed$CHROM, truth$fixed$CHROM[common & observed])) {
        evaluation_abort("Each evaluated chromosome needs shared observed markers in every reference.")
      }
    }
  }
  paths <- unique(c(truth_vcf, unlist(lapply(candidates, function(x) c(x$ref, x$map)), use.names = FALSE)))
  input_hashes <- lapply(paths, function(path) {
    list(path = normalize_path_safe(path, must_exist = TRUE), sha256 = file_sha256(path))
  })
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    evaluation_abort("Cannot create output_dir.")
  }
  output_dir <- normalize_path_safe(output_dir, must_exist = TRUE)
  assert_writable(output_dir)
  for (name in names(plans)) {
    plan <- plans[[name]]
    plan$input_path <- file.path(output_dir, paste0(name, ".vcf"))
    plan$manifest_path <- file.path(output_dir, paste0(name, "_heldout.tsv"))
    write_genotype_vcf(evaluation_input(truth, plan, design), plan$input_path)
    data.table::fwrite(evaluation_manifest(truth, plan), plan$manifest_path, sep = "\t")
    plan$input_sha256 <- file_sha256(plan$input_path)
    plan$manifest_sha256 <- file_sha256(plan$manifest_path)
    plans[[name]] <- plan
  }
  state <- evaluation_execute(
    truth, plans, candidates, baseline, output_dir, field,
    nthreads, memory_gb, beagle_jar, java, min_improvement, max_group_regression
  )
  state$warnings <- preflight_warnings
  tracked <- c(input_hashes, unlist(lapply(plans, function(x) {
    list(
      list(path = x$input_path, sha256 = x$input_sha256),
      list(path = x$manifest_path, sha256 = x$manifest_sha256)
    )
  }), recursive = FALSE))
  if (any(!vapply(tracked, function(x) {
    isTRUE(tryCatch(identical(file_sha256(x$path), x$sha256), error = function(e) FALSE))
  }, logical(1)))) {
    state$errors[[length(state$errors) + 1L]] <- list(
      code = "evaluation_input_changed",
      message = "An input or held-out design changed during evaluation; no recommendation is valid.", where = list()
    )
    state$selection$decision <- "incomplete"
    state$selection$recommended_candidate <- NA_character_
  }
  manifests <- lapply(plans, function(x) {
    x[c(
      "partition", "replicate", "seed", "input_path",
      "input_sha256", "manifest_path", "manifest_sha256"
    )]
  })
  experiment <- list(
    schema_version = "1.0.0", design = design, field = field,
    baseline = baseline, candidates = candidates, samples = samples, inputs = input_hashes,
    manifests = manifests, repeats = repeats, seed = seed, fraction = if (design == "sporadic") fraction else NULL,
    min_improvement = min_improvement, max_group_regression = max_group_regression,
    max_runs = max_runs, nthreads = nthreads, memory_gb = memory_gb
  )
  new_beagler_evaluation(state, experiment, output_dir, the_call)
}

#' @keywords internal
#' @noRd
evaluation_execute <- function(truth, plans, candidates, baseline, output_dir, field,
                               nthreads, memory_gb, beagle_jar, java,
                               min_improvement, max_group_regression) {
  runs <- list()
  errors <- list()
  metrics <- list()
  groups <- list()
  marker_paths <- character()
  execute_partition <- function(partition, ids) {
    for (id in ids) {
      candidate <- candidates[[id]]
      all_cells <- list()
      runtime <- 0
      success <- TRUE
      for (name in names(plans)) {
        plan <- plans[[name]]
        if (plan$partition != partition) next
        job_id <- paste(id, name, sep = "_")
        directory <- file.path(output_dir, "jobs", job_id)
        dir.create(directory, recursive = TRUE)
        attempt <- tryCatch(
          {
            if (!identical(file_sha256(plan$input_path), plan$input_sha256) ||
                  !identical(file_sha256(plan$manifest_path), plan$manifest_sha256)) {
              evaluation_abort("Paired input or eligibility manifest was modified.")
            }
            run <- run_beagle(plan$input_path, file.path(directory, "imputed"),
              ref = candidate$ref,
              map = candidate$map, mode = "impute", nthreads = nthreads, memory_gb = memory_gb,
              seed = plan$seed, beagle_jar = beagle_jar, java = java,
              extra_args = candidate$extra_args, run_dir = directory
            )
            runs[[job_id]] <<- run
            runtime <- runtime + if (is.finite(run$runtime_sec)) run$runtime_sec else 0
            if (run$status != 0L || length(run$output_paths) != 1L) evaluation_abort("Beagle job failed.")
            evaluation_cells(truth, read_genotype_vcf(run$output_paths), plan, field)
          },
          error = function(e) {
            success <<- FALSE
            errors[[length(errors) + 1L]] <<- list(
              code = "evaluation_job_failed",
              message = conditionMessage(e), where = list(
                candidate = id, partition = partition,
                replicate = plan$replicate, job_dir = directory
              )
            )
            evaluation_missing_cells(truth, plan)
          }
        )
        all_cells[[name]] <- attempt
      }
      cells <- data.table::rbindlist(all_cells)
      tables <- evaluation_tables(cells, truth)
      key <- paste(id, partition, sep = "_")
      metrics[[key]] <<- cbind(tibble::tibble(
        candidate = id, partition = partition,
        status = if (success) "ok" else "failed", runtime_sec = runtime
      ), tables$metrics)
      groups[[key]] <<- cbind(tibble::tibble(candidate = id, partition = partition), tables$by_group)
      marker_path <- file.path(output_dir, paste0(key, "_markers.tsv"))
      data.table::fwrite(tables$per_marker, marker_path, sep = "\t", na = "NA")
      marker_paths[key] <<- marker_path
    }
  }
  execute_partition("tune", names(candidates))
  tuning <- as.data.frame(data.table::rbindlist(metrics))
  tuning_groups <- as.data.frame(data.table::rbindlist(groups))
  frozen <- evaluation_select(
    tuning, tuning_groups, names(candidates), baseline,
    min_improvement, max_group_regression
  )
  selection_path <- file.path(output_dir, "selection.json")
  evaluation_write_json(list(
    schema_version = "1.0.0", baseline = baseline,
    tuning = frozen, min_improvement = min_improvement, max_group_regression = max_group_regression,
    tuning_metrics = tuning, tuning_groups = tuning_groups
  ), selection_path)
  selection_sha256 <- file_sha256(selection_path)
  if (frozen$decision != "incomplete") {
    ids <- unique(c(baseline, if (!is.na(frozen$candidate)) frozen$candidate))
    execute_partition("test", ids)
  }
  metrics <- tibble::as_tibble(data.table::rbindlist(metrics))
  groups <- tibble::as_tibble(data.table::rbindlist(groups))
  decision <- if (frozen$decision == "incomplete") {
    "incomplete"
  } else {
    evaluation_confirm(
      metrics[metrics$partition == "test", ], groups[groups$partition == "test", ],
      frozen$candidate, baseline, min_improvement, max_group_regression
    )
  }
  selection_intact <- tryCatch(identical(file_sha256(selection_path), selection_sha256), error = function(e) FALSE)
  if (!selection_intact) {
    errors[[length(errors) + 1L]] <- list(
      code = "evaluation_selection_changed",
      message = "The frozen selection record changed during confirmation.", where = list()
    )
    decision <- "incomplete"
  }
  list(
    metrics = metrics, by_group = groups, runs = runs, errors = errors,
    per_marker_paths = marker_paths, selection = list(
      decision = decision,
      candidate = frozen$candidate,
      recommended_candidate = if (decision == "confirmed") frozen$candidate else NA_character_,
      baseline = baseline, frozen_path = selection_path, frozen_sha256 = selection_sha256
    )
  )
}

#' @keywords internal
#' @noRd
evaluation_missing_cells <- function(truth, plan) {
  index <- which(plan$eligible, arr.ind = TRUE)
  values <- genotype_counts(truth)[, plan$columns, drop = FALSE]
  tibble::tibble(
    marker = index[, 1], sample_id = plan$samples$sample_id[index[, 2]],
    group = plan$samples$group[index[, 2]], truth = values[index], prediction = NA_real_, gt = NA_real_
  )
}
