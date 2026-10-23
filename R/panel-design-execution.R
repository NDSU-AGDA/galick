#' Fixed scoring masks; panels affect engine inputs, never denominators.
#' @keywords internal
#' @noRd
panel_plans <- function(input, repeats, seed, output_dir) {
  plans <- list()
  for (partition in c("design", "tune", "test")) {
    metadata <- input$samples[input$samples$partition == partition, , drop = FALSE]
    columns <- match(metadata$sample_id, input$truth$samples)
    eligible <- matrix(FALSE, nrow(input$counts), length(columns))
    eligible[input$scoring, ] <- !is.na(input$counts[input$scoring, columns, drop = FALSE])
    for (replicate in seq_len(repeats)) {
      name <- paste(partition, replicate, sep = "_")
      offset <- match(partition, c("design", "tune", "test")) - 1L
      plan <- list(
        partition = partition, replicate = replicate,
        seed = as.integer((as.double(seed) + offset * repeats + replicate - 1) %% .Machine$integer.max),
        columns = columns, samples = metadata, eligible = eligible
      )
      path <- file.path(output_dir, paste0(name, "_scoring.tsv"))
      data.table::fwrite(evaluation_manifest(input$truth, plan), path, sep = "\t")
      plan$manifest_path <- path
      plan$manifest_sha256 <- file_sha256(path)
      plans[[name]] <- plan
    }
  }
  plans
}

#' @keywords internal
#' @noRd
panel_integrity <- function(tracked) {
  all(vapply(tracked, function(x) {
    isTRUE(tryCatch(identical(file_sha256(x$path), x$sha256), error = function(e) FALSE))
  }, logical(1)))
}

#' @keywords internal
#' @noRd
panel_safe <- function(candidate, baseline, margin) {
  index <- match(baseline$groups$group, candidate$groups$group)
  candidate$metrics$completion_rate == 1 && is.finite(candidate$metrics$rmse) && !anyNA(index) &&
    all(is.finite(candidate$groups$rmse[index])) &&
    all(candidate$groups$rmse[index] - baseline$groups$rmse <= margin)
}

#' @keywords internal
#' @noRd
panel_execute <- function(input, constraints, ld, baseline, initial, plans, config, tracked) {
  state <- new.env(parent = emptyenv())
  state$runs <- list()
  state$metrics <- list()
  state$groups <- list()
  state$marker_paths <- character()
  state$history <- list()
  state$cache <- list()
  state$errors <- list()
  state$tracked <- tracked
  state$candidate <- baseline
  state$selection <- list(decision = "incomplete", candidate = NA_character_, recommended_candidate = NA_character_)
  for (plan in plans) {
    state$tracked[[length(state$tracked) + 1L]] <- list(
      path = plan$manifest_path, sha256 = plan$manifest_sha256
    )
  }
  freeze <- function(value, name) {
    path <- file.path(config$output_dir, name)
    evaluation_write_json(value, path)
    state$tracked[[length(state$tracked) + 1L]] <- list(path = path, sha256 = file_sha256(path))
    path
  }
  score <- function(rows, partition, id) {
    key <- paste(partition, paste(rows, collapse = ","), sep = ":")
    if (!panel_integrity(state$tracked)) stop("Input or frozen artifact changed during panel design.")
    if (!is.null(state$cache[[key]])) {
      return(state$cache[[key]])
    }
    cells <- list()
    runtime <- 0
    for (name in names(plans)) {
      plan <- plans[[name]]
      if (plan$partition != partition) next
      if (length(state$runs) >= config$max_runs) stop("Panel job budget exhausted.")
      job_id <- paste(id, name, sep = "_")
      directory <- file.path(config$output_dir, "jobs", job_id)
      dir.create(directory, recursive = TRUE)
      plan$observed <- seq_len(nrow(input$truth$gt)) %in% rows
      path <- file.path(directory, "input.vcf")
      write_genotype_vcf(evaluation_input(input$truth, plan, "array"), path)
      state$tracked[[length(state$tracked) + 1L]] <- list(path = path, sha256 = file_sha256(path))
      run <- run_beagle(path, file.path(directory, "imputed"),
        ref = config$reference_vcf,
        map = config$map, mode = "impute", nthreads = config$nthreads, memory_gb = config$memory_gb,
        seed = plan$seed, beagle_jar = config$beagle_jar, java = config$java,
        extra_args = config$extra_args, run_dir = directory
      )
      state$runs[[job_id]] <- run
      if (!isTRUE(run$status == 0L) || length(run$output_paths) != 1L) stop("Beagle panel job failed.")
      if (!panel_integrity(state$tracked)) stop("Input or frozen artifact changed during a panel job.")
      prediction <- read_genotype_vcf(run$output_paths)
      cells[[name]] <- evaluation_cells(input$truth, prediction, plan, config$field)
      runtime <- runtime + if (is.finite(run$runtime_sec)) run$runtime_sec else 0
      state$tracked[[length(state$tracked) + 1L]] <- list(
        path = run$output_paths, sha256 = file_sha256(run$output_paths)
      )
    }
    tables <- evaluation_tables(data.table::rbindlist(cells), input$truth)
    metrics <- cbind(
      tibble::tibble(candidate = id, partition = partition, status = "ok", runtime_sec = runtime),
      tables$metrics
    )
    groups <- cbind(tibble::tibble(candidate = id, partition = partition), tables$by_group)
    marker_path <- file.path(config$output_dir, paste0(id, "_", partition, "_markers.tsv"))
    data.table::fwrite(tables$per_marker, marker_path, sep = "\t", na = "NA")
    state$metrics[[length(state$metrics) + 1L]] <- metrics
    state$groups[[length(state$groups) + 1L]] <- groups
    state$marker_paths[paste(id, partition, sep = "_")] <- marker_path
    result <- list(metrics = metrics, groups = groups)
    state$cache[[key]] <- result
    result
  }
  record <- function(rows, id, scored, accepted) {
    path <- file.path(config$output_dir, paste0(id, "_panel.tsv"))
    data.table::fwrite(panel_markers(input, rows), path, sep = "\t")
    state$tracked[[length(state$tracked) + 1L]] <- list(path = path, sha256 = file_sha256(path))
    state$history[[length(state$history) + 1L]] <- tibble::tibble(
      panel = id, proxy = panel_proxy(rows, ld), rmse = scored$metrics$rmse,
      completion_rate = scored$metrics$completion_rate, accepted = accepted, path = path
    )
  }
  tryCatch(
    {
      base <- score(baseline, "design", "baseline")
      record(baseline, "baseline", base, TRUE)
      if (base$metrics$completion_rate != 1) stop("Baseline has missing design predictions.")
      current <- base
      incumbent <- baseline
      visited <- paste(baseline, collapse = ",")
      consider <- function(rows, id) {
        scored <- score(rows, "design", id)
        accepted <- panel_safe(scored, base, config$max_group_regression) &&
          scored$metrics$rmse < current$metrics$rmse
        record(rows, id, scored, accepted)
        if (accepted) {
          incumbent <<- rows
          current <<- scored
        }
        visited <<- c(visited, paste(rows, collapse = ","))
      }
      if (!identical(initial, baseline)) consider(initial, "ld_seed")
      for (i in seq_len(config$max_proposals)) {
        proposal <- panel_proposal(input, constraints, ld, incumbent, visited)
        if (is.null(proposal)) break
        consider(proposal, paste0("proposal", i))
      }
      state$candidate <- incumbent
      freeze(list(
        panel = panel_markers(input, incumbent), baseline = panel_markers(input, baseline),
        extra_args = config$extra_args, field = config$field,
        min_improvement = config$min_improvement, max_group_regression = config$max_group_regression
      ), "challenger.json")
      tune_base <- score(baseline, "tune", "baseline")
      if (identical(incumbent, baseline)) {
        selected <- list(
          decision = if (tune_base$metrics$completion_rate == 1) "no_clear_winner" else "incomplete",
          candidate = NA_character_
        )
      } else {
        tune_candidate <- score(incumbent, "tune", "challenger")
        selected <- evaluation_select(
          rbind(tune_base$metrics, tune_candidate$metrics),
          rbind(tune_base$groups, tune_candidate$groups), c("baseline", "challenger"), "baseline",
          config$min_improvement, config$max_group_regression
        )
      }
      freeze(selected, "selection.json")
      decision <- selected$decision
      if (decision != "incomplete") {
        test_base <- score(baseline, "test", "baseline")
        metrics <- test_base$metrics
        groups <- test_base$groups
        if (!is.na(selected$candidate)) {
          test_candidate <- score(incumbent, "test", "challenger")
          metrics <- rbind(metrics, test_candidate$metrics)
          groups <- rbind(groups, test_candidate$groups)
        }
        decision <- evaluation_confirm(
          metrics, groups, selected$candidate, "baseline",
          config$min_improvement, config$max_group_regression
        )
      }
      state$selection <- list(
        decision = decision, candidate = selected$candidate,
        recommended_candidate = if (decision == "confirmed") "challenger" else NA_character_
      )
    },
    error = function(e) {
      state$errors[[length(state$errors) + 1L]] <- list(
        code = "panel_design_execution_failed",
        message = conditionMessage(e), where = list()
      )
    }
  )
  if (!panel_integrity(state$tracked)) {
    state$errors[[length(state$errors) + 1L]] <- list(
      code = "panel_design_input_changed",
      message = "An input or frozen artifact changed during panel design.", where = list()
    )
    state$selection$decision <- "incomplete"
    state$selection$recommended_candidate <- NA_character_
  }
  state
}
