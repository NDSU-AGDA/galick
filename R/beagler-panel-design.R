#' @keywords internal
#' @noRd
new_beagler_panel_design <- function(state, input, baseline, config, call) {
  output <- config$output_dir
  metrics <- tibble::as_tibble(data.table::rbindlist(state$metrics))
  groups <- tibble::as_tibble(data.table::rbindlist(state$groups))
  history <- tibble::as_tibble(data.table::rbindlist(state$history))
  if (!ncol(metrics)) {
    metrics <- tibble::tibble(
      candidate = character(), partition = character(),
      status = character(), runtime_sec = numeric(), rmse = numeric(), completion_rate = numeric()
    )
  }
  if (!ncol(groups)) {
    groups <- tibble::tibble(
      candidate = character(), partition = character(),
      group = character(), rmse = numeric(), completion_rate = numeric()
    )
  }
  if (!ncol(history)) {
    history <- tibble::tibble(
      panel = character(), proxy = numeric(), rmse = numeric(),
      completion_rate = numeric(), accepted = logical(), path = character()
    )
  }
  paths <- file.path(output, c("metrics.tsv", "by_group.tsv", "search_history.tsv"))
  for (i in seq_along(paths)) data.table::fwrite(list(metrics, groups, history)[[i]], paths[i], sep = "\t", na = "NA")
  object <- new_beagler_result(
    errors = issues_tibble(state$errors), warnings = input$warnings,
    summary = list(
      decision = state$selection$decision, n_markers = length(baseline), n_runs = length(state$runs),
      n_ineligible = length(input$removed), output_dir = output
    ), call = call
  )
  object$baseline_panel <- panel_markers(input, baseline)
  object$candidate_panel <- panel_markers(input, state$candidate)
  object["recommended_panel"] <- list(
    if (state$selection$decision == "confirmed") object$candidate_panel else NULL
  )
  object$metrics <- metrics
  object$by_group <- groups
  object$search_history <- history
  object$selection <- state$selection
  object$runs <- state$runs
  object$per_marker_paths <- state$marker_paths
  for (name in c("baseline_panel", "candidate_panel")) {
    path <- file.path(output, paste0(name, ".tsv"))
    data.table::fwrite(object[[name]], path, sep = "\t")
    paths <- c(paths, path)
  }
  config$schema_version <- "1.0.0"
  config$samples <- input$samples
  config$scoring_markers <- panel_markers(input, input$scoring)
  config$artifacts <- c(state$tracked, lapply(
    c(paths, state$marker_paths),
    function(path) list(path = path, sha256 = file_sha256(path))
  ))
  config$selection <- state$selection
  config$errors <- object$errors
  config$run_records <- lapply(state$runs, function(run) {
    list(
      command = run$command, status = run$status,
      versions = run$version, path = run$run_record,
      sha256 = if (length(run$run_record) == 1L && !is.na(run$run_record) && file.exists(run$run_record)) {
        file_sha256(run$run_record)
      } else {
        NULL
      }
    )
  })
  config$package_version <- as.character(utils::packageVersion("beagler"))
  config$r_version <- as.character(getRversion())
  evaluation_write_json(config, file.path(output, "panel_design.json"))
  class(object) <- c("beagler_panel_design", class(object))
  object
}

#' @export
print.beagler_panel_design <- function(x, ...) {
  cli::cli_h1("Marker panel design: {x$selection$decision}")
  cli::cli_text("Markers: {x$summary$n_markers}; runs: {x$summary$n_runs}; errors: {nrow(x$errors)}")
  if (nrow(x$metrics)) {
    print(as.data.frame(x)[c("candidate", "partition", "rmse", "completion_rate")], row.names = FALSE)
  }
  invisible(x)
}

#' @export
as.data.frame.beagler_panel_design <- function(x, ...) {
  as.data.frame(x$metrics)
}
