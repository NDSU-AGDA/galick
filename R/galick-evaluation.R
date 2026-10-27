#' @keywords internal
#' @noRd
evaluation_write_json <- function(value, path) {
  json <- jsonlite::toJSON(value,
    auto_unbox = TRUE, pretty = TRUE,
    null = "null", na = "null", digits = NA
  )
  write_exact_lines(json, path)
}

#' @keywords internal
#' @noRd
new_galick_evaluation <- function(state, experiment, output_dir, call) {
  record_path <- file.path(output_dir, "evaluation.json")
  warnings <- state$warnings
  if (any(state$metrics$n_missing > 0)) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "evaluation_missing_predictions",
      message = "Missing predictions remain in coverage counts and prevent selection.", where = list()
    )
  }
  object <- new_galick_result(
    errors = issues_tibble(state$errors), warnings = issues_tibble(warnings),
    summary = list(
      decision = state$selection$decision, candidate = state$selection$recommended_candidate,
      n_runs = length(state$runs), output_dir = output_dir, evaluation_record = record_path
    ), call = call
  )
  for (name in c("metrics", "by_group", "selection", "runs", "per_marker_paths")) object[[name]] <- state[[name]]
  object$design <- experiment
  class(object) <- c("galick_evaluation", class(object))
  for (name in c("metrics", "by_group")) {
    data.table::fwrite(object[[name]], file.path(output_dir, paste0(name, ".tsv")), sep = "\t", na = "NA")
  }
  files <- unique(c(
    unname(state$per_marker_paths), file.path(output_dir, c("metrics.tsv", "by_group.tsv")),
    state$selection$frozen_path
  ))
  artifacts <- lapply(files, function(path) list(path = path, sha256 = file_sha256(path)))
  evaluation_write_json(list(
    schema_version = "1.0.0", design = experiment, selection = state$selection,
    metrics = state$metrics, by_group = state$by_group, artifacts = artifacts,
    run_records = lapply(state$runs, function(run) {
      record_exists <- length(run$run_record) == 1L && !is.na(run$run_record) && file.exists(run$run_record)
      list(
        status = run$status, path = run$run_record, command = run$command,
        versions = run$version, sha256 = if (record_exists) file_sha256(run$run_record) else NULL
      )
    }), package_version = as.character(utils::packageVersion("galick")),
    r_version = as.character(getRversion()), errors = object$errors, warnings = object$warnings
  ), record_path)
  object
}

#' @export
print.galick_evaluation <- function(x, ...) {
  cli::cli_h1("Imputation evaluation: {x$selection$decision}")
  cli::cli_text("Runs: {x$summary$n_runs}; errors: {nrow(x$errors)}; warnings: {nrow(x$warnings)}")
  print(as.data.frame(x)[c("candidate", "partition", "status", "rmse", "completion_rate")], row.names = FALSE)
  for (i in seq_len(min(nrow(x$errors), 3L))) {
    cli::cli_alert_danger("{x$errors$code[i]}: {x$errors$message[i]}")
  }
  invisible(x)
}

#' @export
as.data.frame.galick_evaluation <- function(x, ...) {
  as.data.frame(x$metrics)
}
