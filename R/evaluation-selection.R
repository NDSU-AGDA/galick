#' Choose only among complete paired evaluations
#' @keywords internal
#' @noRd
evaluation_select <- function(metrics, groups, ids, baseline, min_improvement, max_group_regression) {
  complete <- nrow(metrics) == length(ids) && setequal(metrics$candidate, ids) &&
    all(metrics$status == "ok" & metrics$completion_rate == 1 & is.finite(metrics$rmse))
  if (!complete) {
    return(list(decision = "incomplete", candidate = NA_character_))
  }
  scores <- stats::setNames(metrics$rmse, metrics$candidate)
  ordered <- names(sort(scores))
  best <- ordered[1]
  peers <- setdiff(ids, best)
  base_groups <- groups[groups$candidate == baseline, , drop = FALSE]
  best_groups <- groups[groups$candidate == best, , drop = FALSE]
  paired <- match(base_groups$group, best_groups$group)
  safe <- !anyNA(paired) && all(is.finite(best_groups$rmse[paired])) &&
    all(best_groups$rmse[paired] - base_groups$rmse <= max_group_regression)
  clear <- best != baseline && all(scores[peers] - scores[best] >= min_improvement) &&
    all(scores[peers] > scores[best]) && safe
  list(
    decision = if (clear) "selected" else "no_clear_winner",
    candidate = if (clear) best else NA_character_
  )
}

#' @keywords internal
#' @noRd
evaluation_confirm <- function(metrics, groups, selected, baseline, min_improvement, max_group_regression) {
  if (is.na(selected)) {
    if (any(metrics$status != "ok") || any(metrics$completion_rate != 1)) {
      return("incomplete")
    }
    return("no_clear_winner")
  }
  choice <- evaluation_select(
    metrics, groups, c(baseline, selected), baseline,
    min_improvement, max_group_regression
  )
  if (choice$decision == "incomplete") {
    return("incomplete")
  }
  if (identical(choice$candidate, selected)) "confirmed" else "not_confirmed"
}
