#' Align predictions against an immutable evaluation design
#' @keywords internal
#' @noRd
evaluation_cells <- function(truth, prediction, plan, field) {
  invalid <- comparison_contract_error(truth, prediction)
  if (!is.null(invalid)) evaluation_abort(invalid$message)
  row <- match(vcf_variant_keys(truth), vcf_variant_keys(prediction))
  column <- match(plan$samples$sample_id, prediction$samples)
  gt <- matrix(NA_real_, nrow(truth$gt), length(column))
  present_rows <- which(!is.na(row))
  present_columns <- which(!is.na(column))
  gt[present_rows, present_columns] <- genotype_counts(prediction)[
    row[present_rows], column[present_columns],
    drop = FALSE
  ]
  values <- gt
  if (field == "DS") {
    tokens <- prediction$format_fields$DS
    if (is.null(tokens)) evaluation_abort("Requested DS is absent; GT fallback is forbidden.")
    dosage <- suppressWarnings(matrix(as.numeric(tokens), nrow(tokens), ncol(tokens)))
    if (any(tokens != "." & (!is.finite(dosage) | dosage < 0 | dosage > 2))) {
      evaluation_abort("DS must be missing or finite and between zero and two.")
    }
    values[, ] <- NA_real_
    values[present_rows, present_columns] <- dosage[row[present_rows], column[present_columns], drop = FALSE]
  }
  index <- which(plan$eligible, arr.ind = TRUE)
  true_values <- genotype_counts(truth)[, plan$columns, drop = FALSE]
  tibble::tibble(
    marker = index[, 1], sample_id = plan$samples$sample_id[index[, 2]],
    group = plan$samples$group[index[, 2]], truth = true_values[index],
    prediction = values[index], gt = gt[index]
  )
}

#' @keywords internal
#' @noRd
evaluation_metrics <- function(cells) {
  predicted <- !is.na(cells$prediction)
  called <- !is.na(cells$gt)
  squared_error <- sum((cells$prediction[predicted] - cells$truth[predicted])^2)
  r2 <- vapply(split(cells, cells$marker), function(x) {
    keep <- !is.na(x$prediction)
    comparison_r2(x$prediction[keep], x$truth[keep])
  }, numeric(1))
  tibble::tibble(
    n_eligible = nrow(cells), n_predicted = sum(predicted),
    n_missing = sum(!predicted), n_gt_compared = sum(called),
    n_concordant = sum(called & cells$gt == cells$truth, na.rm = TRUE),
    squared_error = squared_error,
    completion_rate = sum(predicted) / nrow(cells),
    concordance_rate = sum(called & cells$gt == cells$truth, na.rm = TRUE) / nrow(cells),
    rmse = if (any(predicted)) sqrt(squared_error / sum(predicted)) else NA_real_,
    dosage_r2 = if (any(is.finite(r2))) mean(r2[is.finite(r2)]) else NA_real_,
    n_undefined_r2 = sum(!is.finite(r2))
  )
}

#' @keywords internal
#' @noRd
evaluation_tables <- function(cells, truth) {
  groups <- lapply(split(cells, cells$group), function(x) {
    cbind(tibble::tibble(group = x$group[1]), evaluation_metrics(x))
  })
  markers <- lapply(split(cells, cells$marker), function(x) {
    i <- x$marker[1]
    cbind(tibble::tibble(
      chrom = truth$fixed$CHROM[i], pos = truth$fixed$POS[i],
      ref = truth$fixed$REF[i], alt = truth$fixed$ALT[i]
    ), evaluation_metrics(x))
  })
  list(
    metrics = evaluation_metrics(cells), by_group = data.table::rbindlist(groups),
    per_marker = data.table::rbindlist(markers)
  )
}
