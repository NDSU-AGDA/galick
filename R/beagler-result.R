#' beagler_result S3 class (internal)
#' @name beagler-result
#' @keywords internal
NULL

#' Construct a beagler_result
#'
#' @param errors A data frame / tibble with columns `code`, `message`, and
#'   optional `where` (list-column). `NULL` is treated as an empty issues
#'   table.
#' @param warnings Same shape as `errors`.
#' @param summary A function-specific named list.
#' @param call The captured call expression (e.g. `match.call()`).
#' @param ok Optional logical override; if `NULL`, derived as `nrow(errors) == 0`.
#'
#' @return A `beagler_result` S3 object.
#' @keywords internal
#' @noRd
new_beagler_result <- function(errors = NULL,
                               warnings = NULL,
                               summary = list(),
                               call = NULL,
                               ok = NULL) {
  errors <- normalize_issues(errors)
  warnings <- normalize_issues(warnings)

  if (is.null(ok)) {
    ok <- nrow(errors) == 0L
  }

  metadata <- list(
    call = call,
    package_version = utils::packageVersion("beagler"),
    time = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
  )

  structure(
    list(
      ok = as.logical(ok),
      errors = errors,
      warnings = warnings,
      summary = as.list(summary),
      metadata = metadata
    ),
    class = c("beagler_result", "list")
  )
}

#' Coerce an issues input to the canonical 3-column tibble
#' @keywords internal
#' @noRd
normalize_issues <- function(x) {
  empty <- tibble::tibble(
    code = character(),
    message = character(),
    where = list()
  )
  if (is.null(x) || (inherits(x, "data.frame") && nrow(x) == 0L)) {
    return(empty)
  }
  if (!inherits(x, "data.frame")) {
    cli::cli_abort(
      "Issues must be a data.frame or tibble; got {.cls {class(x)[1]}}.",
      class = "beagler_error_invalid_issues"
    )
  }
  if (!all(c("code", "message") %in% names(x))) {
    cli::cli_abort(
      "Issues must have at least {.field code} and {.field message} columns.",
      class = "beagler_error_invalid_issues"
    )
  }
  if (!"where" %in% names(x)) {
    x$where <- replicate(nrow(x), list(), simplify = FALSE)
  }
  tibble::as_tibble(x[c("code", "message", "where")])
}

#' @export
print.beagler_result <- function(x, ...) {
  header <- if (x$ok) "beagler_result: OK" else "beagler_result: NOT OK"
  cli::cli_h1(header)
  cli::cli_text(
    "Errors: {nrow(x$errors)}; Warnings: {nrow(x$warnings)}"
  )

  if (nrow(x$errors) > 0L) {
    cli::cli_h3("Errors (top 3)")
    for (i in seq_len(min(3L, nrow(x$errors)))) {
      cli::cli_alert_danger(
        "{x$errors$code[i]}: {x$errors$message[i]}"
      )
    }
    if (nrow(x$errors) > 3L) {
      cli::cli_alert_info("... and {nrow(x$errors) - 3L} more.")
    }
  }

  if (nrow(x$warnings) > 0L) {
    cli::cli_h3("Warnings (top 3)")
    for (i in seq_len(min(3L, nrow(x$warnings)))) {
      cli::cli_alert_warning(
        "{x$warnings$code[i]}: {x$warnings$message[i]}"
      )
    }
    if (nrow(x$warnings) > 3L) {
      cli::cli_alert_info("... and {nrow(x$warnings) - 3L} more.")
    }
  }

  if (length(x$summary) > 0L) {
    cli::cli_h3("Summary")
    for (k in names(x$summary)) {
      cli::cli_text("{k}: {x$summary[[k]]}")
    }
  }

  invisible(x)
}

#' @export
summary.beagler_result <- function(object, ...) {
  print(object)
  if (nrow(object$errors) > 3L) {
    cli::cli_h3("All errors")
    print(object$errors)
  }
  if (nrow(object$warnings) > 3L) {
    cli::cli_h3("All warnings")
    print(object$warnings)
  }
  invisible(object)
}

#' @export
format.beagler_result <- function(x, ...) {
  sprintf(
    "<beagler_result: %s, %d errors, %d warnings>",
    if (x$ok) "ok" else "not ok",
    nrow(x$errors),
    nrow(x$warnings)
  )
}

#' @export
as.data.frame.beagler_result <- function(x, ...) {
  e <- x$errors
  w <- x$warnings

  if (nrow(e) == 0L && nrow(w) == 0L) {
    return(data.frame(
      severity = character(),
      code = character(),
      message = character(),
      where = character(),
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(e) > 0L) e$severity <- "error"
  if (nrow(w) > 0L) w$severity <- "warning"
  combined <- rbind(e, w)
  combined <- combined[, c("severity", "code", "message", "where"), drop = FALSE]

  flat <- flatten_where(combined$where)
  if (!is.null(flat)) {
    combined$where <- NULL
    combined <- cbind(combined, flat)
  } else {
    combined$where <- vapply(combined$where, function(w) {
      if (length(w) == 0L) return("")
      jsonlite::toJSON(w, auto_unbox = TRUE)
    }, character(1))
  }

  as.data.frame(combined, stringsAsFactors = FALSE)
}

#' Try to flatten a `where` list-column into scalar `where_*` columns
#'
#' Returns `NULL` if rows disagree on key set or any value is non-scalar.
#' Callers fall back to JSON stringification.
#'
#' @keywords internal
#' @noRd
flatten_where <- function(where_list) {
  if (length(where_list) == 0L) return(NULL)
  all_keys <- unique(unlist(lapply(where_list, names)))
  if (length(all_keys) == 0L) return(NULL)

  same_keys <- vapply(where_list, function(w) setequal(names(w), all_keys), logical(1))
  if (!all(same_keys)) return(NULL)

  all_scalar <- vapply(where_list, function(w) {
    all(vapply(w, function(v) length(v) == 1L && is.atomic(v), logical(1)))
  }, logical(1))
  if (!all(all_scalar)) return(NULL)

  cols <- lapply(all_keys, function(k) {
    vapply(where_list, function(w) {
      val <- w[[k]]
      if (is.null(val)) NA else as.character(val)
    }, character(1))
  })
  names(cols) <- paste0("where_", all_keys)
  as.data.frame(cols, stringsAsFactors = FALSE)
}
