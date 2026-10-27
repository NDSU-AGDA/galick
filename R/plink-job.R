#' Run a PLINK bridge job without exposing stale or partial output
#' @keywords internal
#' @noRd
run_plink_bridge <- function(plink, args, outputs) {
  parent <- normalize_path_safe(dirname(outputs[1]), TRUE)
  stage <- tempfile("galick_plink_", tmpdir = parent)
  if (!dir.create(stage)) {
    cli::cli_abort("Cannot create a PLINK staging directory.", class = "galick_error_plink_write_failed")
  }
  keep_stage <- FALSE
  on.exit(if (!keep_stage) unlink(stage, recursive = TRUE), add = TRUE)
  prefix <- file.path(stage, "output")
  extensions <- if (length(outputs) == 3L) c(".bed", ".bim", ".fam") else ".vcf"
  staged <- paste0(prefix, extensions)
  result <- run_cli(plink, c(args, "--out", prefix))
  valid <- result$status == 0L && all(file.exists(staged)) &&
    all(file.info(staged)$size > 0)
  if (valid && identical(extensions, ".vcf")) {
    valid <- tryCatch(validate_vcf(staged, target = "general")$ok, error = function(e) FALSE)
  }
  if (!valid) {
    if (result$status == 0L) result$status <- 1L
    return(result)
  }

  # Publish the complete set, restoring previous files if any move fails.
  existing <- outputs[file.exists(outputs)]
  backups <- file.path(stage, paste0("previous", seq_along(existing)))
  moved <- logical(length(existing))
  published <- logical(length(outputs))
  committed <- FALSE
  on.exit(
    {
      if (!committed) {
        unlink(outputs[published])
        for (i in which(moved)) {
          if (!file.rename(backups[i], existing[i])) {
            keep_stage <- TRUE
            cli::cli_abort("Could not restore PLINK output; backup retained at {.file {backups[i]}}.",
              class = "galick_error_plink_write_failed"
            )
          }
        }
      }
    },
    add = TRUE,
    after = FALSE
  )
  for (i in seq_along(existing)) {
    moved[i] <- file.rename(existing[i], backups[i])
    if (!moved[i]) {
      cli::cli_abort("Could not stage previous PLINK output.", class = "galick_error_plink_write_failed")
    }
  }
  for (i in seq_along(outputs)) {
    published[i] <- file.rename(staged[i], outputs[i])
    if (!published[i]) {
      cli::cli_abort("Could not publish PLINK output.", class = "galick_error_plink_write_failed")
    }
  }
  committed <- TRUE
  result
}
