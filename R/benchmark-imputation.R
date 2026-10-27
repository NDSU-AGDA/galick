#' Benchmark imputation by mask -> impute -> compare
#'
#' Convenience orchestrator: treats `vcf` as the truth, masks a fraction
#' of its non-missing GT cells with [mask_genotypes()], runs
#' [run_beagle()] on the masked file (phasing-only by default), and
#' compares the imputed output to the truth via
#' [compare_imputed_vs_true()].
#'
#' Requires Java + Beagle (resolved internally). The seed controls both
#' masking and Beagle. Random-cell masking does not establish accuracy for
#' population, family or array transfer. Keep truth samples out of references.
#'
#' @param vcf Path to the input VCF (treated as the truth).
#' @param output_dir Directory for the masked + imputed files.
#' @param fraction Fraction of non-missing GT cells to mask.
#' @param seed Optional integer seed for reproducibility.
#' @param mode Beagle mode passed to [run_beagle()] (`"impute"` or
#'   `"phase"`).
#' @param field Prediction field to evaluate: GT counts or expected DS dosage.
#' @param ... Extra args forwarded to [run_beagle()] (e.g. `nthreads`,
#'   `memory_gb`, `ref`, `map`).
#'
#' @return A list with three elements: `mask` (a `galick_result`),
#'   `run` (a `galick_run`), and `comparison` (a `galick_result`).
#' @export
benchmark_imputation <- function(vcf,
                                 output_dir,
                                 fraction = 0.1,
                                 seed = NULL,
                                 mode = c("impute", "phase"),
                                 field = c("GT", "DS"),
                                 ...) {
  mode <- match.arg(mode)
  field <- match.arg(field)
  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf) || !file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a path to an existing VCF file.",
      class = "galick_error_benchmark_imputation_invalid_input"
    )
  }
  if (length(output_dir) != 1L || !is.character(output_dir) || is.na(output_dir)) {
    cli::cli_abort(
      "{.arg output_dir} must be a single non-NA character string.",
      class = "galick_error_benchmark_imputation_invalid_input"
    )
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  masked_vcf <- file.path(output_dir, "masked.vcf")
  imputed_prefix <- file.path(output_dir, "imputed")
  imputed_vcf <- paste0(imputed_prefix, ".vcf.gz")

  mask_result <- mask_genotypes(
    vcf = vcf,
    output_vcf = masked_vcf,
    fraction = fraction,
    seed = seed,
    overwrite = TRUE
  )

  run_result <- run_beagle(
    gt = masked_vcf,
    out = imputed_prefix,
    mode = mode,
    seed = seed,
    overwrite = TRUE,
    ...
  )

  comparison <- if (run_result$status == 0L && file.exists(imputed_vcf)) {
    compare_imputed_vs_true(
      imputed_vcf = imputed_vcf,
      true_vcf = vcf,
      masked_vcf = masked_vcf,
      field = field
    )
  } else {
    new_galick_result(
      errors = issues_tibble(list(list(
        code = "benchmark_imputation_run_failed",
        message = sprintf(
          "Beagle run did not produce expected output (status = %s).",
          run_result$status
        ),
        where = list()
      ))),
      warnings = issues_tibble(list()),
      summary = empty_compare_summary(),
      call = match.call()
    )
  }

  list(
    mask = mask_result,
    run = run_result,
    comparison = comparison
  )
}
