#' Convert a VCF to PLINK BED/BIM/FAM format
#'
#' Shells out to PLINK 1.9 or 2.0 with `--vcf <vcf> --make-bed --out <prefix>`.
#' Returns a `galick_run` object capturing the command, status, runtime,
#' and output paths.
#'
#' @param vcf Path to a VCF (`.vcf` or `.vcf.gz`).
#' @param plink_prefix Output PLINK prefix. PLINK writes
#'   `<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam`, plus a log file.
#' @param plink Optional path to `plink`; auto-discovered if `NULL`.
#' @param overwrite Refuse to clobber any existing
#'   `<plink_prefix>.{bed,bim,fam}` unless `TRUE`.
#'
#' @return A `galick_run` S3 object.
#' @examplesIf identical(Sys.getenv("GALICK_RUN_EXAMPLES"), "true")
#' vcf <- system.file("extdata", "tiny_founders.vcf", package = "galick")
#' directory <- tempfile("plink_")
#' dir.create(directory)
#' result <- vcf_to_plink(vcf, file.path(directory, "calls"))
#' result$status
#' unlink(directory, recursive = TRUE)
#' @export
vcf_to_plink <- function(vcf, plink_prefix, plink = NULL, overwrite = FALSE) {
  the_call <- match.call()

  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf) || !file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a path to an existing VCF file.",
      class = "galick_error_vcf_to_plink_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "galick_error_vcf_to_plink_invalid_input"
    )
  }
  if (length(plink_prefix) != 1L || !is.character(plink_prefix) || is.na(plink_prefix)) {
    cli::cli_abort(
      "{.arg plink_prefix} must be a single non-NA character string.",
      class = "galick_error_vcf_to_plink_invalid_input"
    )
  }

  expected <- paste0(plink_prefix, c(".bed", ".bim", ".fam"))
  existing <- expected[file.exists(expected)]
  if (length(existing) > 0L && !overwrite) {
    cli::cli_abort(
      c(
        "Output PLINK files already exist for prefix {.file {plink_prefix}}:",
        i = paste(existing, collapse = ", "),
        i = "Use {.code overwrite = TRUE} to clobber."
      ),
      class = "galick_error_vcf_to_plink_output_exists"
    )
  }

  plink_path <- find_plink(plink)
  result <- run_plink_bridge(plink_path, c("--vcf", normalize_path_safe(vcf, TRUE), "--make-bed"), expected)

  output_paths <- if (result$status == 0L) expected else character()
  versions <- list(
    galick = as.character(utils::packageVersion("galick")),
    r = as.character(getRversion()),
    os = R.version$os,
    plink = "PLINK"
  )

  new_galick_run(
    status = result$status,
    runtime_sec = result$runtime_sec,
    command = result$command,
    version = versions,
    output_paths = output_paths,
    stdout_log = NA_character_,
    stderr_log = NA_character_,
    run_record = NA_character_,
    call = the_call
  )
}
