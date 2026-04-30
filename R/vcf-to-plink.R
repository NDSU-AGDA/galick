#' Convert a VCF to PLINK BED/BIM/FAM format
#'
#' Shells out to PLINK 1.9 or 2.0 with `--vcf <vcf> --make-bed --out <prefix>`.
#' Returns a `beagler_run` object capturing the command, status, runtime,
#' and output paths.
#'
#' @param vcf Path to a VCF (`.vcf` or `.vcf.gz`).
#' @param plink_prefix Output PLINK prefix. PLINK writes
#'   `<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam`, plus a log file.
#' @param plink Optional path to `plink`; auto-discovered if `NULL`.
#' @param overwrite Refuse to clobber any existing
#'   `<plink_prefix>.{bed,bim,fam}` unless `TRUE`.
#'
#' @return A `beagler_run` S3 object.
#' @export
vcf_to_plink <- function(vcf, plink_prefix, plink = NULL, overwrite = FALSE) {
  the_call <- match.call()

  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf) || !file.exists(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a path to an existing VCF file.",
      class = "beagler_error_vcf_to_plink_invalid_input"
    )
  }
  if (!grepl("\\.vcf$|\\.vcf\\.gz$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_vcf_to_plink_invalid_input"
    )
  }
  if (length(plink_prefix) != 1L || !is.character(plink_prefix) || is.na(plink_prefix)) {
    cli::cli_abort(
      "{.arg plink_prefix} must be a single non-NA character string.",
      class = "beagler_error_vcf_to_plink_invalid_input"
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
      class = "beagler_error_vcf_to_plink_output_exists"
    )
  }

  plink_path <- find_plink(plink)
  args <- c("--vcf", vcf, "--make-bed", "--out", plink_prefix)
  result <- run_cli(plink_path, args = args)

  output_paths <- expected[file.exists(expected)]
  versions <- list(
    beagler = as.character(utils::packageVersion("beagler")),
    r = as.character(getRversion()),
    os = R.version$os,
    plink = "PLINK"
  )

  new_beagler_run(
    status = result$status,
    runtime_sec = result$runtime_sec,
    command = c(plink_path, args),
    version = versions,
    output_paths = output_paths,
    stdout_log = NA_character_,
    stderr_log = NA_character_,
    run_record = NA_character_,
    call = the_call
  )
}
