#' Convert PLINK BED/BIM/FAM to a VCF file
#'
#' Shells out to PLINK 1.9 or 2.0 with `--bfile <prefix> --recode vcf --out
#' <out_prefix>`. PLINK writes `<out_prefix>.vcf`. The `vcf` argument is
#' the *full* desired VCF path (must end in `.vcf`); the function strips
#' the extension to derive PLINK's `--out` prefix.
#'
#' @param plink_prefix Input PLINK prefix (no extension). The function
#'   expects `<prefix>.{bed,bim,fam}` to exist.
#' @param vcf Output VCF path, ending in `.vcf`.
#' @param plink Optional path to `plink`; auto-discovered if `NULL`.
#' @param overwrite Refuse to clobber existing `vcf` unless `TRUE`.
#'
#' @return A `beagler_run` S3 object.
#' @export
plink_to_vcf <- function(plink_prefix, vcf, plink = NULL, overwrite = FALSE) {
  the_call <- match.call()

  if (length(plink_prefix) != 1L || !is.character(plink_prefix) || is.na(plink_prefix)) {
    cli::cli_abort(
      "{.arg plink_prefix} must be a single non-NA character string.",
      class = "beagler_error_plink_to_vcf_invalid_input"
    )
  }
  required <- paste0(plink_prefix, c(".bed", ".bim", ".fam"))
  missing_files <- required[!file.exists(required)]
  if (length(missing_files) > 0L) {
    cli::cli_abort(
      c(
        "Missing PLINK input files for prefix {.file {plink_prefix}}:",
        i = paste(missing_files, collapse = ", ")
      ),
      class = "beagler_error_plink_to_vcf_invalid_input"
    )
  }
  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a single non-NA character string.",
      class = "beagler_error_plink_to_vcf_invalid_input"
    )
  }
  if (!grepl("\\.vcf$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} (PLINK writes plain .vcf).",
      class = "beagler_error_plink_to_vcf_invalid_input"
    )
  }
  if (file.exists(vcf) && !overwrite) {
    cli::cli_abort(
      "Output VCF already exists: {.file {vcf}} (use {.code overwrite = TRUE}).",
      class = "beagler_error_plink_to_vcf_output_exists"
    )
  }

  plink_path <- find_plink(plink)
  out_prefix <- sub("\\.vcf$", "", vcf, ignore.case = TRUE)
  args <- c("--bfile", plink_prefix, "--recode", "vcf",
            "--out", out_prefix)
  result <- run_cli(plink_path, args = args)

  output_paths <- if (file.exists(vcf)) vcf else character()

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
