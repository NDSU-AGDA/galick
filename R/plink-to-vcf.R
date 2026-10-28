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
#' @return A `galick_run` S3 object.
#' @examplesIf nzchar(Sys.which("plink"))
#' vcf <- system.file("extdata", "tiny_founders.vcf", package = "galick")
#' directory <- tempfile("plink_roundtrip_")
#' dir.create(directory)
#' prefix <- file.path(directory, "calls")
#' converted <- vcf_to_plink(vcf, prefix)
#' if (identical(converted$status, 0L)) {
#'   restored <- plink_to_vcf(prefix, file.path(directory, "restored.vcf"))
#'   restored$status
#' }
#' unlink(directory, recursive = TRUE)
#' @export
plink_to_vcf <- function(plink_prefix, vcf, plink = NULL, overwrite = FALSE) {
  the_call <- match.call()

  if (length(plink_prefix) != 1L || !is.character(plink_prefix) || is.na(plink_prefix)) {
    cli::cli_abort(
      "{.arg plink_prefix} must be a single non-NA character string.",
      class = "galick_error_plink_to_vcf_invalid_input"
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
      class = "galick_error_plink_to_vcf_invalid_input"
    )
  }
  if (length(vcf) != 1L || !is.character(vcf) || is.na(vcf)) {
    cli::cli_abort(
      "{.arg vcf} must be a single non-NA character string.",
      class = "galick_error_plink_to_vcf_invalid_input"
    )
  }
  if (!grepl("\\.vcf$", vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg vcf} must end in {.val .vcf} (PLINK writes plain .vcf).",
      class = "galick_error_plink_to_vcf_invalid_input"
    )
  }
  if (file.exists(vcf) && !overwrite) {
    cli::cli_abort(
      "Output VCF already exists: {.file {vcf}} (use {.code overwrite = TRUE}).",
      class = "galick_error_plink_to_vcf_output_exists"
    )
  }

  plink_path <- find_plink(plink)
  args <- c(
    "--bfile", normalize_path_safe(plink_prefix, FALSE), "--recode", "vcf"
  )
  result <- run_plink_bridge(plink_path, args, vcf)

  output_paths <- if (result$status == 0L) vcf else character()

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
