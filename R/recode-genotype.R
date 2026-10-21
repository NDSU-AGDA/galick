#' Recode genotype tokens to VCF GT strings
#'
#' Internal helper. Takes a vector of tokens plus its encoding name and
#' returns a character vector of VCF GT strings (`"0/0"`, `"0/1"`, `"1/1"`,
#' or `"./."` for missing).
#'
#' For non-IUPAC encodings, `ref` and `alt` are ignored. For IUPAC, both are
#' required and must be distinct biallelic SNP alleles from `{A, C, G, T}`.
#'
#' Out-of-alphabet tokens raise `beagler_error_invalid_genotype_token` so the
#' caller (typically `csv_to_vcf()`) can wrap them as per-marker errors.
#'
#' @param x A vector of tokens (any atomic type; coerced to character).
#' @param encoding One of `"AHB"`, `"012"`, `"centered"`, `"AABB"`, `"VCF"`,
#'   `"REFHET"`, `"IUPAC"`.
#' @param ref,alt Reference and alternate alleles for IUPAC decoding.
#'   Ignored for other encodings.
#' @param ahb_ref Either `"A"` or `"B"`. With `"A"` (default), `A` codes
#'   homo-REF and `B` codes homo-ALT. Pass `"B"` to swap.
#' @param aabb_ref_homo Either `"AA"` or `"BB"`. Mirror of `ahb_ref` for
#'   AABB encoding.
#' @param missing_codes Tokens treated as missing.
#'
#' @return Character vector of `"0/0"` / `"0/1"` / `"1/1"` / `"./."`.
#' @keywords internal
#' @noRd
recode_to_vcf_gt <- function(x,
                             encoding,
                             ref = NULL,
                             alt = NULL,
                             ahb_ref = c("A", "B"),
                             aabb_ref_homo = c("AA", "BB"),
                             missing_codes = c(NA, ".", "--", "", "N", "NN", "./.")) {
  ahb_ref <- match.arg(ahb_ref)
  aabb_ref_homo <- match.arg(aabb_ref_homo)

  s <- toupper(trimws(as.character(x)))
  miss_norm <- toupper(trimws(as.character(missing_codes)))
  is_missing <- is.na(s) | !nzchar(s) | s %in% miss_norm

  out <- character(length(s))
  if (length(out) > 0L) out[is_missing] <- "./."

  if (any(!is_missing)) {
    to_recode <- s[!is_missing]
    out[!is_missing] <- switch(encoding,
      AHB = recode_ahb(to_recode, ahb_ref),
      `012` = recode_012(to_recode),
      centered = recode_centered(to_recode),
      AABB = recode_aabb(to_recode, aabb_ref_homo),
      VCF = recode_vcf(to_recode),
      REFHET = recode_refhet(to_recode),
      IUPAC = recode_iupac(to_recode, ref = ref, alt = alt),
      cli::cli_abort(
        "Unknown encoding: {.val {encoding}}",
        class = "beagler_error_invalid_genotype_token"
      )
    )
  }

  out
}

#' @keywords internal
#' @noRd
recode_with_table <- function(x, table, label) {
  unrecognized <- setdiff(unique(x), names(table))
  if (length(unrecognized) > 0L) {
    cli::cli_abort(
      "Unrecognized {label} token(s): {.val {unrecognized}}",
      class = "beagler_error_invalid_genotype_token"
    )
  }
  unname(table[x])
}

#' @keywords internal
#' @noRd
recode_ahb <- function(x, ahb_ref) {
  table <- if (ahb_ref == "A") {
    c(A = "0/0", H = "0/1", B = "1/1")
  } else {
    c(A = "1/1", H = "0/1", B = "0/0")
  }
  recode_with_table(x, table, "AHB")
}

#' @keywords internal
#' @noRd
recode_012 <- function(x) {
  recode_with_table(x, c("0" = "0/0", "1" = "0/1", "2" = "1/1"), "012")
}

#' @keywords internal
#' @noRd
recode_centered <- function(x) {
  recode_with_table(x, c("-1" = "0/0", "0" = "0/1", "1" = "1/1"), "centered")
}

#' @keywords internal
#' @noRd
recode_aabb <- function(x, aabb_ref_homo) {
  table <- if (aabb_ref_homo == "AA") {
    c(AA = "0/0", AB = "0/1", BB = "1/1")
  } else {
    c(AA = "1/1", AB = "0/1", BB = "0/0")
  }
  recode_with_table(x, table, "AABB")
}

#' @keywords internal
#' @noRd
recode_vcf <- function(x) {
  recode_with_table(
    x,
    c("0/0" = "0/0", "0/1" = "0/1", "1/0" = "0/1", "1/1" = "1/1"),
    "VCF"
  )
}

#' @keywords internal
#' @noRd
recode_refhet <- function(x) {
  recode_with_table(x, c(REF = "0/0", HET = "0/1", ALT = "1/1"), "REFHET")
}

#' IUPAC decoding requires per-marker REF / ALT
#' @keywords internal
#' @noRd
recode_iupac <- function(x, ref, alt) {
  if (length(ref) != 1L || length(alt) != 1L || is.na(ref) || is.na(alt)) {
    cli::cli_abort(
      "IUPAC encoding requires scalar non-NA {.arg ref} and {.arg alt}.",
      class = "beagler_error_invalid_genotype_token"
    )
  }
  ref <- toupper(ref)
  alt <- toupper(alt)
  if (!ref %in% c("A", "C", "G", "T") ||
        !alt %in% c("A", "C", "G", "T") ||
        ref == alt) {
    cli::cli_abort(
      "IUPAC needs distinct biallelic SNP {.arg ref}/{.arg alt} from A/C/G/T.",
      class = "beagler_error_invalid_genotype_token"
    )
  }

  pair_to_iupac <- c(
    AC = "M", CA = "M",
    AG = "R", GA = "R",
    AT = "W", TA = "W",
    CG = "S", GC = "S",
    CT = "Y", TC = "Y",
    GT = "K", TG = "K"
  )
  het_code <- unname(pair_to_iupac[paste0(ref, alt)])

  table_map <- c(
    stats::setNames("0/0", ref),
    stats::setNames("1/1", alt),
    stats::setNames("0/1", het_code)
  )

  recode_with_table(x, table_map, sprintf("IUPAC (ref=%s, alt=%s)", ref, alt))
}
