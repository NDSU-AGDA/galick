#' Genotype encoding detection (internal)
#' @name encoding-genotype
#' @keywords internal
NULL

# Encoding alphabets for v0.1.0. Tokens are the canonical (trim+upper)
# representation; the detector normalizes input before matching.
#
# IUPAC for v0.1.0 deliberately excludes H and B so it does not collide with
# AHB. M2-spec sign-off pins this; revisit if/when polyploid IUPAC ships.
#
# `012` and `centered` overlap on {0, 1}. Disambiguation rule: `-1` is unique
# to `centered`, so the detector prefers `centered` when `-1` is observed and
# `012` otherwise.
ENCODING_ALPHABETS <- list(
  AHB = c("A", "H", "B"),
  `012` = c("0", "1", "2"),
  centered = c("-1", "0", "1"),
  AABB = c("AA", "AB", "BB"),
  VCF = c("0/0", "0/1", "1/0", "1/1"),
  REFHET = c("REF", "HET", "ALT"),
  IUPAC = c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "N")
)

#' Detect the genotype encoding of a vector of tokens
#'
#' Strips whitespace, uppercases, and removes missing-code tokens before
#' matching. Returns one of the seven named encodings if exactly one fits
#' (smallest-alphabet-that-fits wins ties), `"mixed"` if every observed
#' token is recognized by some encoding but no single encoding contains
#' all of them, or `"unknown"` if any observed token is unrecognized or
#' the input is empty after stripping.
#'
#' @keywords internal
#' @noRd
detect_genotype_encoding <- function(x,
                                     missing_codes = c(NA, ".", "--", "", "N", "NN", "./.")) {
  observed <- normalize_observed_tokens(x, missing_codes)
  if (length(observed) == 0L) {
    return("unknown")
  }

  fits <- vapply(
    ENCODING_ALPHABETS,
    function(alpha) all(observed %in% alpha),
    logical(1)
  )

  if (any(fits)) {
    fitting <- names(fits)[fits]
    # 012 vs centered: -1 is unique to centered, so its presence pins
    # the encoding; otherwise default to 012.
    if ("centered" %in% fitting && "012" %in% fitting) {
      if ("-1" %in% observed) {
        fitting <- setdiff(fitting, "012")
      } else {
        fitting <- setdiff(fitting, "centered")
      }
    }
    sizes <- vapply(ENCODING_ALPHABETS, length, integer(1))
    return(fitting[which.min(sizes[fitting])])
  }

  all_recognized_tokens <- unique(unlist(ENCODING_ALPHABETS, use.names = FALSE))
  if (all(observed %in% all_recognized_tokens)) {
    return("mixed")
  }

  "unknown"
}

#' Trim, uppercase, drop missing-coded entries, return unique tokens
#' @keywords internal
#' @noRd
normalize_observed_tokens <- function(x, missing_codes) {
  s <- toupper(trimws(as.character(x)))
  miss_norm <- toupper(trimws(as.character(missing_codes)))
  is_missing <- is.na(s) | !nzchar(s) | s %in% miss_norm
  unique(s[!is_missing])
}
