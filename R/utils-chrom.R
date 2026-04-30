#' Chromosome label utilities (internal)
#' @name utils-chrom
#' @keywords internal
NULL

#' Detect the style of a vector of chromosome labels
#'
#' Returns one of `"chr_prefixed"`, `"numeric"`, `"letter"`, `"mixed"`, or
#' `"unknown"`. `"mixed"` is reserved for cases where every label matches one
#' of the recognized styles, but the styles span more than one. `"unknown"`
#' covers anything else, including wheat-style labels like `"1A"`.
#'
#' @keywords internal
#' @noRd
chrom_style <- function(x) {
  s <- unique(as.character(x))
  s <- s[!is.na(s) & nzchar(s)]
  if (length(s) == 0L) return("unknown")

  is_chr <- grepl("^chr[0-9A-Za-z]+$", s, ignore.case = TRUE)
  is_num <- grepl("^[0-9]+$", s) & !is_chr
  is_let <- grepl("^[A-Za-z]+$", s) & !is_chr

  recognized <- is_chr | is_num | is_let
  if (!all(recognized)) return("unknown")

  if (all(is_chr)) return("chr_prefixed")
  if (all(is_num)) return("numeric")
  if (all(is_let)) return("letter")

  "mixed"
}

#' Alias for `chrom_style()` to avoid shadowing in callers that take a
#' `chrom_style` parameter (e.g. `validate_marker_map()`).
#' @keywords internal
#' @noRd
detect_chrom_style <- function(x) chrom_style(x)

#' Normalize chromosome labels toward a target style
#'
#' `target = "numeric"` strips any leading `chr`/`Chr`/`CHR`. `target = "chr_prefixed"`
#' strips and re-adds a lowercase `chr` prefix. `target = "as_is"` returns input
#' unchanged (after `as.character()`).
#'
#' @keywords internal
#' @noRd
normalize_chrom_label <- function(x, target = c("numeric", "chr_prefixed", "as_is")) {
  target <- match.arg(target)
  s <- as.character(x)
  if (target == "as_is") return(s)

  na_mask <- is.na(s)
  stripped <- sub("^[Cc][Hh][Rr]", "", s)
  out <- if (target == "numeric") stripped else paste0("chr", stripped)
  out[na_mask] <- NA_character_
  out
}
