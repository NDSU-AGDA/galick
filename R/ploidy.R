#' Identify subgenome labels in chromosome names
#'
#' Extracts suffix labels from names such as `chr1A`, `chr1B`, and `chr1D`.
#' Labels describe naming only and do not establish species or locus ploidy.
#'
#' @param vcf Path to a VCF file.
#' @return A character vector of observed labels, empty when none match.
#' @examples
#' vcf <- system.file("extdata", "tiny_founders.vcf", package = "galick")
#' detect_subgenome_labels(vcf) # No subgenome suffix in chromosome "1".
#' @export
detect_subgenome_labels <- function(vcf) {
  x <- if (inherits(vcf, "galick_genotype_vcf")) {
    vcf
  } else {
    read_genotype_vcf(vcf)
  }
  chrom <- unique(as.character(x$fixed$CHROM))
  if (length(chrom) == 0L) {
    return(character())
  }

  labels <- sub("^(?i:chr)", "", chrom, perl = TRUE)
  has_subgenome <- grepl("^[0-9]+[A-Za-z](?:[._-].*)?$", labels)
  if (!any(has_subgenome)) {
    return(character())
  }

  subgenomes <- toupper(sub(
    "^[0-9]+([A-Za-z])(?:[._-].*)?$",
    "\\1",
    labels[has_subgenome]
  ))
  unique(subgenomes)
}

#' Detect ploidy encoded in VCF genotype fields
#'
#' Counts allele slots in every non-shorthand GT call and returns the maximum
#' observed ploidy. Separators and allele values do not affect the count. A
#' bare `.` is skipped because it may mean haploid missing or whole-call
#' missing; if no call establishes ploidy, `NA_integer_` is returned.
#'
#' @param vcf Path to a VCF file.
#' @return A single integer maximum observed genotype ploidy.
#' @examples
#' vcf <- system.file("extdata", "tiny_founders.vcf", package = "galick")
#' detect_genotype_ploidy(vcf)
#' @export
detect_genotype_ploidy <- function(vcf) {
  x <- if (inherits(vcf, "galick_genotype_vcf")) {
    vcf
  } else {
    read_genotype_vcf(vcf)
  }
  calls <- as.character(x$gt)
  calls <- calls[!is.na(calls) & nzchar(calls) & calls != "."]
  if (length(calls) == 0L) {
    return(NA_integer_)
  }

  ploidies <- vapply(strsplit(calls, "[/|]"), length, integer(1))
  as.integer(max(ploidies))
}
