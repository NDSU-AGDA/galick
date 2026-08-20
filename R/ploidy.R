#' Detect chromosome-derived species ploidy
#'
#' Infers subgenome ploidy from chromosome names such as `chr1A`, `chr1B`,
#' and `chr1D`. This describes chromosome naming, not the number of allele
#' slots stored in GT. Ordinary chromosome names are reported as diploid;
#' an empty VCF returns `NA_integer_`.
#'
#' @param vcf Path to a VCF file.
#' @return A single integer chromosome-derived ploidy estimate.
#' @export
detect_chrom_ploidy <- function(vcf) {
  x <- if (inherits(vcf, "beagler_genotype_vcf")) {
    vcf
  } else {
    read_genotype_vcf(vcf)
  }
  chrom <- unique(as.character(x$fixed$CHROM))
  if (length(chrom) == 0L) return(NA_integer_)

  labels <- sub("^(?i:chr)", "", chrom, perl = TRUE)
  has_subgenome <- grepl("^[0-9]+[A-Za-z](?:[._-].*)?$", labels)
  if (!any(has_subgenome)) return(2L)

  subgenomes <- toupper(sub(
    "^[0-9]+([A-Za-z])(?:[._-].*)?$",
    "\\1",
    labels[has_subgenome]
  ))
  as.integer(2L * length(unique(subgenomes)))
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
#' @export
detect_genotype_ploidy <- function(vcf) {
  x <- if (inherits(vcf, "beagler_genotype_vcf")) {
    vcf
  } else {
    read_genotype_vcf(vcf)
  }
  calls <- as.character(x$gt)
  calls <- calls[!is.na(calls) & nzchar(calls) & calls != "."]
  if (length(calls) == 0L) return(NA_integer_)

  ploidies <- vapply(strsplit(calls, "[/|]"), length, integer(1))
  as.integer(max(ploidies))
}
