#' Read standard IUPAC HapMap calls into the genotype representation
#'
#' HapMap does not retain VCF phase sets or arbitrary genotype ploidy.
#'
#' @keywords internal
#' @noRd
read_genotype_hapmap <- function(path) {
  hmp <- read_hapmap(path)
  fixed <- data.table::data.table(
    CHROM = hmp$metadata$chrom,
    POS = as.numeric(hmp$metadata$pos),
    ID = hmp$metadata$marker_id,
    REF = hmp$metadata$ref,
    ALT = hmp$metadata$alt,
    QUAL = ".",
    FILTER = "PASS",
    INFO = "."
  )
  fixed <- refresh_vcf_ordinals(fixed)
  gt <- hmp$gt_matrix
  if (length(gt) > 0L) {
    for (i in seq_len(nrow(gt))) {
      gt[i, ] <- vapply(
        gt[i, ],
        hapmap_call_to_gt,
        character(1),
        ref = fixed$REF[i],
        alt = fixed$ALT[i]
      )
    }
  }
  meta <- c(
    "##fileformat=VCFv4.2",
    sprintf("##contig=<ID=%s>", unique(fixed$CHROM))
  )
  new_genotype_vcf(fixed, gt, hmp$samples, meta)
}

#' @keywords internal
#' @noRd
hapmap_call_to_gt <- function(call, ref, alt) {
  recode_to_vcf_gt(call, encoding = "IUPAC", ref = ref, alt = alt)
}
