#' Read lossless HapMap allele calls into the genotype representation
#'
#' This adapter reads HapMap calls emitted by `export_genotypes()`, where each
#' allele is written explicitly and `/` or `|` retains the VCF phase separator.
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
  if (is.na(call) || !nzchar(call) || call == ".") return(".")
  slots <- strsplit(call, "[/|]")[[1]]
  separators <- regmatches(call, gregexpr("[/|]", call))[[1]]
  alleles <- vcf_alleles(ref, alt)
  indices <- vapply(slots, function(slot) {
    if (slot == ".") return(".")
    index <- match(slot, alleles) - 1L
    if (is.na(index)) {
      cli::cli_abort(
        "HapMap call allele {.val {slot}} is not in the marker allele list.",
        class = "beagler_error_hapmap_unresolvable_allele"
      )
    }
    as.character(index)
  }, character(1))
  interleave_gt(indices, separators)
}
