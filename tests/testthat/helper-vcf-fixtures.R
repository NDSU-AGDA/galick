scientific_vcf <- function(samples, records, meta = character()) {
  path <- tempfile(fileext = ".vcf")
  withr::defer(unlink(path), envir = parent.frame())
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    meta,
    paste(c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
      "INFO", "FORMAT", samples
    ), collapse = "\t"),
    records
  ), path)
  path
}
