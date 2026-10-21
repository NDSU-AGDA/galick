simulation_fixture <- function(ploidy = 2L, biparental = FALSE) {
  path <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  samples <- paste0("f", 1:4)
  map <- data.frame(
    chrom = rep(c("1", "2"), each = 6), pos = rep(seq(100, 600, 100), 2),
    cm = rep(seq(0, 50, 10), 2)
  )
  calls <- vapply(seq_len(nrow(map)), function(i) {
    hap <- if (biparental) {
      c(rep(0, ploidy), rep(1, ploidy), rep(0, ploidy), rep(1, ploidy))
    } else {
      as.integer((seq_len(ploidy * 4) + i %/% 3) %% 3 > 0)
    }
    vapply(split(hap, rep(1:4, each = ploidy)), paste, character(1), collapse = "|")
  }, character(4))
  records <- vapply(seq_len(nrow(map)), function(i) {
    paste(c(map$chrom[i], map$pos[i], paste0("rs", i), "A", "G", ".", "PASS", ".", "GT", calls[, i]),
      collapse = "\t"
    )
  }, character(1))
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=2>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples),
      collapse = "\t"
    ), records
  ), path)
  list(path = path, map = map)
}
