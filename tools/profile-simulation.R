pkgload::load_all(".")
directory <- tempfile("galick_profile_")
dir.create(directory)
withr::defer(unlink(directory, recursive = TRUE), envir = globalenv())
fixed <- data.table::data.table(
  CHROM = rep(c("1", "2"), each = 500),
  POS = rep(seq(100, 50000, 100), 2), ID = paste0("m", 1:1000), REF = "A", ALT = "G",
  QUAL = ".", FILTER = "PASS", INFO = "."
)
calls <- outer(seq_len(1000), seq_len(20), function(i, j) c("0|0", "0|1", "1|0", "1|1")[(i %/% 20 + j) %% 4 + 1])
founders <- galick:::new_genotype_vcf(fixed, calls, paste0("f", 1:20), "##fileformat=VCFv4.2")
path <- file.path(directory, "founders.vcf")
galick:::write_genotype_vcf(founders, path)
map <- data.frame(chrom = fixed$CHROM, pos = fixed$POS, cm = rep(seq(0, 100, length.out = 500), 2))
timing <- system.time(result <- simulate_vcf(path, map, file.path(directory, "simulation"),
  n_target = 100L, n_reference = 50L, population_size = 200L, generations = 3L, seed = 142L
))
print(timing)
cat("Estimated allele cells:", result$model$estimated_allele_cells, "\n")
cat("Output bytes:", sum(file.info(list.files(result$summary$output_dir, full.names = TRUE))$size), "\n")
cat("Pedigree rows:", nrow(data.table::fread(result$paths[["pedigree"]])), "\n")
