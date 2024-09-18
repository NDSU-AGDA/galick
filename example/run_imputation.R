library(beagler)
genomic_file_path <- Sys.getenv("FILE_PATH")
run_beagle(
 gt = genomic_file_path,
 out = "imputed_2021_NDSU_AYT",
 n_threads = 1,
 impute = TRUE,
 burnin = 3,
 iterations = 12,
 phase_states = 280
)

# To unzipe the file use
# gunzip imputed_genotypes.vcf.gz
# or
# gunzip -c imputed_2021_NDSU_AYT.vcf.gz > imputed_2021_NDSU_AYT.vcf

