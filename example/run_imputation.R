library(beagler)
genomic_file_path <- "/Users/didier.murilloflorez/Desktop/Genomic_selection_fles/2021_NDSU_AYT.vcf"
run_beagle(
 gt = genomic_file_path,
 out = "imputed_2021_NDSU_AYT",
 n_threads = 1,
 impute = TRUE,
 burnin = NULL,
 iterations = 10,
 phase_states = 280
)

# To unzipe the file use
# gunzip imputed_genotypes.vcf.gz
# or
# gunzip -c imputed_2021_NDSU_AYT.vcf.gz > imputed_2021_NDSU_AYT.vcf



