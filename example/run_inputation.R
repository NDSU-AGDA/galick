library(beagler)
genomic_file_path <- "/Users/didier.murilloflorez/Desktop/Genomic_selection_fles/2021_NDSU_AYT.vcf"
run_beagle(gt = genomic_file_path, out = "imputed_2021_NDSU_AYT", n_threads = 1, impute = TRUE)
