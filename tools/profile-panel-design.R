# Development performance fixture, not evidence of biological accuracy.
pkgload::load_all(".", quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
directory <- if (length(args)) args[1] else tempfile("panel_profile_")
if (dir.exists(directory)) stop("Use a new output directory.")
dir.create(directory, recursive = TRUE)
samples <- paste0("sample", seq_len(120))
j <- seq_along(samples) - 1L
pattern <- rbind(j %% 3, (j %/% 3) %% 3, (j %/% 9) %% 3, j %% 3, j %% 3, j %% 3)
values <- pattern[rep(seq_len(6), each = 20), , drop = FALSE]
values <- rbind(values, values)
gt <- matrix(c("0/0", "0/1", "1/1")[values + 1L], nrow(values))
fixed <- data.table::data.table(
  CHROM = rep(c("1", "2"), each = 120), POS = rep(seq_len(120) * 1000, 2),
  ID = paste0("m", seq_len(240)), REF = "A", ALT = "G", QUAL = ".", FILTER = ".", INFO = "."
)
truth <- galick:::new_genotype_vcf(
  fixed, gt, samples, c("##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=2>")
)
truth_path <- file.path(directory, "truth.vcf")
galick:::write_genotype_vcf(truth, truth_path)
reference <- truth
reference$samples <- paste0("reference", seq_along(samples))
reference$gt[, ] <- gsub("/", "|", reference$gt, fixed = TRUE)
colnames(reference$gt) <- reference$samples
reference_path <- file.path(directory, "reference.vcf")
galick:::write_genotype_vcf(reference, reference_path)
metadata <- data.frame(
  sample_id = samples, partition = rep(c("design", "tune", "test"), each = 40),
  group = rep(rep(c("A", "B"), each = 20), 3), block = samples
)
identities <- data.frame(chrom = fixed$CHROM, pos = fixed$POS, ref = fixed$REF, alt = fixed$ALT)
elapsed <- system.time(result <- design_marker_panel(
  truth_path, reference_path, metadata,
  candidate_markers = identities[seq_len(240) %% 3 != 0, ],
  scoring_markers = identities[seq_len(240) %% 3 == 0, ],
  n_markers = 20L, output_dir = file.path(directory, "experiment"),
  repeats = 1L, max_proposals = 2L, max_runs = 8L, memory_gb = 1L,
  extra_args = c("iterations=2", "burnin=1")
))
stopifnot(result$ok)
profile <- list(
  elapsed_sec = unname(elapsed["elapsed"]), summary = result$summary,
  n_truth_samples = 120L, n_reference_samples = 120L, n_candidates = 160L, n_scoring = 80L,
  r_version = as.character(getRversion()), package_version = as.character(utils::packageVersion("galick"))
)
jsonlite::write_json(profile, file.path(directory, "profile.json"), auto_unbox = TRUE, pretty = TRUE)
print(profile)
