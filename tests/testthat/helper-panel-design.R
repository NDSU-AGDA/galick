panel_fixture <- function() {
  samples <- paste0("p", seq_len(120))
  j <- seq_along(samples) - 1L
  values <- rbind(j %% 3, (j %/% 3) %% 3, (j %/% 9) %% 3, j %% 3, j %% 3, j %% 3)
  values <- rbind(values, values)
  gt <- matrix(c("0/0", "0/1", "1/1")[values + 1L], nrow(values))
  fixed <- data.table::data.table(
    CHROM = rep(c("1", "2"), each = 6L),
    POS = rep(seq(100, 600, 100), 2L), ID = paste0("m", seq_len(12)), REF = "A", ALT = "G",
    QUAL = ".", FILTER = ".", INFO = "."
  )
  truth <- new_genotype_vcf(
    fixed, gt, samples,
    c("##fileformat=VCFv4.2", "##contig=<ID=1>", "##contig=<ID=2>")
  )
  directory <- withr::local_tempdir(.local_envir = parent.frame())
  path <- file.path(directory, "truth.vcf")
  write_genotype_vcf(truth, path)
  ref <- truth
  ref$samples <- paste0("ref", seq_along(samples))
  ref$gt[, ] <- gsub("/", "|", gt, fixed = TRUE)
  colnames(ref$gt) <- ref$samples
  reference <- file.path(directory, "reference.vcf")
  write_genotype_vcf(ref, reference)
  metadata <- data.frame(
    sample_id = samples, partition = rep(c("design", "tune", "test"), each = 40),
    group = rep(rep(c("A", "B"), each = 20), 3), block = samples
  )
  input <- list(truth = truth)
  list(
    path = path, reference = reference, truth = truth, samples = metadata,
    candidates = panel_markers(input, c(1:4, 7:10)), scoring = panel_markers(input, c(5:6, 11:12))
  )
}

panel_arguments <- function(f) {
  list(
    truth_vcf = f$path, reference_vcf = f$reference, sample_data = f$samples,
    candidate_markers = f$candidates, scoring_markers = f$scoring, n_markers = 2L,
    output_dir = withr::local_tempdir(.local_envir = parent.frame()),
    field = "DS", repeats = 1L, max_proposals = 2L, max_runs = 8L
  )
}

panel_mock_runner <- function(truth, behavior = "normal", observer = NULL) {
  force(truth)
  force(behavior)
  function(gt, out, ref, map, mode, nthreads, memory_gb, seed, beagle_jar, java, extra_args, run_dir) {
    input <- read_genotype_vcf(gt)
    partition <- strsplit(basename(run_dir), "_", fixed = TRUE)[[1]][2]
    if (!is.null(observer)) observer(input, seed, run_dir)
    x <- truth
    x$samples <- input$samples
    x$gt <- truth$gt[, match(input$samples, truth$samples), drop = FALSE]
    good <- all(input$fixed$POS == 100)
    if (behavior == "tie") good <- FALSE
    if (behavior == "reversal" && partition == "test") good <- !good
    dosage <- genotype_counts(x)
    dosage[is.na(dosage)] <- 0
    if (!good) dosage <- (dosage + 1) %% 3
    x$gt[, ] <- c("0/0", "0/1", "1/1")[dosage + 1]
    x$format_fields <- list(DS = matrix(as.character(dosage), nrow(dosage)))
    if (behavior == "group_harm") {
      adjusted <- genotype_counts(truth)[, match(input$samples, truth$samples), drop = FALSE]
      affected <- if (good) seq_len(ncol(adjusted)) > 20 else seq_len(ncol(adjusted)) <= 20
      adjusted[, affected] <- if (good) {
        ifelse(adjusted[, affected] < 2, adjusted[, affected] + 0.2, adjusted[, affected] - 0.2)
      } else {
        (adjusted[, affected] + 1) %% 3
      }
      x$format_fields$DS <- matrix(as.character(adjusted), nrow(adjusted))
    }
    x$formats <- rep("GT:DS", nrow(x$gt))
    x$meta <- c(x$meta, '##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">')
    if (behavior == "missing" && good) x$format_fields$DS[1:6, ] <- "."
    if (behavior == "missing_baseline" && !good) x$format_fields$DS[, ] <- "."
    if (behavior == "absent_ds") {
      x$format_fields <- list()
      x$formats[] <- "GT"
    }
    output <- paste0(out, ".vcf")
    write_genotype_vcf(x, output)
    new_beagler_run(
      status = if (behavior == "failed") 1L else 0L,
      runtime_sec = 1, output_paths = output
    )
  }
}
