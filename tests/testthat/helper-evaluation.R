evaluation_fixture <- function() {
  samples <- paste0("s", 1:8)
  calls <- outer(seq_len(12), seq_along(samples), function(i, j) (i + j) %% 3)
  gt <- matrix(c("0/0", "0/1", "1/1")[calls + 1], nrow(calls), ncol(calls))
  records <- vapply(seq_len(nrow(gt)), function(i) {
    paste(c("1", i * 100, paste0("rs", i), "A", "G", ".", "PASS", ".", "GT", gt[i, ]), collapse = "\t")
  }, character(1))
  path <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  writeLines(c(
    "##fileformat=VCFv4.2", "##contig=<ID=1>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"),
    records
  ), path)
  metadata <- data.frame(
    sample_id = samples, partition = rep(c("tune", "test"), each = 4),
    group = rep(c("A", "B"), 4), block = paste0("family", rep(1:4, each = 2))
  )
  list(
    path = path, samples = metadata, truth = read_genotype_vcf(path),
    candidates = list(base = list(extra_args = "iterations=1"), good = list(extra_args = "iterations=2"))
  )
}

evaluation_reference <- function(fixture) {
  x <- fixture$truth
  x$samples <- paste0("r", seq_along(x$samples))
  colnames(x$gt) <- x$samples
  x$gt[, ] <- gsub("/", "|", x$gt, fixed = TRUE)
  path <- withr::local_tempfile(fileext = ".vcf", .local_envir = parent.frame())
  write_genotype_vcf(x, path)
  path
}

evaluation_mock_runner <- function(truth, mode = "normal", observer = NULL) {
  runner <- function(gt, out, ref, map, mode, nthreads, memory_gb, seed, beagle_jar, java, extra_args, run_dir) {
    input <- read_genotype_vcf(gt)
    partition <- if (grepl("_test_", run_dir, fixed = TRUE)) "test" else "tune"
    good <- "iterations=2" %in% extra_args
    if (!is.null(observer)) observer(input, seed, run_dir)
    x <- truth
    x$fixed <- data.table::copy(truth$fixed)
    x$samples <- input$samples
    x$gt <- truth$gt[, match(input$samples, truth$samples), drop = FALSE]
    x$format_fields <- list()
    x$formats <- rep("GT", nrow(x$gt))
    wrong <- !good
    if (behavior == "test_reversal" && partition == "test") wrong <- good
    if (wrong) {
      values <- genotype_counts(x)
      x$gt[, ] <- c("0/1", "1/1", "0/0")[values + 1]
    }
    x$format_fields$DS <- matrix(vapply(as.character(x$gt), gt_to_dosage, numeric(1)), nrow(x$gt), ncol(x$gt))
    x$format_fields$DS[, ] <- as.character(x$format_fields$DS)
    x$formats <- rep("GT:DS", nrow(x$gt))
    x$meta <- c(x$meta, '##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">')
    if (behavior == "missing" && good) x$gt[, ] <- "./."
    if (behavior == "missing" && good) x$format_fields$DS[, ] <- "."
    if (behavior == "absent_ds") {
      x$format_fields <- list()
      x$formats <- rep("GT", nrow(x$gt))
    }
    if (behavior == "swap" && good) {
      x$fixed$REF <- "G"
      x$fixed$ALT <- "A"
    }
    if (behavior == "drop" && good) {
      x$fixed <- x$fixed[-1]
      x$gt <- x$gt[-1, , drop = FALSE]
      x$format_fields$DS <- x$format_fields$DS[-1, , drop = FALSE]
      x$formats <- x$formats[-1]
    }
    if (behavior == "permuted") {
      index <- rev(seq_along(x$samples))
      x$samples <- x$samples[index]
      x$gt <- x$gt[, index, drop = FALSE]
      x$format_fields <- lapply(x$format_fields, function(value) value[, index, drop = FALSE])
    }
    output <- paste0(out, ".vcf")
    write_genotype_vcf(x, output)
    record <- write_run_record_json(run_dir,
      command = c("java", "test.jar"),
      versions = list(galick = "0.1.0", r = "4.1.0", os = "test"), runtime_sec = 1,
      status = if (behavior == "failed" && good) 1L else 0L,
      started_at = "2026-10-21T00:00:00Z", ended_at = "2026-10-21T00:00:01Z"
    )
    new_galick_run(
      status = if (behavior == "failed" && good) 1L else 0L,
      runtime_sec = 1, output_paths = output, run_record = record
    )
  }
  behavior <- mode
  runner
}
