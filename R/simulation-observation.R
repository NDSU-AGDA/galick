#' @keywords internal
#' @noRd
simulation_noise_settings <- function(missingness, error_rate, panel, founders) {
  allowed <- c("rate", "sample_sd", "marker_sd")
  bad_names <- is.null(names(missingness)) || anyDuplicated(names(missingness)) ||
    any(!names(missingness) %in% allowed)
  if (!is.list(missingness) || (length(missingness) && bad_names)) {
    simulation_abort("missingness accepts rate, sample_sd and marker_sd only.")
  }
  settings <- utils::modifyList(list(rate = 0.05, sample_sd = 0, marker_sd = 0), missingness)
  simulation_number(settings$rate, "missingness rate", 0, 1)
  simulation_number(settings$sample_sd, "sample_sd", 0, 10)
  simulation_number(settings$marker_sd, "marker_sd", 0, 10)
  simulation_number(error_rate, "error_rate", 0, 1)
  keep <- rep(TRUE, nrow(founders$vcf$gt))
  if (!is.null(panel)) {
    fields <- c("chrom", "pos", "ref", "alt")
    if (!is.data.frame(panel) || !all(fields %in% names(panel)) || !nrow(panel) || anyNA(panel[fields]) ||
          !is.numeric(panel$pos) || any(!is.finite(panel$pos) | panel$pos < 1 | panel$pos != floor(panel$pos))) {
      simulation_abort("observed_markers requires valid chrom, pos, ref and alt columns.")
    }
    keys <- paste(panel$chrom, format_vcf_pos(panel$pos), panel$ref, panel$alt, sep = "\r")
    if (anyDuplicated(keys) || any(!keys %in% vcf_variant_keys(founders$vcf))) {
      simulation_abort("Observed marker identities must match founders exactly once.")
    }
    keep <- vcf_variant_keys(founders$vcf) %in% keys
  }
  list(missingness = settings, error_rate = error_rate, keep = keep)
}

#' @keywords internal
#' @noRd
simulation_observe <- function(truth, settings, ploidy) {
  x <- truth
  x$fixed <- data.table::copy(truth$fixed[which(settings$keep)])
  x$gt <- truth$gt[settings$keep, , drop = FALSE]
  original <- x$gt
  x$gt[, ] <- gsub("|", "/", x$gt, fixed = TRUE)
  x$format_fields <- list()
  x$formats <- rep("GT", nrow(x$gt))
  x$meta <- x$meta[!grepl("^##FORMAT=", x$meta)]
  m <- settings$missingness
  sample_effect <- stats::rnorm(ncol(x$gt), sd = m$sample_sd)
  marker_effect <- stats::rnorm(nrow(x$gt), sd = m$marker_sd)
  effects <- outer(marker_effect, sample_effect, "+")
  if (m$rate %in% c(0, 1)) {
    probability <- matrix(m$rate, nrow(x$gt), ncol(x$gt))
  } else {
    mid <- stats::qlogis(m$rate)
    span <- max(abs(effects)) + 40
    intercept <- stats::uniroot(function(z) mean(stats::plogis(z + effects)) - m$rate,
      c(mid - span, mid + span),
      tol = 1e-10
    )$root
    probability <- stats::plogis(intercept + effects)
  }
  missing <- matrix(stats::runif(length(x$gt)) < probability, nrow(x$gt), ncol(x$gt))
  errors <- matrix(stats::runif(length(x$gt)) < settings$error_rate, nrow(x$gt), ncol(x$gt)) & !missing
  for (index in which(errors)) {
    alleles <- strsplit(x$gt[index], "/", fixed = TRUE)[[1]]
    copy <- sample.int(ploidy, 1L)
    alleles[copy] <- as.character(1L - as.integer(alleles[copy]))
    x$gt[index] <- paste(alleles, collapse = "/")
  }
  x$gt[missing] <- paste(rep(".", ploidy), collapse = "/")
  changed <- which(missing | errors, arr.ind = TRUE)
  events <- data.frame(
    sample_id = x$samples[changed[, 2]], chrom = x$fixed$CHROM[changed[, 1]],
    pos = x$fixed$POS[changed[, 1]], event = ifelse(missing[changed], "missing", "error"),
    truth = original[changed], observed = x$gt[changed], stringsAsFactors = FALSE
  )
  list(
    vcf = x, events = events, summary = list(
      n_panel_cells = length(x$gt),
      n_off_panel_cells = sum(!settings$keep) * ncol(x$gt), n_missing = sum(missing),
      expected_missing_rate = mean(probability), realized_missing_rate = mean(missing),
      n_called = sum(!missing), n_errors = sum(errors),
      realized_error_rate = if (any(!missing)) sum(errors) / sum(!missing) else NA_real_
    ),
    by_sample = data.frame(
      sample_id = x$samples, expected_missing_rate = colMeans(probability),
      realized_missing_rate = colMeans(missing)
    ),
    by_marker = data.frame(
      chrom = x$fixed$CHROM, pos = x$fixed$POS,
      expected_missing_rate = rowMeans(probability), realized_missing_rate = rowMeans(missing)
    )
  )
}

#' @keywords internal
#' @noRd
simulation_diagnostics <- function(vcf, metadata, role, ploidy) {
  counts <- genotype_counts(vcf)
  by_marker <- list()
  for (population in unique(metadata$population)) {
    columns <- match(metadata$sample_id[metadata$population == population], vcf$samples)
    dosage <- counts[, columns, drop = FALSE]
    called <- !is.na(dosage)
    n_called <- rowSums(called)
    frequency <- rowSums(dosage, na.rm = TRUE) / (ploidy * n_called)
    heterozygosity <- rowSums(dosage > 0 & dosage < ploidy, na.rm = TRUE) / n_called
    frequency[n_called == 0] <- NA_real_
    heterozygosity[n_called == 0] <- NA_real_
    by_marker[[population]] <- data.frame(
      role = role, population = population,
      chrom = vcf$fixed$CHROM, pos = vcf$fixed$POS, n_called = n_called,
      allele_frequency = frequency, heterozygosity = heterozygosity
    )
  }
  list(by_marker = data.table::rbindlist(by_marker), counts = counts)
}

#' @keywords internal
#' @noRd
simulation_ld <- function(vcf, metadata, map, role) {
  counts <- genotype_counts(vcf)
  output <- list()
  edges <- c(0, 1, 5, 10, 50, Inf)
  for (population in unique(metadata$population)) {
    columns <- match(metadata$sample_id[metadata$population == population], vcf$samples)
    pairs <- list()
    for (chrom in unique(vcf$fixed$CHROM)) {
      rows <- which(vcf$fixed$CHROM == chrom)
      for (lag in seq_len(min(10L, length(rows) - 1L))) {
        first <- rows[seq_len(length(rows) - lag)]
        second <- rows[seq_len(length(rows) - lag) + lag]
        r2 <- vapply(seq_along(first), function(i) {
          a <- counts[first[i], columns]
          b <- counts[second[i], columns]
          use <- !is.na(a) & !is.na(b)
          comparison_r2(a[use], b[use])
        }, numeric(1))
        pairs[[length(pairs) + 1L]] <- data.frame(distance_cm = map$cm[second] - map$cm[first], r2 = r2)
      }
    }
    pairs <- data.table::rbindlist(pairs)
    for (i in seq_len(length(edges) - 1L)) {
      use <- if (nrow(pairs)) pairs$distance_cm >= edges[i] & pairs$distance_cm < edges[i + 1L] else logical()
      values <- if (nrow(pairs)) pairs$r2[use] else numeric()
      defined <- is.finite(values)
      output[[length(output) + 1L]] <- data.frame(
        role = role, population = population,
        lower_cm = edges[i], upper_cm = edges[i + 1L], n_pairs = length(values), n_defined = sum(defined),
        mean_dosage_r2 = if (any(defined)) mean(values[defined]) else NA_real_
      )
    }
  }
  data.table::rbindlist(output)
}
