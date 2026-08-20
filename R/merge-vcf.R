#' Merge VCF callsets by samples or by markers
#'
#' `by_sample` combines different sample columns for the same ordered marker
#' set. `by_marker` combines different marker rows for the same sample set.
#' Both modes parse records before merging, preserve duplicate records with an
#' occurrence ordinal, reindex file-relative allele numbers, and drop REF
#' mismatches with an explicit report entry.
#'
#' @param vcfs Character vector of at least two `.vcf` or `.vcf.gz` files.
#' @param output_vcf Output path ending in `.vcf` or `.vcf.gz`.
#' @param mode `"by_sample"` or `"by_marker"`.
#' @param sample_collision For `by_sample`, one of `"error"`, `"suffix"`, or
#'   `"coalesce"`. Coalescing accepts identical calls or fills a fully missing
#'   call from the other file; conflicting called genotypes are refused.
#' @param overwrite Refuse to replace `output_vcf` unless `TRUE`.
#'
#' @return A `beagler_result` whose summary is the merge report.
#' @export
merge_vcf <- function(vcfs,
                      output_vcf,
                      mode = c("by_sample", "by_marker"),
                      sample_collision = c("error", "suffix", "coalesce"),
                      overwrite = FALSE) {
  the_call <- match.call()
  mode <- match.arg(mode)
  sample_collision <- match.arg(sample_collision)
  if (!is.character(vcfs) || length(vcfs) < 2L || anyNA(vcfs)) {
    cli::cli_abort(
      "{.arg vcfs} must contain at least two VCF paths.",
      class = "beagler_error_merge_invalid_input"
    )
  }
  if (any(!file.exists(vcfs))) {
    cli::cli_abort(
      "VCF input file(s) not found: {.file {vcfs[!file.exists(vcfs)]}}.",
      class = "beagler_error_merge_invalid_input"
    )
  }
  assert_export_output(output_vcf)
  if (!grepl("\\.vcf(?:\\.gz)?$", output_vcf, ignore.case = TRUE)) {
    cli::cli_abort(
      "{.arg output_vcf} must end in .vcf or .vcf.gz.",
      class = "beagler_error_merge_invalid_input"
    )
  }
  assert_output_available(output_vcf, overwrite)

  inputs <- lapply(vcfs, read_genotype_vcf)
  input_names <- make.unique(basename(vcfs))
  names(inputs) <- input_names
  sites_in <- stats::setNames(
    vapply(inputs, function(x) nrow(x$fixed), integer(1)),
    input_names
  )
  samples_in <- stats::setNames(
    vapply(inputs, function(x) length(x$samples), integer(1)),
    input_names
  )

  unresolved <- drop_unresolvable_sites(inputs, mode)
  inputs <- unresolved$inputs
  ref_mismatch <- drop_ref_mismatch_sites(inputs, mode)
  inputs <- ref_mismatch$inputs
  dropped <- unique(rbind(unresolved$dropped, ref_mismatch$dropped))

  combined <- if (mode == "by_sample") {
    merge_vcf_by_sample(inputs, sample_collision)
  } else {
    merge_vcf_by_marker(inputs)
  }
  combined$vcf$meta <- unique(unlist(lapply(inputs, `[[`, "meta")))
  combined$vcf$contigs <- vcf_contig_order(
    combined$vcf$meta,
    combined$vcf$fixed$CHROM
  )
  combined$vcf$source <- normalizePath(vcfs, mustWork = TRUE)

  written <- write_genotype_vcf(
    combined$vcf,
    output_vcf,
    overwrite = overwrite
  )
  warnings <- list()
  if (nrow(dropped) > 0L) {
    warnings[[1L]] <- list(
      code = "merge_sites_dropped",
      message = sprintf("%d site record(s) were dropped; see sites_dropped.",
                        nrow(dropped)),
      where = list(reasons = unique(dropped$reason))
    )
  }

  new_beagler_result(
    warnings = issues_tibble(warnings),
    summary = list(
      operation = "merge",
      mode = mode,
      input_paths = normalizePath(vcfs, mustWork = TRUE),
      output_paths = normalizePath(written, mustWork = TRUE),
      sites_in = sites_in,
      sites_out = nrow(combined$vcf$fixed),
      sites_dropped = dropped,
      samples_in = samples_in,
      samples_out = length(combined$vcf$samples),
      collisions = combined$collisions,
      genotype_ploidy = detect_genotype_ploidy(combined$vcf),
      chrom_ploidy = detect_chrom_ploidy(combined$vcf)
    ),
    call = the_call
  )
}

#' Drop calls that cannot be resolved against their file's REF/ALT
#' @keywords internal
#' @noRd
drop_unresolvable_sites <- function(inputs, mode) {
  bad_keys <- character()
  bad_by_file <- vector("list", length(inputs))
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    bad <- logical(nrow(x$fixed))
    if (ncol(x$gt) > 0L) {
      for (i in seq_len(nrow(x$fixed))) {
        bad[i] <- any(!vapply(
          x$gt[i, ], valid_gt_for_alt, logical(1), alt = x$fixed$ALT[i]
        ))
      }
    }
    bad_by_file[[j]] <- bad
    if (mode == "by_sample" && any(bad)) {
      bad_keys <- c(bad_keys, site_alignment_key(x$fixed[bad]))
    }
  }
  bad_keys <- unique(bad_keys)

  reports <- list()
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    bad <- if (mode == "by_sample") {
      site_alignment_key(x$fixed) %in% bad_keys
    } else {
      bad_by_file[[j]]
    }
    if (any(bad)) {
      reports[[length(reports) + 1L]] <- drop_rows_report(
        x$fixed[bad], "unresolvable alleles"
      )
      inputs[[j]] <- subset_genotype_vcf(x, !bad)
    }
  }
  list(inputs = inputs, dropped = bind_drop_reports(reports))
}

#' Drop REF disagreements instead of selecting a reference allele
#' @keywords internal
#' @noRd
drop_ref_mismatch_sites <- function(inputs, mode) {
  rows <- data.table::rbindlist(lapply(seq_along(inputs), function(j) {
    fixed <- inputs[[j]]$fixed
    if (nrow(fixed) == 0L) return(NULL)
    data.table::data.table(
      file = j,
      CHROM = fixed$CHROM,
      POS = fixed$POS,
      REF = fixed$REF,
      position_ordinal = fixed$position_ordinal
    )
  }))
  if (nrow(rows) == 0L) {
    return(list(inputs = inputs, dropped = empty_drop_report()))
  }

  if (mode == "by_sample") {
    grouped <- split(rows, site_alignment_key(rows))
    is_bad <- vapply(grouped, function(group) {
      length(unique(group$file)) > 1L && length(unique(group$REF)) > 1L
    }, logical(1))
    bad_keys <- names(grouped)[is_bad]
  } else {
    by_position <- split(rows, paste(rows$CHROM, rows$POS, sep = "\r"))
    bad_positions <- vapply(by_position, function(group) {
      refs <- split(group$REF, group$file)
      length(refs) > 1L && length(unique(vapply(
        refs,
        function(x) paste(sort(unique(x)), collapse = ","),
        character(1)
      ))) > 1L
    }, logical(1))
    bad_keys <- names(by_position)[bad_positions]
  }

  reports <- list()
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    key <- if (mode == "by_sample") {
      site_alignment_key(x$fixed)
    } else {
      paste(x$fixed$CHROM, x$fixed$POS, sep = "\r")
    }
    bad <- key %in% bad_keys
    if (any(bad)) {
      reports[[length(reports) + 1L]] <- drop_rows_report(
        x$fixed[bad], "REF mismatch"
      )
      inputs[[j]] <- subset_genotype_vcf(x, !bad)
    }
  }
  list(inputs = inputs, dropped = bind_drop_reports(reports))
}

#' @keywords internal
#' @noRd
merge_vcf_by_sample <- function(inputs, sample_collision) {
  keys <- lapply(inputs, function(x) site_alignment_key(x$fixed))
  reference_keys <- sort(keys[[1]])
  same_markers <- vapply(keys, function(key) {
    identical(sort(key), reference_keys)
  }, logical(1))
  if (!all(same_markers)) {
    cli::cli_abort(
      paste0(
        "{.code mode = \"by_sample\"} requires the same marker keys in ",
        "every file (CHROM, POS, and position ordinal)."
      ),
      class = "beagler_error_merge_marker_contract"
    )
  }

  aligned <- lapply(seq_along(inputs), function(j) {
    subset_genotype_vcf(inputs[[j]], match(keys[[1]], keys[[j]]))
  })
  n_sites <- nrow(aligned[[1]]$fixed)
  fixed <- data.table::copy(aligned[[1]]$fixed)
  reindexed <- vector("list", length(aligned))
  for (j in seq_along(aligned)) {
    reindexed[[j]] <- aligned[[j]]$gt
  }

  for (i in seq_len(n_sites)) {
    refs <- vapply(aligned, function(x) x$fixed$REF[i], character(1))
    if (length(unique(refs)) != 1L) {
      cli::cli_abort(
        "Internal REF reconciliation failed at merged record {i}.",
        class = "beagler_error_merge_ref_mismatch"
      )
    }
    union_alt <- unique(unlist(lapply(aligned, function(x) {
      setdiff(vcf_alleles(x$fixed$REF[i], x$fixed$ALT[i]), refs[1])
    })))
    fixed$ALT[i] <- if (length(union_alt) == 0L) {
      "."
    } else {
      paste(union_alt, collapse = ",")
    }
    for (j in seq_along(aligned)) {
      reindexed[[j]][i, ] <- vapply(
        aligned[[j]]$gt[i, ],
        reindex_gt_call,
        character(1),
        old_ref = aligned[[j]]$fixed$REF[i],
        old_alt = aligned[[j]]$fixed$ALT[i],
        new_ref = refs[1],
        new_alt = fixed$ALT[i]
      )
    }
  }

  sample_result <- combine_sample_matrices(
    reindexed,
    lapply(aligned, `[[`, "samples"),
    sample_collision
  )
  fixed <- refresh_vcf_ordinals(fixed)
  vcf <- new_genotype_vcf(
    fixed,
    sample_result$gt,
    sample_result$samples,
    aligned[[1]]$meta
  )
  list(vcf = vcf, collisions = sample_result$collisions)
}

#' @keywords internal
#' @noRd
merge_vcf_by_marker <- function(inputs) {
  reference_samples <- inputs[[1]]$samples
  same_samples <- vapply(inputs, function(x) {
    setequal(x$samples, reference_samples) &&
      length(x$samples) == length(reference_samples)
  }, logical(1))
  if (!all(same_samples)) {
    cli::cli_abort(
      paste0(
        "{.code mode = \"by_marker\"} requires the same sample names in ",
        "every file."
      ),
      class = "beagler_error_merge_sample_contract"
    )
  }

  marker_keys <- lapply(inputs, canonical_marker_keys)
  duplicated_across_files <- anyDuplicated(unlist(marker_keys)) > 0L
  if (duplicated_across_files) {
    all_keys <- unlist(marker_keys, use.names = FALSE)
    key_files <- rep(seq_along(marker_keys), lengths(marker_keys))
    duplicated_across_files <- any(vapply(
      split(key_files, all_keys),
      function(x) data.table::uniqueN(x) > 1L,
      logical(1)
    ))
  }
  if (duplicated_across_files) {
    cli::cli_abort(
      "{.code mode = \"by_marker\"} requires non-overlapping marker records.",
      class = "beagler_error_merge_marker_contract"
    )
  }

  fixed <- data.table::rbindlist(
    lapply(inputs, `[[`, "fixed"),
    use.names = TRUE
  )
  gt <- do.call(rbind, lapply(inputs, function(x) {
    x$gt[, match(reference_samples, x$samples), drop = FALSE]
  }))
  colnames(gt) <- reference_samples
  fixed <- refresh_vcf_ordinals(fixed)
  list(
    vcf = new_genotype_vcf(fixed, gt, reference_samples, inputs[[1]]$meta),
    collisions = empty_collision_report()
  )
}

#' @keywords internal
#' @noRd
combine_sample_matrices <- function(matrices, samples, policy) {
  all_samples <- unlist(samples, use.names = FALSE)
  duplicates <- unique(all_samples[duplicated(all_samples)])
  if (length(duplicates) > 0L && policy == "error") {
    cli::cli_abort(
      "Duplicate sample name(s): {.val {duplicates}}.",
      class = "beagler_error_merge_sample_collision"
    )
  }

  if (policy %in% c("error", "suffix")) {
    output_samples <- if (policy == "suffix") {
      make.unique(all_samples, sep = ".")
    } else {
      all_samples
    }
    gt <- do.call(cbind, matrices)
    colnames(gt) <- output_samples
    changed <- output_samples != all_samples
    collisions <- data.frame(
      sample = all_samples[changed],
      resolution = rep("suffix", sum(changed)),
      output_sample = output_samples[changed],
      stringsAsFactors = FALSE
    )
    return(list(gt = gt, samples = output_samples, collisions = collisions))
  }

  output_samples <- unique(all_samples)
  gt <- matrix(
    NA_character_,
    nrow = nrow(matrices[[1]]),
    ncol = length(output_samples),
    dimnames = list(NULL, output_samples)
  )
  seen <- logical(length(output_samples))
  for (j in seq_along(matrices)) {
    for (k in seq_along(samples[[j]])) {
      target <- match(samples[[j]][k], output_samples)
      incoming <- matrices[[j]][, k]
      if (!seen[target]) {
        gt[, target] <- incoming
        seen[target] <- TRUE
      } else {
        current <- gt[, target]
        current_missing <- vapply(current, is_fully_missing_gt, logical(1))
        incoming_missing <- vapply(incoming, is_fully_missing_gt, logical(1))
        conflict <- !current_missing & !incoming_missing & current != incoming
        if (any(conflict)) {
          cli::cli_abort(
            paste0(
              "Conflicting calls while coalescing sample ",
              "{.val {samples[[j]][k]}}."
            ),
            class = "beagler_error_merge_coalesce_conflict"
          )
        }
        replace <- current_missing & !incoming_missing
        gt[replace, target] <- incoming[replace]
      }
    }
  }
  collisions <- data.frame(
    sample = duplicates,
    resolution = rep("coalesce", length(duplicates)),
    output_sample = duplicates,
    stringsAsFactors = FALSE
  )
  list(gt = gt, samples = output_samples, collisions = collisions)
}

#' @keywords internal
#' @noRd
reindex_gt_call <- function(gt, old_ref, old_alt, new_ref, new_alt) {
  if (is.na(gt) || !nzchar(gt)) return(gt)
  slots <- strsplit(gt, "[/|]")[[1]]
  separators <- regmatches(gt, gregexpr("[/|]", gt))[[1]]
  old_alleles <- vcf_alleles(old_ref, old_alt)
  new_alleles <- vcf_alleles(new_ref, new_alt)
  rewritten <- vapply(slots, function(slot) {
    if (slot == ".") return(".")
    allele <- old_alleles[as.numeric(slot) + 1L]
    index <- match(allele, new_alleles) - 1L
    if (is.na(index)) {
      cli::cli_abort(
        "An allele could not be reindexed during VCF merge.",
        class = "beagler_error_unresolvable_alleles"
      )
    }
    as.character(index)
  }, character(1))
  interleave_gt(rewritten, separators)
}

#' @keywords internal
#' @noRd
interleave_gt <- function(slots, separators) {
  if (length(separators) == 0L) return(slots)
  paste0(
    paste0(c(rbind(slots[-length(slots)], separators)), collapse = ""),
    slots[length(slots)]
  )
}

#' @keywords internal
#' @noRd
is_fully_missing_gt <- function(gt) {
  if (is.na(gt) || !nzchar(gt)) return(TRUE)
  all(strsplit(gt, "[/|]")[[1]] == ".")
}

#' @keywords internal
#' @noRd
site_alignment_key <- function(fixed) {
  paste(
    fixed$CHROM,
    format_vcf_pos(fixed$POS),
    fixed$position_ordinal,
    sep = "\r"
  )
}

#' @keywords internal
#' @noRd
canonical_marker_keys <- function(x) {
  canonical_alt <- vapply(x$fixed$ALT, function(alt) {
    if (alt == ".") {
      "."
    } else {
      paste(sort(strsplit(alt, ",", fixed = TRUE)[[1]]), collapse = ",")
    }
  }, character(1))
  base <- paste(
    x$fixed$CHROM,
    format_vcf_pos(x$fixed$POS),
    x$fixed$REF,
    canonical_alt,
    sep = "\r"
  )
  ordinal <- stats::ave(seq_along(base), base, FUN = seq_along)
  paste(base, ordinal, sep = "\r")
}

#' @keywords internal
#' @noRd
refresh_vcf_ordinals <- function(fixed) {
  fixed <- data.table::copy(fixed)
  fixed$record_ordinal <- ordinal_within(
    fixed,
    c("CHROM", "POS", "REF", "ALT")
  )
  fixed$position_ordinal <- ordinal_within(fixed, c("CHROM", "POS"))
  fixed
}

#' @keywords internal
#' @noRd
subset_genotype_vcf <- function(x, rows) {
  x$fixed <- x$fixed[rows]
  x$gt <- x$gt[rows, , drop = FALSE]
  x
}

#' @keywords internal
#' @noRd
new_genotype_vcf <- function(fixed, gt, samples, meta) {
  structure(
    list(
      fixed = fixed,
      gt = gt,
      meta = meta,
      samples = samples,
      contigs = vcf_contig_order(meta, fixed$CHROM),
      source = character()
    ),
    class = c("beagler_genotype_vcf", "list")
  )
}

#' @keywords internal
#' @noRd
drop_rows_report <- function(fixed, reason) {
  data.frame(
    chrom = fixed$CHROM,
    pos = fixed$POS,
    record_ordinal = fixed$position_ordinal,
    reason = reason,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
bind_drop_reports <- function(reports) {
  if (length(reports) == 0L) return(empty_drop_report())
  unique(do.call(rbind, reports))
}
