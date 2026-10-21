#' Merge VCF callsets by samples or by markers
#'
#' `by_sample` combines different sample columns for the same marker
#' set. `by_marker` combines different marker rows for the same sample set.
#' Both modes preserve phase-set boundaries and validate reference dictionaries.
#' Sample merges reconcile compatible alleles, recompute AC/AN, union failed
#' filters, and discard other INFO and QUAL values with a report. Merged FORMAT
#' retains GT/PS only; omitted fields are reported. Independent phase blocks
#' receive distinct identifiers. Identical duplicate records retain their
#' occurrence order; other repeated-coordinate mappings must be unambiguous.
#'
#' @param vcfs Character vector of at least two `.vcf` or `.vcf.gz` files.
#' @param output_vcf Output path ending in `.vcf` or `.vcf.gz`.
#' @param mode `"by_sample"` or `"by_marker"`.
#' @param sample_collision For `by_sample`, one of `"error"`, `"suffix"`, or
#'   `"coalesce"`. Coalescing accepts identical calls or fills a fully missing
#'   call from the other file; conflicting called genotypes are refused.
#'   Phased calls in colliding samples are refused because independent phase
#'   orientations cannot be reconciled by selecting one call.
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
  meta <- merge_vcf_metadata(inputs)
  discarded <- merge_discarded_fields(inputs, mode)
  if (mode == "by_sample" && sample_collision == "coalesce") {
    all_samples <- unlist(lapply(inputs, `[[`, "samples"))
    duplicates <- unique(all_samples[duplicated(all_samples)])
    phased_collision <- any(vapply(inputs, function(x) {
      cols <- x$samples %in% duplicates
      any(grepl("|", x$gt[, cols, drop = FALSE], fixed = TRUE))
    }, logical(1)))
    if (phased_collision) {
      cli::cli_abort(
        "Coalescing phased samples requires explicit phase reconciliation.",
        class = "beagler_error_merge_phase_conflict"
      )
    }
  }
  sites_in <- stats::setNames(
    vapply(inputs, function(x) nrow(x$fixed), integer(1)),
    input_names
  )
  samples_in <- stats::setNames(
    vapply(inputs, function(x) length(x$samples), integer(1)),
    input_names
  )

  if (mode == "by_sample") inputs <- align_sample_inputs(inputs)
  unresolved <- drop_unresolvable_sites(inputs, mode)
  inputs <- unresolved$inputs
  ref_mismatch <- drop_ref_mismatch_sites(inputs, mode)
  inputs <- ref_mismatch$inputs
  dropped <- rbind(unresolved$dropped, ref_mismatch$dropped)
  if (all(vapply(inputs, function(x) nrow(x$fixed) == 0L, logical(1)))) {
    return(new_beagler_result(
      errors = issues_tibble(list(list(
        code = "merge_no_records_left",
        message = "No records remain; no output was written.", where = list()
      ))),
      summary = list(
        output_paths = character(), sites_in = sites_in,
        sites_out = 0L, sites_dropped = dropped
      ), call = the_call
    ))
  }
  inputs <- namespace_phase_sets(inputs)

  combined <- if (mode == "by_sample") {
    merge_vcf_by_sample(inputs, sample_collision)
  } else {
    merge_vcf_by_marker(inputs)
  }
  meta <- meta[!grepl("^##FORMAT=<", meta)]
  if (mode == "by_sample") {
    meta <- meta[!grepl("^##INFO=<", meta)]
    meta <- c(
      meta,
      '##INFO=<ID=AC,Number=A,Type=Integer,Description="Called ALT counts">',
      '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total called alleles">'
    )
  }
  if ("PS" %in% names(combined$vcf$format_fields)) {
    meta <- c(meta, paste0(
      "##FORMAT=<ID=PS,Number=1,Type=Integer,",
      'Description="Input-specific phase set">'
    ))
  }
  combined$vcf$meta <- meta
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
  annotations_changed <- length(discarded$format) || length(discarded$info) ||
    mode == "by_sample"
  if (annotations_changed) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "merge_annotations_reconciled",
      message = paste(
        "Merged FORMAT retains GT/PS; sample merges recompute AC/AN",
        "and omit other INFO and QUAL."
      ),
      where = list(
        format_fields = discarded$format, info_fields = discarded$info,
        qual_removed = mode == "by_sample"
      )
    )
  }
  if (nrow(dropped) > 0L) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "merge_sites_dropped",
      message = sprintf(
        "%d site record(s) were dropped; see sites_dropped.",
        nrow(dropped)
      ),
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
      subgenome_labels = detect_subgenome_labels(combined$vcf),
      dropped_format_fields = discarded$format,
      dropped_info_fields = discarded$info
    ),
    call = the_call
  )
}

#' Drop calls that cannot be resolved against their file's REF/ALT
#' @keywords internal
#' @noRd
drop_unresolvable_sites <- function(inputs, mode) {
  bad_by_file <- vector("list", length(inputs))
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    bad <- logical(nrow(x$fixed))
    if (ncol(x$gt) > 0L) {
      for (i in seq_len(nrow(x$fixed))) {
        bad[i] <- any(!vapply(
          x$gt[i, ], valid_gt_for_alt, logical(1),
          alt = x$fixed$ALT[i]
        ))
      }
    }
    bad_by_file[[j]] <- bad
  }

  reports <- list()
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    bad <- if (mode == "by_sample") {
      Reduce(`|`, bad_by_file)
    } else {
      bad_by_file[[j]]
    }
    if (any(bad)) {
      reports[[length(reports) + 1L]] <- drop_rows_report(
        x$fixed[bad], "unresolvable alleles", x$source
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
  if (mode == "by_marker" || nrow(inputs[[1]]$fixed) == 0L) {
    return(list(inputs = inputs, dropped = empty_drop_report()))
  }
  bad <- vapply(seq_len(nrow(inputs[[1]]$fixed)), function(i) {
    length(unique(vapply(inputs, function(x) x$fixed$REF[i], ""))) > 1L
  }, logical(1))
  reports <- list()
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    if (any(bad)) {
      reports[[length(reports) + 1L]] <- drop_rows_report(
        x$fixed[bad], "REF mismatch", x$source
      )
      inputs[[j]] <- subset_genotype_vcf(x, !bad)
    }
  }
  list(inputs = inputs, dropped = bind_drop_reports(reports))
}

#' @keywords internal
#' @noRd
merge_vcf_by_sample <- function(inputs, sample_collision) {
  aligned <- inputs
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
  phase_result <- combine_sample_matrices(
    lapply(aligned, phase_set_matrix),
    lapply(aligned, `[[`, "samples"), sample_collision
  )
  fixed <- reconcile_merged_fixed(fixed, sample_result$gt, aligned)
  fixed <- refresh_vcf_ordinals(fixed)
  vcf <- new_genotype_vcf(
    fixed,
    sample_result$gt,
    sample_result$samples,
    aligned[[1]]$meta,
    format_fields = if (any(phase_result$gt != ".")) {
      list(PS = phase_result$gt)
    } else {
      list()
    }
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
  ps <- do.call(rbind, lapply(inputs, function(x) {
    phase_set_matrix(x)[, match(reference_samples, x$samples), drop = FALSE]
  }))
  fixed <- refresh_vcf_ordinals(fixed)
  list(
    vcf = new_genotype_vcf(fixed, gt, reference_samples, inputs[[1]]$meta,
      format_fields = if (any(ps != ".")) list(PS = ps) else list()
    ),
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
        same_call <- canonical_unphased_gt(current) ==
          canonical_unphased_gt(incoming)
        conflict <- !current_missing & !incoming_missing & !same_call
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
  if (is.na(gt) || !nzchar(gt)) {
    return(gt)
  }
  slots <- strsplit(gt, "[/|]")[[1]]
  separators <- regmatches(gt, gregexpr("[/|]", gt))[[1]]
  old_alleles <- vcf_alleles(old_ref, old_alt)
  new_alleles <- vcf_alleles(new_ref, new_alt)
  rewritten <- vapply(slots, function(slot) {
    if (slot == ".") {
      return(".")
    }
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
  if (length(separators) == 0L) {
    return(slots)
  }
  paste0(
    paste0(c(rbind(slots[-length(slots)], separators)), collapse = ""),
    slots[length(slots)]
  )
}

#' @keywords internal
#' @noRd
is_fully_missing_gt <- function(gt) {
  if (is.na(gt) || !nzchar(gt)) {
    return(TRUE)
  }
  all(strsplit(gt, "[/|]")[[1]] == ".")
}

#' @keywords internal
#' @noRd
align_sample_inputs <- function(inputs) {
  reference <- inputs[[1]]$fixed
  positions <- lapply(inputs, function(x) {
    paste(x$fixed$CHROM, format_vcf_pos(x$fixed$POS), sep = "\r")
  })
  canonical <- lapply(inputs, canonical_marker_keys)
  groups <- lapply(positions, function(keys) split(seq_along(keys), keys))
  for (j in seq_along(inputs)) {
    if (!identical(sort(positions[[1]]), sort(positions[[j]]))) {
      cli::cli_abort(
        "Sample merges require the same positions and record multiplicities.",
        class = "beagler_error_merge_marker_contract"
      )
    }
    index <- integer(nrow(reference))
    for (group in seq_along(groups[[1]])) {
      a <- groups[[1]][[group]]
      b <- groups[[j]][[group]]
      ids_a <- reference$ID[a]
      ids_b <- inputs[[j]]$fixed$ID[b]
      usable_ids <- !any(ids_a %in% c(".", "")) && !anyDuplicated(ids_a) &&
        !anyDuplicated(ids_b) && setequal(ids_a, ids_b)
      if (length(a) == 1L) {
        index[a] <- b
      } else if (usable_ids) {
        index[a] <- b[match(ids_a, ids_b)]
      } else if (setequal(canonical[[1]][a], canonical[[j]][b])) {
        index[a] <- b[match(canonical[[1]][a], canonical[[j]][b])]
      } else {
        cli::cli_abort(
          paste(
            "Repeated-coordinate marker alignment is ambiguous;",
            "normalize or supply unique matching IDs."
          ),
          class = "beagler_error_merge_ambiguous_markers"
        )
      }
    }
    inputs[[j]] <- subset_genotype_vcf(inputs[[j]], index)
  }
  inputs
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
  x$format_fields <- lapply(x$format_fields, function(field) {
    field[rows, , drop = FALSE]
  })
  if (!is.null(x$formats)) x$formats <- x$formats[rows]
  x
}

#' @keywords internal
#' @noRd
new_genotype_vcf <- function(fixed, gt, samples, meta, format_fields = list()) {
  structure(
    list(
      fixed = fixed,
      gt = gt,
      format_fields = format_fields,
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
drop_rows_report <- function(fixed, reason, input_path) {
  data.frame(
    input_path = input_path,
    record_index = fixed$source_record,
    id = fixed$ID,
    ref = fixed$REF,
    alt = fixed$ALT,
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
  if (length(reports) == 0L) {
    return(empty_drop_report())
  }
  do.call(rbind, reports)
}

#' @keywords internal
#' @noRd
phase_set_matrix <- function(x) {
  if (!is.null(x$format_fields$PS)) {
    return(x$format_fields$PS)
  }
  matrix(".", nrow(x$gt), ncol(x$gt), dimnames = dimnames(x$gt))
}

#' @keywords internal
#' @noRd
canonical_unphased_gt <- function(calls) {
  vapply(calls, function(call) {
    if (grepl("|", call, fixed = TRUE)) {
      return(call)
    }
    paste(sort(strsplit(call, "/", fixed = TRUE)[[1]]), collapse = "/")
  }, character(1))
}
