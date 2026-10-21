#' Validate an evaluation experiment before starting any jobs
#' @keywords internal
#' @noRd
evaluation_abort <- function(message) {
  cli::cli_abort(message, class = "beagler_error_evaluation_invalid_input")
}

#' @keywords internal
#' @noRd
evaluation_number <- function(x, name, lower, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x > upper || (integer && x != floor(x))) {
    evaluation_abort(paste0(name, " is outside its permitted range."))
  }
}

#' @keywords internal
#' @noRd
evaluation_candidates <- function(candidates, baseline) {
  ids <- names(candidates)
  if (!is.list(candidates) || length(candidates) < 2L || length(candidates) > 4L ||
        length(ids) != length(candidates) || anyNA(ids) || anyDuplicated(ids) ||
        any(!grepl("^[A-Za-z][A-Za-z0-9_-]*$", ids))) {
    evaluation_abort("Provide two to four uniquely named candidate lists (letters, digits, _ or -).")
  }
  if (!is.character(baseline) || length(baseline) != 1L || is.na(baseline) || !baseline %in% ids) {
    evaluation_abort("baseline must name one candidate.")
  }
  allowed <- c("ne", "em", "err", "burnin", "iterations", "phase-states", "imp-states", "window", "overlap")
  lapply(candidates, function(x) {
    if (!is.list(x) || (length(x) && (is.null(names(x)) || anyDuplicated(names(x)) ||
                                        any(!names(x) %in% c("ref", "map", "extra_args"))))) {
      evaluation_abort("Candidates accept only ref, map and extra_args.")
    }
    args <- if (is.null(x$extra_args)) character() else x$extra_args
    if (!is.character(args) || anyNA(args) || any(!grepl("^[A-Za-z][A-Za-z0-9_-]*=[^\r\n]+$", args))) {
      evaluation_abort("Candidate extra_args must be key=value strings.")
    }
    keys <- sub("=.*$", "", args)
    if (anyDuplicated(keys) || any(!keys %in% allowed)) {
      evaluation_abort("Candidate arguments must be unique, supported model settings; exclusions are forbidden.")
    }
    values <- sub("^[^=]+=", "", args)
    for (i in seq_along(keys)) {
      if (keys[i] == "em") {
        if (!values[i] %in% c("true", "false")) evaluation_abort("em must be true or false.")
      } else {
        number <- suppressWarnings(as.numeric(values[i]))
        integer <- keys[i] %in% c("ne", "burnin", "iterations", "phase-states", "imp-states")
        minimum <- if (keys[i] %in% c("err", "overlap")) 0 else .Machine$double.eps
        evaluation_number(number, keys[i], minimum, integer = integer)
      }
    }
    window <- if ("window" %in% keys) as.numeric(values[match("window", keys)]) else 40
    overlap <- if ("overlap" %in% keys) as.numeric(values[match("overlap", keys)]) else 2
    if (window < 1.1 * overlap) evaluation_abort("window must be at least 1.1 times overlap.")
    if ("ne" %in% keys && !"em=false" %in% args) {
      evaluation_abort("Comparing a fixed ne requires em=false.")
    }
    for (name in c("ref", "map")) {
      value <- x[[name]]
      if (!is.null(value) && (!is.character(value) || length(value) != 1L || is.na(value) ||
                                !file.exists(value) || dir.exists(value))) {
        evaluation_abort(paste0("Candidate ", name, " must be an existing file."))
      }
    }
    list(ref = x$ref, map = x$map, extra_args = args)
  })
}

#' @keywords internal
#' @noRd
evaluation_samples <- function(data, samples) {
  if (!is.data.frame(data) || anyDuplicated(names(data)) ||
        !all(c("sample_id", "partition") %in% names(data)) || !nrow(data)) {
    evaluation_abort("sample_data requires sample_id and partition columns.")
  }
  data <- as.data.frame(
    data[intersect(c("sample_id", "partition", "group", "block"), names(data))],
    stringsAsFactors = FALSE
  )
  for (name in intersect(c("sample_id", "partition", "group", "block"), names(data))) {
    data[[name]] <- as.character(data[[name]])
    if (anyNA(data[[name]]) || any(!nzchar(trimws(data[[name]])))) {
      evaluation_abort("Sample identifiers, partitions, groups and blocks cannot be missing or empty.")
    }
  }
  if (anyDuplicated(data$sample_id) || any(!data$sample_id %in% samples) ||
        !setequal(unique(data$partition), c("tune", "test"))) {
    evaluation_abort("Assign unique truth sample IDs to nonempty tune and test partitions.")
  }
  if ("block" %in% names(data) && any(vapply(
    split(data$partition, data$block),
    function(x) length(unique(x)) != 1L, logical(1)
  ))) {
    evaluation_abort("A declared block cannot cross partitions.")
  }
  if (!"group" %in% names(data)) data$group <- "all"
  data[order(match(data$sample_id, samples)), , drop = FALSE]
}

#' @keywords internal
#' @noRd
evaluation_plan <- function(truth, samples, design, observed_markers, fraction, repeats, seed) {
  observed <- rep(TRUE, nrow(truth$gt))
  if (design == "array") {
    columns <- c("chrom", "pos", "ref", "alt")
    if (!is.data.frame(observed_markers) || !all(columns %in% names(observed_markers)) ||
          !nrow(observed_markers) || anyNA(observed_markers[columns])) {
      evaluation_abort("Array design requires observed_markers with chrom, pos, ref and alt.")
    }
    if (!is.numeric(observed_markers$pos) || any(!is.finite(observed_markers$pos)) ||
          any(observed_markers$pos < 1 | observed_markers$pos != floor(observed_markers$pos))) {
      evaluation_abort("Observed marker positions must be positive whole numbers.")
    }
    keys <- paste(observed_markers$chrom, format_vcf_pos(observed_markers$pos),
      observed_markers$ref, observed_markers$alt,
      sep = "\r"
    )
    if (anyDuplicated(keys) || any(!keys %in% vcf_variant_keys(truth))) {
      evaluation_abort("Observed marker identities must be unique, exact matches to truth.")
    }
    observed <- vcf_variant_keys(truth) %in% keys
    if (!setequal(truth$fixed$CHROM, truth$fixed$CHROM[observed])) {
      evaluation_abort("Each evaluated chromosome requires at least one observed marker.")
    }
  } else if (!is.null(observed_markers)) {
    evaluation_abort("observed_markers is only valid for array design.")
  }
  plans <- list()
  for (partition in c("tune", "test")) {
    metadata <- samples[samples$partition == partition, , drop = FALSE]
    columns <- match(metadata$sample_id, truth$samples)
    known <- !is.na(genotype_counts(truth)[, columns, drop = FALSE])
    for (replicate in seq_len(repeats)) {
      offset <- if (partition == "test") repeats else 0
      job_seed <- as.integer((as.double(seed) + replicate - 1 + offset) %% .Machine$integer.max)
      eligible <- known
      if (design == "sporadic") {
        eligible[, ] <- FALSE
        withr::local_seed(job_seed)
        # Sample within reporting groups so every group has an explicit denominator.
        for (group in unique(metadata$group)) {
          available <- known
          available[, metadata$group != group] <- FALSE
          indices <- which(available)
          count <- as.integer(round(length(indices) * fraction))
          if (!count) evaluation_abort("Every partition/group needs at least one held-out called genotype.")
          eligible[indices[sample.int(length(indices), count)]] <- TRUE
        }
      } else {
        eligible[observed, ] <- FALSE
      }
      if (any(vapply(unique(metadata$group), function(group) {
        sum(eligible[, metadata$group == group, drop = FALSE]) == 0L
      }, logical(1)))) {
        evaluation_abort("Every partition/group needs eligible held-out truth.")
      }
      plans[[paste(partition, replicate, sep = "_")]] <- list(
        partition = partition, replicate = replicate, seed = job_seed, samples = metadata,
        columns = columns, eligible = eligible, observed = observed
      )
    }
  }
  plans
}

#' @keywords internal
#' @noRd
evaluation_input <- function(truth, plan, design) {
  x <- truth
  x$fixed <- data.table::copy(truth$fixed)
  x$samples <- plan$samples$sample_id
  x$gt <- truth$gt[, plan$columns, drop = FALSE]
  x$gt[, ] <- gsub("|", "/", x$gt, fixed = TRUE)
  if (design == "sporadic") x$gt[plan$eligible] <- "./."
  keep <- if (design == "array") plan$observed else rep(TRUE, nrow(x$gt))
  x$fixed <- x$fixed[which(keep)]
  x$gt <- x$gt[keep, , drop = FALSE]
  x$fixed$INFO <- "."
  x$fixed$QUAL <- "."
  x$fixed$FILTER <- "."
  x$format_fields <- list()
  x$formats <- rep("GT", nrow(x$gt))
  x$meta <- x$meta[grepl("^##fileformat=|^##contig=<|^##reference=", x$meta)]
  x
}

#' @keywords internal
#' @noRd
evaluation_manifest <- function(truth, plan) {
  index <- which(plan$eligible, arr.ind = TRUE)
  tibble::tibble(
    sample_id = plan$samples$sample_id[index[, 2]],
    chrom = truth$fixed$CHROM[index[, 1]], pos = truth$fixed$POS[index[, 1]],
    ref = truth$fixed$REF[index[, 1]], alt = truth$fixed$ALT[index[, 1]]
  )
}
