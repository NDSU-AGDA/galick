#' Parse structured VCF metadata without splitting quoted commas
#' @keywords internal
#' @noRd
vcf_meta_attributes <- function(line) {
  inner <- sub("^[^=]+=<", "", line)
  inner <- sub(">$", "", inner)
  separator <- ',(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)'
  fields <- strsplit(inner, separator, perl = TRUE)[[1]]
  stats::setNames(sub("^[^=]+=", "", fields), sub("=.*$", "", fields))
}

#' Reconcile reference dictionaries and typed header declarations
#' @keywords internal
#' @noRd
merge_vcf_metadata <- function(inputs) {
  all_meta <- unlist(lapply(inputs, `[[`, "meta"), use.names = FALSE)
  references <- unique(all_meta[startsWith(all_meta, "##reference=")])
  if (length(references) > 1L) {
    cli::cli_abort("Input VCF reference assemblies disagree.",
      class = "beagler_error_merge_header_conflict"
    )
  }
  structured <- grepl("^##(contig|INFO|FORMAT|FILTER|ALT)=<", all_meta)
  combined <- list()
  for (line in unique(all_meta[structured])) {
    kind <- sub("^##([^=]+)=.*$", "\\1", line)
    attrs <- vcf_meta_attributes(line)
    if (!"ID" %in% names(attrs)) {
      cli::cli_abort("Structured VCF metadata has no ID: {line}.",
        class = "beagler_error_merge_header_conflict"
      )
    }
    key <- paste(kind, attrs[["ID"]], sep = ":")
    old <- combined[[key]]
    if (!is.null(old)) {
      compared <- if (kind == "contig") {
        intersect(names(old), names(attrs))
      } else {
        intersect(c("Number", "Type"), intersect(names(old), names(attrs)))
      }
      if (any(old[compared] != attrs[compared])) {
        cli::cli_abort("Incompatible VCF header definitions for {.val {key}}.",
          class = "beagler_error_merge_header_conflict"
        )
      }
      attrs <- c(old, attrs[setdiff(names(attrs), names(old))])
    }
    combined[[key]] <- attrs
  }
  typed <- vapply(names(combined), function(key) {
    kind <- sub(":.*$", "", key)
    attrs <- combined[[key]]
    fields <- paste(names(attrs), attrs, sep = "=", collapse = ",")
    paste0("##", kind, "=<", fields, ">")
  }, character(1))
  unique(c(all_meta[!structured], unname(typed)))
}

#' Report non-GT fields deliberately omitted from merged callsets
#' @keywords internal
#' @noRd
merge_discarded_fields <- function(inputs, mode) {
  formats <- unique(unlist(lapply(inputs, function(x) names(x$format_fields))))
  info <- if (mode == "by_sample") {
    unique(unlist(lapply(inputs, function(x) {
      tokens <- unlist(strsplit(x$fixed$INFO, ";", fixed = TRUE))
      sub("=.*$", "", tokens[tokens != "."])
    })))
  } else {
    character()
  }
  list(format = setdiff(formats, "PS"), info = setdiff(info, c("AC", "AN")))
}

#' Give independent input phase blocks distinct integer identifiers
#' @keywords internal
#' @noRd
namespace_phase_sets <- function(inputs, preserve_fields = FALSE) {
  next_id <- 0L
  for (j in seq_along(inputs)) {
    x <- inputs[[j]]
    old <- x$format_fields$PS
    ps <- matrix(".", nrow(x$gt), ncol(x$gt), dimnames = dimnames(x$gt))
    for (k in seq_along(x$samples)) {
      phased <- grepl("|", x$gt[, k], fixed = TRUE)
      if (!any(phased)) next
      label <- if (is.null(old)) rep(".", nrow(x$gt)) else old[, k]
      keys <- paste(x$fixed$CHROM[phased], label[phased], sep = "\r")
      ids <- match(keys, unique(keys))
      ps[phased, k] <- as.character(next_id + ids)
      next_id <- next_id + length(unique(keys))
    }
    if (preserve_fields) {
      if (any(ps != ".")) {
        inputs[[j]]$format_fields$PS <- ps
        phased_rows <- rowSums(ps != ".") > 0L
        for (row in which(phased_rows)) {
          fields <- strsplit(inputs[[j]]$formats[row], ":", fixed = TRUE)[[1]]
          inputs[[j]]$formats[row] <- paste(unique(c(fields, "PS")), collapse = ":")
        }
      }
    } else {
      inputs[[j]]$format_fields <- if (any(ps != ".")) list(PS = ps) else list()
    }
  }
  inputs
}

#' Recompute called allele counts; retain failures from every merged record
#' @keywords internal
#' @noRd
reconcile_merged_fixed <- function(fixed, gt, aligned) {
  fixed$QUAL <- "."
  for (i in seq_len(nrow(fixed))) {
    alleles <- unlist(strsplit(gt[i, ], "[/|]"))
    alleles <- alleles[!is.na(alleles) & alleles != "."]
    n_alt <- length(vcf_alleles(fixed$REF[i], fixed$ALT[i])) - 1L
    counts <- vapply(seq_len(n_alt), function(k) {
      sum(alleles == as.character(k))
    }, integer(1))
    fixed$INFO[i] <- paste0(
      if (n_alt > 0L) paste0("AC=", paste(counts, collapse = ","), ";") else "",
      "AN=", length(alleles)
    )
    original_filters <- vapply(aligned, function(x) x$fixed$FILTER[i], "")
    filters <- unique(unlist(strsplit(original_filters, ";", fixed = TRUE)))
    failures <- setdiff(filters, c(".", "PASS"))
    fixed$FILTER[i] <- if (length(failures)) {
      paste(failures, collapse = ";")
    } else if ("." %in% filters) {
      "."
    } else {
      "PASS"
    }
  }
  fixed
}
