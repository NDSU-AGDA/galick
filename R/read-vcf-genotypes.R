#' Read a VCF into the genotype export representation (internal)
#'
#' The representation keeps fixed fields in a `data.table`, GT calls in a
#' character matrix, header metadata verbatim, and sample names separately.
#' Additional FORMAT fields are retained in named matrices alongside GT.
#'
#' @import data.table
#' @keywords internal
#' @noRd
read_genotype_vcf <- function(path) {
  parsed <- parse_vcf_raw(path)
  n_expected <- length(parsed$header_cols)

  if (n_expected < 9L) {
    cli::cli_abort(
      "VCF must have the nine fixed columns through FORMAT.",
      class = "galick_error_genotype_vcf_invalid"
    )
  }
  if (anyDuplicated(parsed$samples)) {
    cli::cli_abort(
      "VCF sample names must be unique within a file.",
      class = "galick_error_genotype_vcf_duplicate_samples"
    )
  }

  bad_shape <- which(parsed$body_col_counts != n_expected)
  if (length(bad_shape) > 0L) {
    cli::cli_abort(
      paste0(
        "VCF record column counts do not match the header at record(s) ",
        "{bad_shape}."
      ),
      class = "galick_error_genotype_vcf_invalid"
    )
  }

  n_samples <- length(parsed$samples)
  format_fields <- list()
  formats <- character(length(parsed$body_split))
  if (length(parsed$body_split) == 0L) {
    fixed <- empty_vcf_fixed()
    gt <- matrix(
      character(),
      nrow = 0L,
      ncol = n_samples,
      dimnames = list(NULL, parsed$samples)
    )
  } else {
    body <- do.call(rbind, parsed$body_split)
    pos <- suppressWarnings(as.numeric(body[, 2]))
    bad_pos <- is.na(pos) | !is.finite(pos) | pos < 1 | pos != floor(pos)
    if (any(bad_pos)) {
      cli::cli_abort(
        paste0(
          "VCF POS must contain positive whole numbers at record(s) ",
          "{which(bad_pos)}."
        ),
        class = "galick_error_genotype_vcf_invalid"
      )
    }

    fixed <- data.table::data.table(
      CHROM = body[, 1],
      POS = pos,
      ID = body[, 3],
      REF = body[, 4],
      ALT = body[, 5],
      QUAL = body[, 6],
      FILTER = body[, 7],
      INFO = body[, 8]
    )
    fixed$record_ordinal <- ordinal_within(
      fixed,
      c("CHROM", "POS", "REF", "ALT")
    )
    fixed$position_ordinal <- ordinal_within(fixed, c("CHROM", "POS"))
    fixed$source_record <- seq_len(nrow(fixed))

    gt <- matrix(NA_character_, nrow = nrow(body), ncol = n_samples)
    formats <- body[, 9]
    all_fields <- unique(unlist(strsplit(formats, ":", fixed = TRUE)))
    for (field in setdiff(all_fields, c("GT", "."))) {
      format_fields[[field]] <- matrix(".", nrow(body), n_samples,
        dimnames = list(NULL, parsed$samples)
      )
    }
    if (n_samples > 0L) {
      for (i in seq_len(nrow(body))) {
        record_fields <- strsplit(body[i, 9], ":", fixed = TRUE)[[1]]
        bad_fields <- anyDuplicated(record_fields) ||
          any(record_fields %in% c("", ".")) || endsWith(body[i, 9], ":")
        if (bad_fields) {
          cli::cli_abort("Malformed or duplicated FORMAT fields at record {i}.",
            class = "galick_error_genotype_vcf_invalid"
          )
        }
        gt_index <- match("GT", record_fields)
        if (is.na(gt_index)) {
          cli::cli_abort(
            "VCF FORMAT has no GT field at record {i}.",
            class = "galick_error_genotype_vcf_missing_gt"
          )
        }
        cells <- body[i, 9L + seq_len(n_samples)]
        cell_parts <- strsplit(cells, ":", fixed = TRUE)
        bad_cells <- any(lengths(cell_parts) > length(record_fields)) ||
          any(endsWith(cells, ":")) || any(unlist(cell_parts) == "")
        if (bad_cells) {
          cli::cli_abort("Malformed sample FORMAT values at record {i}.",
            class = "galick_error_genotype_vcf_invalid"
          )
        }
        gt[i, ] <- vapply(cell_parts, function(parts) {
          if (length(parts) < gt_index) NA_character_ else parts[gt_index]
        }, character(1))
        for (field in setdiff(record_fields, "GT")) {
          index <- match(field, record_fields)
          format_fields[[field]][i, ] <- vapply(
            cell_parts, function(parts) {
              if (length(parts) < index) "." else parts[index]
            }, character(1)
          )
        }
      }
    }
    colnames(gt) <- parsed$samples
  }

  structure(
    list(
      fixed = fixed,
      gt = gt,
      format_fields = format_fields,
      formats = formats,
      meta = parsed$meta,
      samples = parsed$samples,
      contigs = vcf_contig_order(parsed$meta, fixed$CHROM),
      source = normalize_path_safe(path, must_exist = TRUE)
    ),
    class = c("galick_genotype_vcf", "list")
  )
}

#' Number records within a set of key columns while retaining input order
#' @keywords internal
#' @noRd
ordinal_within <- function(data, columns) {
  if (nrow(data) == 0L) {
    return(integer())
  }
  key <- do.call(
    paste,
    c(lapply(columns, function(column) data[[column]]), sep = "\r")
  )
  as.integer(stats::ave(seq_len(nrow(data)), key, FUN = seq_along))
}

#' @keywords internal
#' @noRd
empty_vcf_fixed <- function() {
  data.table::data.table(
    CHROM = character(),
    POS = numeric(),
    ID = character(),
    REF = character(),
    ALT = character(),
    QUAL = character(),
    FILTER = character(),
    INFO = character(),
    source_record = integer(),
    record_ordinal = integer(),
    position_ordinal = integer()
  )
}

#' Extract contig IDs from metadata, followed by unlisted body contigs
#' @keywords internal
#' @noRd
vcf_contig_order <- function(meta, body_chrom = character()) {
  contig_lines <- meta[grepl("^##contig=<", meta)]
  ids <- sub("^.*ID=([^,>]+).*$", "\\1", contig_lines)
  ids <- gsub('^"|"$', "", ids)
  ids <- ids[nzchar(ids) & ids != contig_lines]
  unique(c(ids, setdiff(unique(as.character(body_chrom)), ids)))
}

#' Return TRUE for syntactically valid, in-range VCF GT calls
#' @keywords internal
#' @noRd
valid_gt_for_alt <- function(gt, alt) {
  if (is.na(gt) || !nzchar(gt)) {
    return(FALSE)
  }
  if (!grepl("^(\\.|[0-9]+)([/|](\\.|[0-9]+))*$", gt)) {
    return(FALSE)
  }
  slots <- strsplit(gt, "[/|]")[[1]]
  called <- slots[slots != "."]
  if (length(called) == 0L) {
    return(TRUE)
  }
  n_alt <- if (identical(alt, ".") || !nzchar(alt)) {
    0L
  } else {
    length(strsplit(alt, ",", fixed = TRUE)[[1]])
  }
  all(as.numeric(called) <= n_alt)
}

#' Assert that every retained GT call can be resolved against REF/ALT
#' @keywords internal
#' @noRd
assert_resolvable_genotypes <- function(vcf) {
  bad <- logical(nrow(vcf$fixed))
  if (nrow(vcf$fixed) > 0L && ncol(vcf$gt) > 0L) {
    for (i in seq_len(nrow(vcf$fixed))) {
      bad[i] <- any(!vapply(
        vcf$gt[i, ],
        valid_gt_for_alt,
        logical(1),
        alt = vcf$fixed$ALT[i]
      ))
    }
  }
  if (any(bad)) {
    cli::cli_abort(
      paste0(
        "GT contains malformed or out-of-range allele indices at record(s) ",
        "{which(bad)}."
      ),
      class = "galick_error_unresolvable_alleles"
    )
  }
  invisible(vcf)
}
