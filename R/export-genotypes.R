#' Export VCF genotypes to common exchange formats
#'
#' Reads GT calls once and exports them without external command-line tools.
#' VCF, HapMap, CSV, and TSV retain every allele slot and the original phase
#' separator. PLINK PED/MAP and dosage output refuse genotype ploidy above two
#' unless `downcode = "first_two"` is explicitly requested.
#'
#' @param vcf Path to a `.vcf` or `.vcf.gz` file.
#' @param output Output file, or an output prefix for `format = "plink"`.
#' @param format One of `"vcf"`, `"plink"`, `"hapmap"`, `"csv"`, `"tsv"`,
#'   or `"dosage"`.
#' @param downcode `NULL` (the default) or `"first_two"`. The latter keeps
#'   only the first two allele slots for formats that cannot represent higher
#'   ploidy. It is never applied implicitly.
#' @param overwrite Refuse to replace output files unless `TRUE`.
#' @param digits Number of significant digits available to numeric table
#'   writers. Values are also written with scientific notation disabled.
#'
#' @return A `beagler_result` whose summary is the export report.
#' @export
export_genotypes <- function(vcf,
                             output,
                             format = c(
                               "vcf", "plink", "hapmap", "csv", "tsv",
                               "dosage"
                             ),
                             downcode = NULL,
                             overwrite = FALSE,
                             digits = 17L) {
  the_call <- match.call()
  format <- match.arg(format)
  assert_export_output(output)
  if (!is.null(downcode) && !identical(downcode, "first_two")) {
    cli::cli_abort(
      "{.arg downcode} must be NULL or {.val first_two}.",
      class = "beagler_error_export_invalid_downcode"
    )
  }
  if (length(digits) != 1L || is.na(digits) || digits < 1L) {
    cli::cli_abort(
      "{.arg digits} must be one positive number.",
      class = "beagler_error_export_invalid_input"
    )
  }

  x <- read_genotype_vcf(vcf)
  assert_resolvable_genotypes(x)
  gt_ploidy <- detect_genotype_ploidy(x)
  lossy_format <- format %in% c("plink", "dosage")
  polyploid_without_rule <- lossy_format && !is.na(gt_ploidy) &&
    gt_ploidy > 2L && is.null(downcode)
  if (polyploid_without_rule) {
    cli::cli_abort(
      paste0(
        "Format {format} has two allele slots; detected genotype ploidy ",
        "{gt_ploidy}. Supply {.code downcode = \"first_two\"} explicitly."
      ),
      class = "beagler_error_export_polyploid_refused"
    )
  }

  paths <- switch(
    format,
    vcf = write_genotype_vcf(x, output, overwrite = overwrite),
    plink = write_plink_ped_map(
      x, output, downcode = downcode, overwrite = overwrite
    ),
    hapmap = write_genotype_hapmap(x, output, overwrite = overwrite),
    csv = write_genotype_delimited(
      x, output, sep = ",", overwrite = overwrite
    ),
    tsv = write_genotype_delimited(
      x, output, sep = "\t", overwrite = overwrite
    ),
    dosage = write_genotype_dosage(
      x,
      output,
      downcode = downcode,
      overwrite = overwrite,
      digits = as.integer(digits)
    )
  )

  warnings <- list()
  if (lossy_format && !is.na(gt_ploidy) && gt_ploidy > 2L) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "export_downcoded_genotypes",
      message = sprintf(
        "Genotypes of ploidy %d were downcoded using first_two.",
        gt_ploidy
      ),
      where = list(format = format)
    )
  }
  if (format == "plink" && any(grepl("|", x$gt, fixed = TRUE))) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "export_plink_phase_dropped",
      message = "PLINK PED does not retain VCF genotype phase.",
      where = list(format = format)
    )
  }

  new_beagler_result(
    warnings = issues_tibble(warnings),
    summary = list(
      operation = "export",
      format = format,
      input_paths = x$source,
      output_paths = normalizePath(paths, mustWork = TRUE),
      sites_in = nrow(x$fixed),
      sites_out = nrow(x$fixed),
      sites_dropped = empty_drop_report(),
      samples_in = length(x$samples),
      samples_out = length(x$samples),
      collisions = empty_collision_report(),
      genotype_ploidy = gt_ploidy,
      chrom_ploidy = detect_chrom_ploidy(x),
      downcode = if (is.null(downcode)) NA_character_ else downcode
    ),
    call = the_call
  )
}

#' @keywords internal
#' @noRd
assert_export_output <- function(output) {
  invalid <- length(output) != 1L || !is.character(output) ||
    is.na(output) || !nzchar(output)
  if (invalid) {
    cli::cli_abort(
      "{.arg output} must be one non-empty path.",
      class = "beagler_error_export_invalid_input"
    )
  }
  invisible(output)
}

#' @keywords internal
#' @noRd
assert_output_available <- function(paths, overwrite) {
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0L && !isTRUE(overwrite)) {
    cli::cli_abort(
      "Output file(s) already exist: {.file {existing}}.",
      class = "beagler_error_export_output_exists"
    )
  }
  parents <- unique(dirname(paths))
  bad_parent <- !dir.exists(parents) | file.access(parents, 2L) != 0L
  if (any(bad_parent)) {
    cli::cli_abort(
      "Output directory is not writable: {.file {parents[bad_parent]}}.",
      class = "beagler_error_export_invalid_input"
    )
  }
  invisible(paths)
}

#' Write a genotype representation as VCF 4.2
#' @keywords internal
#' @noRd
write_genotype_vcf <- function(x, output, overwrite = FALSE) {
  if (!grepl("\\.vcf(?:\\.gz)?$", output, ignore.case = TRUE)) {
    cli::cli_abort(
      "VCF output must end in {.val .vcf} or {.val .vcf.gz}.",
      class = "beagler_error_export_extension"
    )
  }
  assert_output_available(output, overwrite)

  row_order <- order_vcf_rows(x)
  fixed <- x$fixed[row_order]
  gt <- x$gt[row_order, , drop = FALSE]
  contigs <- vcf_contig_order(x$meta, fixed$CHROM)
  contigs <- contigs[contigs %in% fixed$CHROM]

  retained_meta <- x$meta[
    !grepl("^##fileformat=|^##contig=<|^##FORMAT=<ID=GT(?:,|>)", x$meta)
  ]
  contig_lines <- vapply(contigs, function(chrom) {
    length_value <- existing_contig_length(x$meta, chrom)
    if (is.na(length_value)) {
      length_value <- max(fixed$POS[fixed$CHROM == chrom])
    }
    sprintf(
      "##contig=<ID=%s,length=%s>",
      chrom,
      format_vcf_pos(length_value)
    )
  }, character(1))

  header <- c(
    "##fileformat=VCFv4.2",
    retained_meta,
    contig_lines,
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    paste(
      c(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "INFO", "FORMAT", x$samples
      ),
      collapse = "\t"
    )
  )

  records <- character(nrow(fixed))
  for (i in seq_len(nrow(fixed))) {
    records[i] <- paste(
      c(
        fixed$CHROM[i], format_vcf_pos(fixed$POS[i]), fixed$ID[i],
        fixed$REF[i], fixed$ALT[i], fixed$QUAL[i], fixed$FILTER[i],
        fixed$INFO[i], "GT", gt[i, ]
      ),
      collapse = "\t"
    )
  }
  write_exact_lines(c(header, records), output)
  output
}

#' @keywords internal
#' @noRd
write_exact_lines <- function(lines, output) {
  con <- if (grepl("\\.gz$", output, ignore.case = TRUE)) {
    gzfile(output, "wb")
  } else {
    file(output, "wb")
  }
  on.exit(close(con))
  writeBin(charToRaw(paste0(paste(lines, collapse = "\n"), "\n")), con)
  invisible(output)
}

#' @keywords internal
#' @noRd
order_vcf_rows <- function(x) {
  contigs <- vcf_contig_order(x$meta, x$fixed$CHROM)
  order(match(x$fixed$CHROM, contigs), x$fixed$POS, seq_len(nrow(x$fixed)))
}

#' @keywords internal
#' @noRd
format_vcf_pos <- function(x) {
  format(x, scientific = FALSE, trim = TRUE, digits = 22L)
}

#' @keywords internal
#' @noRd
existing_contig_length <- function(meta, chrom) {
  lines <- meta[grepl("^##contig=<", meta)]
  if (length(lines) == 0L) return(NA_real_)
  ids <- sub("^.*ID=([^,>]+).*$", "\\1", lines)
  ids <- gsub('^"|"$', "", ids)
  line <- lines[match(chrom, ids)]
  if (length(line) == 0L || is.na(line) || !grepl("[,<]length=", line)) {
    return(NA_real_)
  }
  value <- sub("^.*[,<]length=([0-9]+).*$", "\\1", line)
  suppressWarnings(as.numeric(value))
}

#' @keywords internal
#' @noRd
write_genotype_delimited <- function(x, output, sep, overwrite = FALSE) {
  expected <- if (identical(sep, ",")) "\\.csv$" else "\\.(?:tsv|txt)$"
  if (!grepl(expected, output, ignore.case = TRUE)) {
    cli::cli_abort(
      "Output extension does not match the requested delimited format.",
      class = "beagler_error_export_extension"
    )
  }
  assert_output_available(output, overwrite)
  table <- as.data.frame(
    x$fixed[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"),
            with = FALSE]
  )
  table$POS <- vapply(table$POS, format_vcf_pos, character(1))
  if (length(x$samples) > 0L) {
    table <- cbind(table, as.data.frame(x$gt, stringsAsFactors = FALSE))
  }
  utils::write.table(
    table,
    output,
    sep = sep,
    row.names = FALSE,
    col.names = TRUE,
    quote = identical(sep, ","),
    na = ".",
    eol = "\n"
  )
  output
}

#' @keywords internal
#' @noRd
write_genotype_hapmap <- function(x, output, overwrite = FALSE) {
  if (!grepl("\\.(?:hmp\\.txt|hmp|txt|tsv)$", output, ignore.case = TRUE)) {
    cli::cli_abort(
      "HapMap output must end in .hmp.txt, .hmp, .txt, or .tsv.",
      class = "beagler_error_export_extension"
    )
  }
  assert_output_available(output, overwrite)

  calls <- x$gt
  if (length(calls) > 0L) {
    for (i in seq_len(nrow(calls))) {
      calls[i, ] <- vapply(
        calls[i, ],
        gt_indices_to_alleles,
        character(1),
        ref = x$fixed$REF[i],
        alt = x$fixed$ALT[i]
      )
    }
  }
  table <- data.frame(
    `rs#` = x$fixed$ID,
    alleles = vapply(seq_len(nrow(x$fixed)), function(i) {
      paste(vcf_alleles(x$fixed$REF[i], x$fixed$ALT[i]), collapse = "/")
    }, character(1)),
    chrom = x$fixed$CHROM,
    pos = vapply(x$fixed$POS, format_vcf_pos, character(1)),
    strand = "+",
    `assembly#` = "NA",
    center = "NA",
    protLSID = "NA",
    assayLSID = "NA",
    panel = "NA",
    QCcode = "NA",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (length(x$samples) > 0L) {
    table <- cbind(table, as.data.frame(calls, stringsAsFactors = FALSE))
  }
  utils::write.table(
    table,
    output,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ".",
    eol = "\n"
  )
  output
}

#' @keywords internal
#' @noRd
vcf_alleles <- function(ref, alt) {
  c(ref, if (identical(alt, ".") || !nzchar(alt)) {
    character()
  } else {
    strsplit(alt, ",", fixed = TRUE)[[1]]
  })
}

#' Replace VCF allele indices with allele strings while retaining separators
#' @keywords internal
#' @noRd
gt_indices_to_alleles <- function(gt, ref, alt) {
  if (is.na(gt) || !nzchar(gt)) return(".")
  slots <- strsplit(gt, "[/|]")[[1]]
  separators <- regmatches(gt, gregexpr("[/|]", gt))[[1]]
  alleles <- vcf_alleles(ref, alt)
  decoded <- vapply(slots, function(slot) {
    if (slot == ".") "." else alleles[as.numeric(slot) + 1L]
  }, character(1))
  if (length(separators) == 0L) return(decoded)
  paste0(
    paste0(
      c(rbind(decoded[-length(decoded)], separators)),
      collapse = ""
    ),
    decoded[length(decoded)]
  )
}

#' @keywords internal
#' @noRd
write_genotype_dosage <- function(x, output, downcode, overwrite, digits) {
  if (!grepl("\\.(?:csv|tsv|txt)$", output, ignore.case = TRUE)) {
    cli::cli_abort(
      "Dosage output must end in .csv, .tsv, or .txt.",
      class = "beagler_error_export_extension"
    )
  }
  assert_output_available(output, overwrite)
  gt <- downcode_gt_matrix(x$gt, downcode)
  dosage <- matrix(
    NA_real_,
    nrow = length(x$samples),
    ncol = nrow(x$fixed),
    dimnames = list(x$samples, marker_names_from_metadata(
      fixed_marker_metadata(x$fixed)
    ))
  )
  if (length(dosage) > 0L) {
    for (i in seq_len(nrow(x$fixed))) {
      dosage[, i] <- vapply(gt[i, ], gt_to_dosage, numeric(1))
    }
  }
  if (ncol(dosage) == 0L) {
    formatted <- matrix(
      character(),
      nrow = nrow(dosage),
      ncol = 0L,
      dimnames = dimnames(dosage)
    )
  } else {
    formatted <- apply(dosage, 2L, function(values) {
      ifelse(
        is.na(values),
        "NA",
        format(values, digits = digits, scientific = FALSE, trim = TRUE)
      )
    })
    formatted <- matrix(
      formatted,
      nrow = nrow(dosage),
      dimnames = dimnames(dosage)
    )
  }
  table <- data.frame(sample = rownames(formatted), formatted,
                      check.names = FALSE, stringsAsFactors = FALSE)
  sep <- if (grepl("\\.csv$", output, ignore.case = TRUE)) "," else "\t"
  utils::write.table(
    table,
    output,
    sep = sep,
    row.names = FALSE,
    col.names = TRUE,
    quote = identical(sep, ","),
    na = "NA",
    eol = "\n"
  )
  output
}

#' @keywords internal
#' @noRd
write_plink_ped_map <- function(x, prefix, downcode, overwrite) {
  prefix <- sub("\\.(?:ped|map)$", "", prefix, ignore.case = TRUE)
  paths <- paste0(prefix, c(".ped", ".map"))
  assert_output_available(paths, overwrite)
  if (any(grepl("[[:space:]]", x$samples))) {
    cli::cli_abort(
      "PLINK sample names cannot contain whitespace.",
      class = "beagler_error_export_plink_sample_name"
    )
  }

  gt <- downcode_gt_matrix(x$gt, downcode)
  marker_ids <- marker_names_from_metadata(fixed_marker_metadata(x$fixed))
  map <- data.frame(
    chrom = x$fixed$CHROM,
    id = marker_ids,
    cm = "0",
    pos = vapply(x$fixed$POS, format_vcf_pos, character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    map, paths[2], sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = FALSE, eol = "\n"
  )

  ped <- matrix("0", nrow = length(x$samples), ncol = 6L + 2L * nrow(x$fixed))
  if (length(x$samples) > 0L) {
    ped[, 1:6] <- cbind(x$samples, x$samples, "0", "0", "0", "-9")
    for (i in seq_len(nrow(x$fixed))) {
      pairs <- lapply(gt[i, ], gt_to_plink_pair,
                      ref = x$fixed$REF[i], alt = x$fixed$ALT[i])
      ped[, 6L + 2L * i - 1L] <- vapply(pairs, `[[`, character(1), 1L)
      ped[, 6L + 2L * i] <- vapply(pairs, `[[`, character(1), 2L)
    }
  }
  utils::write.table(
    ped, paths[1], sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = FALSE, eol = "\n"
  )
  paths
}

#' @keywords internal
#' @noRd
downcode_gt_matrix <- function(gt, downcode) {
  if (!identical(downcode, "first_two")) return(gt)
  out <- gt
  out[] <- vapply(as.character(gt), function(call) {
    if (is.na(call) || !nzchar(call) || call == ".") return(call)
    slots <- strsplit(call, "[/|]")[[1]]
    if (length(slots) <= 2L) return(call)
    separator <- regmatches(call, regexpr("[/|]", call))
    paste(slots[1:2], collapse = separator)
  }, character(1))
  out
}

#' @keywords internal
#' @noRd
gt_to_plink_pair <- function(gt, ref, alt) {
  if (is.na(gt) || !nzchar(gt) || gt == ".") return(c("0", "0"))
  slots <- strsplit(gt, "[/|]")[[1]]
  if (length(slots) == 1L) slots <- rep(slots, 2L)
  if (length(slots) != 2L) {
    cli::cli_abort(
      "PLINK PED requires exactly two allele slots after downcoding.",
      class = "beagler_error_export_plink_ploidy"
    )
  }
  alleles <- vcf_alleles(ref, alt)
  vapply(slots, function(slot) {
    if (slot == ".") "0" else alleles[as.numeric(slot) + 1L]
  }, character(1))
}

#' @keywords internal
#' @noRd
fixed_marker_metadata <- function(fixed) {
  data.frame(
    chrom = fixed$CHROM,
    pos = fixed$POS,
    id = fixed$ID,
    ref = fixed$REF,
    alt = fixed$ALT,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
empty_drop_report <- function() {
  data.frame(
    chrom = character(),
    pos = numeric(),
    record_ordinal = integer(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
empty_collision_report <- function() {
  data.frame(
    sample = character(),
    resolution = character(),
    output_sample = character(),
    stringsAsFactors = FALSE
  )
}
