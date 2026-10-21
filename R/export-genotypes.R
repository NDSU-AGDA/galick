#' Export VCF genotypes to common exchange formats
#'
#' VCF preserves FORMAT fields and reference metadata. CSV and TSV retain GT
#' allele slots but not FORMAT/PS phase-block metadata. HapMap uses diploid
#' IUPAC SNP calls; PED requires biallelic
#' markers and at most two allele slots. These formats report phase loss and
#' conversion of partial calls to missing. Compressed VCF output requires bgzip.
#'
#' @param vcf Path to a `.vcf` or `.vcf.gz` file.
#' @param output Output file, or an output prefix for `format = "plink"`.
#' @param format One of `"vcf"`, `"plink"`, `"hapmap"`, `"csv"`, `"tsv"`,
#'   or `"dosage"`.
#' @param downcode Must be `NULL`. Arbitrary allele truncation is unsupported.
#' @param overwrite Refuse to replace output files unless `TRUE`. Replacing a
#'   VCF removes its old `.tbi`/`.csi` indexes; rebuild them for the new file.
#' @param digits Number of significant digits available to numeric table
#'   writers. Values are also written with scientific notation disabled.
#' @param field For dosage output, `"GT"` writes ALT copy counts at any ploidy;
#'   `"DS"` writes expected ALT dosage without falling back to GT. Both require
#'   biallelic markers. Partial GT calls become NA.
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
                             digits = 17L,
                             field = c("GT", "DS")) {
  the_call <- match.call()
  format <- match.arg(format)
  field <- match.arg(field)
  assert_export_output(output)
  if (!is.null(downcode)) {
    cli::cli_abort(
      "Allele truncation is unsupported; keep {.arg downcode} NULL.",
      class = "beagler_error_export_invalid_downcode"
    )
  }
  valid_digits <- is.numeric(digits) && length(digits) == 1L &&
    !is.na(digits) && is.finite(digits) && digits >= 1L && digits <= 22L &&
    digits == floor(digits)
  if (!valid_digits) {
    cli::cli_abort(
      "{.arg digits} must be a whole number between 1 and 22.",
      class = "beagler_error_export_invalid_input"
    )
  }

  x <- read_genotype_vcf(vcf)
  assert_resolvable_genotypes(x)
  gt_ploidy <- detect_genotype_ploidy(x)
  lossy_format <- format %in% c("plink", "hapmap")
  if (lossy_format && !is.na(gt_ploidy) && gt_ploidy > 2L) {
    cli::cli_abort(
      paste0(
        "Format {format} cannot represent genotype ploidy ",
        "{gt_ploidy}. Use VCF or complete allele-count output."
      ),
      class = "beagler_error_export_polyploid_refused"
    )
  }
  multiallelic <- any(grepl(",", x$fixed$ALT, fixed = TRUE))
  if (format %in% c("plink", "hapmap", "dosage") && multiallelic) {
    cli::cli_abort(
      "This export requires biallelic records; split multiallelic sites first.",
      class = "beagler_error_export_multiallelic"
    )
  }
  if (format == "hapmap") {
    alleles <- c(x$fixed$REF, x$fixed$ALT)
    if (any(!alleles %in% c("A", "C", "G", "T"))) {
      cli::cli_abort("HapMap export requires biallelic A/C/G/T SNPs.",
        class = "beagler_error_export_hapmap_alleles"
      )
    }
    slots <- lengths(strsplit(as.character(x$gt), "[/|]"))
    if (any(slots != 2L & as.character(x$gt) != ".")) {
      cli::cli_abort("HapMap export requires diploid calls or missing values.",
        class = "beagler_error_export_hapmap_ploidy"
      )
    }
  }

  paths <- switch(format,
    vcf = write_genotype_vcf(x, output, overwrite = overwrite),
    plink = write_plink_ped_map(
      x, output,
      overwrite = overwrite
    ),
    hapmap = write_genotype_hapmap(x, output, overwrite = overwrite),
    csv = write_genotype_delimited(
      x, output,
      sep = ",", overwrite = overwrite
    ),
    tsv = write_genotype_delimited(
      x, output,
      sep = "\t", overwrite = overwrite
    ),
    dosage = write_genotype_dosage(
      x,
      output,
      overwrite = overwrite,
      digits = as.integer(digits),
      field = field
    )
  )

  warnings <- list()
  partial <- vapply(as.character(x$gt), function(call) {
    slots <- strsplit(call, "[/|]")[[1]]
    any(slots == ".") && any(slots != ".")
  }, logical(1))
  if (lossy_format && any(partial)) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "export_partial_calls_missing",
      message = "Partially missing calls were written as fully missing.",
      where = list(format = format, n_calls = sum(partial))
    )
  }
  if (lossy_format && any(grepl("|", x$gt, fixed = TRUE))) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "export_phase_dropped",
      message = "The selected format does not retain genotype phase.",
      where = list(format = format)
    )
  }
  if (format %in% c("csv", "tsv") && length(x$format_fields)) {
    warnings[[length(warnings) + 1L]] <- list(
      code = "export_format_fields_omitted",
      message = paste(
        "Delimited output retains GT only, not other FORMAT fields",
        "or phase-set metadata."
      ),
      where = list(format_fields = names(x$format_fields))
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
      subgenome_labels = detect_subgenome_labels(x),
      field = if (format == "dosage") field else NA_character_
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
  if (any(dir.exists(paths))) {
    cli::cli_abort("Output paths must be files, not directories.",
      class = "beagler_error_export_invalid_input"
    )
  }
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
  paths <- c(output, paste0(output, c(".tbi", ".csi")))
  assert_output_available(paths, overwrite)

  row_order <- order_vcf_rows(x)
  fixed <- x$fixed[row_order]
  gt <- x$gt[row_order, , drop = FALSE]
  contigs <- vcf_contig_order(x$meta, fixed$CHROM)

  retained_meta <- x$meta[
    !grepl("^##fileformat=|^##contig=<|^##FORMAT=<ID=GT(?:,|>)", x$meta)
  ]
  original_contigs <- x$meta[grepl("^##contig=<", x$meta)]
  original_ids <- vapply(original_contigs, function(line) {
    gsub('^"|"$', "", vcf_meta_attributes(line)[["ID"]])
  }, character(1))
  contig_lines <- vapply(contigs, function(chrom) {
    index <- match(chrom, original_ids)
    if (is.na(index)) {
      sprintf("##contig=<ID=%s>", chrom)
    } else {
      original_contigs[index]
    }
  }, character(1))
  version <- x$meta[startsWith(x$meta, "##fileformat=")]

  header <- c(
    if (length(version)) version[1] else "##fileformat=VCFv4.2",
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
    extra <- names(x$format_fields)
    if (!is.null(x$formats)) {
      fields <- strsplit(x$formats[row_order[i]], ":", fixed = TRUE)[[1]]
      extra <- intersect(extra, fields)
    }
    calls <- gt[i, ]
    for (field in extra) {
      calls <- paste(calls, x$format_fields[[field]][row_order[i], ], sep = ":")
    }
    records[i] <- paste(
      c(
        fixed$CHROM[i], format_vcf_pos(fixed$POS[i]), fixed$ID[i],
        fixed$REF[i], fixed$ALT[i], fixed$QUAL[i], fixed$FILTER[i],
        fixed$INFO[i], paste(c("GT", extra), collapse = ":"), calls
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
  compressed <- grepl("\\.gz$", output, ignore.case = TRUE)
  bgzip <- if (compressed) Sys.which("bgzip") else ""
  if (compressed && !nzchar(bgzip)) {
    cli::cli_abort(
      "BGZF output requires bgzip on PATH; use a plain .vcf otherwise.",
      class = "beagler_error_bgzip_not_found"
    )
  }
  tmp <- tempfile("beagler_", tmpdir = dirname(output), fileext = ".vcf")
  on.exit(unlink(c(tmp, paste0(tmp, ".gz"))), add = TRUE)
  writeLines(lines, tmp, useBytes = TRUE)
  source <- tmp
  if (compressed) {
    result <- run_cli(bgzip, c("-f", tmp))
    source <- paste0(tmp, ".gz")
    if (result$status != 0L || !file.exists(source)) {
      cli::cli_abort("bgzip failed: {result$stderr}",
        class = "beagler_error_bgzip_failed"
      )
    }
  }
  # Stage the old file and indexes together; restore them on failure.
  existing <- c(output, paste0(output, c(".tbi", ".csi")))
  existing <- existing[file.exists(existing)]
  backups <- vapply(existing, function(path) {
    tempfile("beagler_backup_", tmpdir = dirname(path))
  }, character(1))
  moved <- logical(length(existing))
  installed <- FALSE
  on.exit(
    {
      if (installed) {
        unlink(backups)
      } else {
        for (i in which(moved)) {
          if (!file.rename(backups[i], existing[i])) {
            warning("Could not restore output; backup retained at ", backups[i])
          }
        }
      }
    },
    add = TRUE
  )
  for (i in seq_along(existing)) {
    moved[i] <- file.rename(existing[i], backups[i])
    if (!moved[i]) {
      cli::cli_abort("Could not stage existing output {.file {existing[i]}}.",
        class = "beagler_error_export_write_failed"
      )
    }
  }
  if (!file.rename(source, output)) {
    cli::cli_abort(
      "Could not replace {.file {output}} with the completed output.",
      class = "beagler_error_export_write_failed"
    )
  }
  installed <- TRUE
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
      with = FALSE
    ]
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
  x <- subset_genotype_vcf(x, order_vcf_rows(x))

  calls <- x$gt
  if (length(calls) > 0L) {
    for (i in seq_len(nrow(calls))) {
      calls[i, ] <- vapply(
        calls[i, ],
        gt_to_hapmap_iupac,
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
    strand = rep("+", nrow(x$fixed)),
    `assembly#` = rep("NA", nrow(x$fixed)),
    center = rep("NA", nrow(x$fixed)),
    protLSID = rep("NA", nrow(x$fixed)),
    assayLSID = rep("NA", nrow(x$fixed)),
    `panelLSID` = rep("NA", nrow(x$fixed)),
    QCcode = rep("NA", nrow(x$fixed)),
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

#' @keywords internal
#' @noRd
write_genotype_dosage <- function(x, output, overwrite, digits, field) {
  if (!grepl("\\.(?:csv|tsv|txt)$", output, ignore.case = TRUE)) {
    cli::cli_abort(
      "Dosage output must end in .csv, .tsv, or .txt.",
      class = "beagler_error_export_extension"
    )
  }
  assert_output_available(output, overwrite)
  gt <- x$gt
  if (field == "DS" && is.null(x$format_fields$DS)) {
    cli::cli_abort("DS was requested but is absent; no GT fallback is applied.",
      class = "beagler_error_export_missing_ds"
    )
  }
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
      if (field == "GT") {
        dosage[, i] <- vapply(gt[i, ], gt_to_dosage, numeric(1))
      } else {
        tokens <- x$format_fields$DS[i, ]
        values <- suppressWarnings(as.numeric(tokens))
        ploidies <- lengths(strsplit(gt[i, ], "[/|]"))
        invalid <- !is.finite(values) | values < 0 |
          (gt[i, ] != "." & values > ploidies)
        if (any(tokens != "." & invalid)) {
          cli::cli_abort("Invalid DS values at record {i}.",
            class = "beagler_error_export_invalid_ds"
          )
        }
        dosage[, i] <- values
      }
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
  table <- data.frame(
    sample = rownames(formatted), formatted,
    check.names = FALSE, stringsAsFactors = FALSE
  )
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
write_plink_ped_map <- function(x, prefix, overwrite) {
  prefix <- sub("\\.(?:ped|map)$", "", prefix, ignore.case = TRUE)
  paths <- paste0(prefix, c(".ped", ".map"))
  assert_output_available(paths, overwrite)
  if (any(grepl("[[:space:]]", x$samples))) {
    cli::cli_abort(
      "PLINK sample names cannot contain whitespace.",
      class = "beagler_error_export_plink_sample_name"
    )
  }

  gt <- x$gt
  marker_ids <- marker_names_from_metadata(fixed_marker_metadata(x$fixed))
  map <- data.frame(
    chrom = x$fixed$CHROM,
    id = marker_ids,
    cm = rep("0", nrow(x$fixed)),
    pos = vapply(x$fixed$POS, format_vcf_pos, character(1)),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    map, paths[2],
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = FALSE, eol = "\n"
  )

  ped <- matrix("0", nrow = length(x$samples), ncol = 6L + 2L * nrow(x$fixed))
  if (length(x$samples) > 0L) {
    ped[, 1:6] <- cbind(x$samples, x$samples, "0", "0", "0", "-9")
    for (i in seq_len(nrow(x$fixed))) {
      pairs <- lapply(gt[i, ], gt_to_plink_pair,
        ref = x$fixed$REF[i], alt = x$fixed$ALT[i]
      )
      ped[, 6L + 2L * i - 1L] <- vapply(pairs, `[[`, character(1), 1L)
      ped[, 6L + 2L * i] <- vapply(pairs, `[[`, character(1), 2L)
    }
  }
  utils::write.table(
    ped, paths[1],
    sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = FALSE, eol = "\n"
  )
  paths
}

#' @keywords internal
#' @noRd
gt_to_hapmap_iupac <- function(gt, ref, alt) {
  slots <- strsplit(gt, "[/|]")[[1]]
  if (any(slots == ".")) {
    return("N")
  }
  alleles <- sort(unique(vcf_alleles(ref, alt)[as.integer(slots) + 1L]))
  if (length(alleles) == 1L) {
    return(alleles)
  }
  codes <- c(AC = "M", AG = "R", AT = "W", CG = "S", CT = "Y", GT = "K")
  unname(codes[paste0(alleles, collapse = "")])
}

#' @keywords internal
#' @noRd
gt_to_plink_pair <- function(gt, ref, alt) {
  if (is.na(gt) || !nzchar(gt) || gt == ".") {
    return(c("0", "0"))
  }
  slots <- strsplit(gt, "[/|]")[[1]]
  if (any(slots == ".")) {
    return(c("0", "0"))
  }
  if (length(slots) == 1L) slots <- rep(slots, 2L)
  if (length(slots) != 2L) {
    cli::cli_abort(
      "PLINK PED requires diploid calls or homozygous-coded haploid calls.",
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
    input_path = character(),
    record_index = integer(),
    id = character(),
    ref = character(),
    alt = character(),
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
