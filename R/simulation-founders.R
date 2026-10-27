#' @keywords internal
#' @noRd
simulation_abort <- function(message, class = "galick_error_simulation_invalid_input") {
  cli::cli_abort(message, class = class)
}

#' @keywords internal
#' @noRd
simulation_number <- function(x, name, lower = 0, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < lower || x > upper || (integer && x != floor(x))) {
    simulation_abort(paste0(name, " is outside its permitted range."))
  }
}

#' @keywords internal
#' @noRd
simulation_founders <- function(path, map, ploidy, founder_data) {
  validate_vcf_path(path, "founder_vcf")
  checked <- validate_vcf(path)
  if (!checked$ok) simulation_abort("Founder VCF failed general validation.")
  x <- read_genotype_vcf(path)
  if (!nrow(x$gt) || !ncol(x$gt) || any(!x$fixed$REF %in% c("A", "C", "G", "T")) ||
        any(!x$fixed$ALT %in% c("A", "C", "G", "T")) || any(x$fixed$REF == x$fixed$ALT) ||
        anyDuplicated(paste(x$fixed$CHROM, x$fixed$POS))) {
    simulation_abort("Founders require unique biallelic A/C/G/T SNP positions.")
  }
  pattern <- paste0("^[01](\\|[01]){", ploidy - 1L, "}$")
  if (any(!grepl(pattern, x$gt))) simulation_abort("Founder GTs must be complete and phased at the declared ploidy.")
  ps <- x$format_fields$PS
  if (!is.null(ps)) {
    for (chrom in unique(x$fixed$CHROM)) {
      for (j in seq_len(ncol(x$gt))) {
        rows <- which(x$fixed$CHROM == chrom)
        heterozygous <- vapply(
          strsplit(x$gt[rows, j], "|", fixed = TRUE),
          function(a) length(unique(a)) > 1L, logical(1)
        )
        labels <- unique(ps[rows[heterozygous], j])
        if (length(labels) > 1L) simulation_abort("Founder phase blocks cannot be stitched across a chromosome.")
      }
    }
  }
  order <- order_vcf_rows(x)
  x$fixed <- data.table::copy(x$fixed[order])
  x$gt <- x$gt[order, , drop = FALSE]
  x$fixed$ID <- marker_names_from_metadata(data.frame(
    id = x$fixed$ID, chrom = x$fixed$CHROM, pos = x$fixed$POS, ref = x$fixed$REF, alt = x$fixed$ALT
  ))
  if (anyDuplicated(x$fixed$ID)) simulation_abort("Founder marker IDs must be unique after missing-ID fallback.")
  if (is.character(map) && length(map) == 1L && !is.na(map) && file.exists(map)) {
    map <- data.table::fread(map, data.table = FALSE)
  }
  if (!is.data.frame(map) || anyDuplicated(names(map)) ||
        !all(c("chrom", "pos", "cm") %in% names(map))) {
    simulation_abort("map requires chrom, pos and cm columns.")
  }
  map <- as.data.frame(map[c("chrom", "pos", "cm")])
  map$chrom <- as.character(map$chrom)
  if (anyNA(map) || any(!nzchar(map$chrom)) || !is.numeric(map$pos) || !is.numeric(map$cm) ||
        any(!is.finite(map$pos) | map$pos < 1 | map$pos != floor(map$pos)) ||
        any(!is.finite(map$cm) | map$cm < 0)) {
    simulation_abort("Map positions and cM must be finite, valid numbers.")
  }
  keys <- paste(map$chrom, map$pos, sep = "\r")
  expected <- paste(x$fixed$CHROM, x$fixed$POS, sep = "\r")
  if (anyDuplicated(keys) || !setequal(keys, expected)) {
    simulation_abort("Map must cover founder positions exactly once.")
  }
  map <- map[match(expected, keys), , drop = FALSE]
  for (chrom in unique(map$chrom)) {
    if (any(diff(map$cm[map$chrom == chrom]) < 0)) {
      simulation_abort("Genetic positions must be nondecreasing within chromosomes.")
    }
  }
  if (is.null(founder_data)) founder_data <- data.frame(sample_id = x$samples, population = "P1")
  if (!is.data.frame(founder_data) || !all(c("sample_id", "population") %in% names(founder_data))) {
    simulation_abort("founder_data requires sample_id and population.")
  }
  founder_data <- as.data.frame(founder_data[c("sample_id", "population")])
  founder_data[] <- lapply(founder_data, as.character)
  if (anyNA(founder_data) || any(!nzchar(trimws(founder_data$population))) ||
        anyDuplicated(founder_data$sample_id) || !setequal(founder_data$sample_id, x$samples)) {
    simulation_abort("Assign every founder exactly once to a nonempty population label.")
  }
  founder_data <- founder_data[match(x$samples, founder_data$sample_id), , drop = FALSE]
  x$fixed$QUAL <- "."
  x$fixed$FILTER <- "."
  x$fixed$INFO <- "."
  x$meta <- x$meta[grepl("^##fileformat=|^##contig=<|^##reference=", x$meta)]
  list(vcf = x, map = map, metadata = founder_data, ploidy = ploidy)
}

#' @keywords internal
#' @noRd
simulation_counts <- function(x, groups, name, lower) {
  if (length(groups) == 1L && is.null(names(x))) names(x) <- groups
  if (!is.numeric(x) || length(x) != length(groups) || is.null(names(x)) ||
        anyDuplicated(names(x)) || !setequal(names(x), groups) || any(!is.finite(x)) ||
        any(x < lower | x != floor(x) | x > .Machine$integer.max)) {
    simulation_abort(paste0(name, " must provide integer counts for every population."))
  }
  stats::setNames(as.double(x[match(groups, names(x))]), groups)
}

#' @keywords internal
#' @noRd
simulation_alleles <- function(founders) {
  x <- founders$vcf
  indices <- split(seq_len(nrow(x$gt)), factor(x$fixed$CHROM, levels = unique(x$fixed$CHROM)))
  lapply(indices, function(rows) {
    values <- array(0L, c(length(rows), founders$ploidy, ncol(x$gt)))
    for (i in seq_along(rows)) {
      values[i, , ] <- matrix(
        as.integer(unlist(strsplit(x$gt[rows[i], ], "|", fixed = TRUE))),
        founders$ploidy, ncol(x$gt)
      )
    }
    values
  })
}
