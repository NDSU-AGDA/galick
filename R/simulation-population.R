#' @keywords internal
#' @noRd
simulation_populations <- function(founders, scenario, parents, n_target, n_reference,
                                   population_size, generations, event_sink) {
  alleles <- simulation_alleles(founders)
  x <- founders$vcf
  maps <- split(founders$map$cm / 100, factor(founders$map$chrom, levels = names(alleles)))
  prefix <- "sim_"
  while (any(startsWith(x$samples, prefix))) prefix <- paste0(prefix, "x")
  pedigree <- list(data.frame(
    sample_id = x$samples, parent_1 = NA_character_, parent_2 = NA_character_,
    population = founders$metadata$population, generation = 0L, lineage = x$samples,
    family = NA_character_, stringsAsFactors = FALSE
  ))
  output <- list()
  groups <- names(n_target)
  for (g in seq_along(groups)) {
    group <- groups[g]
    index <- if (scenario == "outbred") which(founders$metadata$population == group) else match(parents, x$samples)
    current <- lapply(alleles, function(a) a[, , index, drop = FALSE])
    ids <- x$samples[index]
    total <- n_target[g] + n_reference[g]
    final_count <- if (scenario == "outbred") population_size else total
    lineages <- sprintf("%sp%03d_line%06d", prefix, g, seq_len(final_count))
    for (generation in seq_len(generations)) {
      count <- if (scenario == "f2" && generation == 1L) max(2L, final_count) else final_count
      if (scenario != "outbred" && generation == 1L) {
        mating <- matrix(rep(1:2, count), ncol = 2L, byrow = TRUE)
      } else if (scenario == "selfed") {
        mating <- cbind(seq_len(count), seq_len(count))
      } else {
        mating <- t(replicate(count, sample.int(length(ids), 2L, replace = FALSE)))
      }
      next_ids <- sprintf("%sp%03d_g%03d_i%06d", prefix, g, generation, seq_len(count))
      result <- simulation_generation(
        current, ids, mating, next_ids, maps, group, generation,
        rep(lineages, length.out = count), event_sink
      )
      pedigree[[length(pedigree) + 1L]] <- result$pedigree
      current <- result$haplotypes
      ids <- next_ids
    }
    selected <- sample.int(length(ids), total, replace = FALSE)
    output[[group]] <- list(
      haplotypes = lapply(current, function(a) a[, , selected, drop = FALSE]),
      samples = result$pedigree[selected, , drop = FALSE]
    )
    output[[group]]$samples$role <- c(rep("target", n_target[g]), rep("reference", n_reference[g]))
  }
  samples <- data.table::rbindlist(lapply(output, `[[`, "samples"))
  combined <- lapply(alleles, function(a) array(0L, c(dim(a)[1:2], nrow(samples))))
  offset <- 0L
  for (group in output) {
    columns <- offset + seq_len(nrow(group$samples))
    for (chrom in names(combined)) combined[[chrom]][, , columns] <- group$haplotypes[[chrom]]
    offset <- max(columns)
  }
  list(haplotypes = combined, samples = samples, pedigree = data.table::rbindlist(pedigree))
}

#' @keywords internal
#' @noRd
simulation_vcf <- function(founders, population, role) {
  columns <- which(population$samples$role == role)
  fixed <- data.table::copy(founders$vcf$fixed)
  gt <- matrix(character(), nrow(fixed), length(columns))
  for (chrom in names(population$haplotypes)) {
    rows <- which(fixed$CHROM == chrom)
    a <- population$haplotypes[[chrom]]
    for (j in seq_along(columns)) {
      gt[rows, j] <- apply(matrix(a[, , columns[j]], nrow = length(rows)), 1, paste, collapse = "|")
    }
  }
  samples <- population$samples$sample_id[columns]
  colnames(gt) <- samples
  meta <- c(founders$vcf$meta, '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Chromosome-wide simulated phase">')
  new_genotype_vcf(fixed, gt, samples, meta,
    format_fields = list(PS = matrix("1", nrow(gt), ncol(gt)))
  )
}
