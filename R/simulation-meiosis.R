#' Haldane gametes, with one chromosome sampled from each random bivalent
#' @keywords internal
#' @noRd
simulation_meiosis <- function(haplotypes, map_m, parent_id = "parent", child_id = "child",
                               chrom = "1", copy_offset = 0L) {
  ploidy <- ncol(haplotypes)
  pairs <- if (ploidy == 2L) {
    matrix(1:2, nrow = 1L)
  } else {
    matrix(sample.int(4L), ncol = 2L, byrow = TRUE)
  }
  gamete <- matrix(0L, nrow(haplotypes), ploidy / 2L)
  events <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    length_m <- max(map_m) - min(map_m)
    crossovers <- sort(stats::runif(stats::rpois(1L, length_m), min(map_m), max(map_m)))
    start <- sample.int(2L, 1L)
    selected <- pairs[i, 1L + (start - 1L + findInterval(map_m, crossovers)) %% 2L]
    gamete[, i] <- haplotypes[cbind(seq_len(nrow(haplotypes)), selected)]
    switches <- pairs[i, 1L + (start - 1L + seq.int(0L, length(crossovers))) %% 2L]
    events[[i]] <- data.frame(
      child_id = child_id, parent_id = parent_id, chrom = chrom,
      child_copy = copy_offset + i, event = c("start", rep("crossover", length(crossovers))),
      cm = c(min(map_m), crossovers) * 100, parent_homolog = switches,
      pair_1 = pairs[i, 1], pair_2 = pairs[i, 2], stringsAsFactors = FALSE
    )
  }
  list(haplotypes = gamete, events = data.table::rbindlist(events))
}

#' @keywords internal
#' @noRd
simulation_generation <- function(previous, parent_ids, mating, ids, maps, population, generation,
                                  lineages, event_sink) {
  ploidy <- dim(previous[[1]])[2]
  children <- lapply(previous, function(x) array(0L, c(dim(x)[1:2], length(ids))))
  for (chrom in names(previous)) {
    for (i in seq_along(ids)) {
      for (side in 1:2) {
        parent <- mating[i, side]
        hap <- matrix(previous[[chrom]][, , parent], nrow = dim(previous[[chrom]])[1], ncol = ploidy)
        gamete <- simulation_meiosis(hap, maps[[chrom]], parent_ids[parent], ids[i], chrom,
          copy_offset = (side - 1L) * ploidy / 2L
        )
        for (copy in seq_len(ncol(gamete$haplotypes))) {
          log <- gamete$events[gamete$events$child_copy == (side - 1L) * ploidy / 2L + copy, ]
          positions <- log$cm / 100
          positions[1] <- min(maps[[chrom]])
          inherited <- log$parent_homolog[findInterval(maps[[chrom]], positions)]
          expected <- hap[cbind(seq_len(nrow(hap)), inherited)]
          if (any(expected != gamete$haplotypes[, copy])) {
            simulation_abort("Transmission ledger failed replay validation.")
          }
        }
        copies <- (side - 1L) * ploidy / 2L + seq_len(ploidy / 2L)
        children[[chrom]][, copies, i] <- gamete$haplotypes
        event_sink(gamete$events)
      }
    }
  }
  pedigree <- data.frame(
    sample_id = ids, parent_1 = parent_ids[mating[, 1]],
    parent_2 = parent_ids[mating[, 2]], population = population, generation = generation,
    lineage = lineages, family = paste(parent_ids[mating[, 1]], parent_ids[mating[, 2]], sep = ":"),
    stringsAsFactors = FALSE
  )
  list(haplotypes = children, pedigree = pedigree)
}
