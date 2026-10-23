#' Maximum feasible interval packing, preserving required anchors.
#' @keywords internal
#' @noRd
panel_capacity <- function(rows, required, positions, spacing) {
  if (spacing == 0) {
    return(length(rows))
  }
  required <- sort(required)
  if (length(required) > 1L && any(diff(positions[required]) < spacing)) {
    return(-Inf)
  }
  available <- setdiff(rows, required)
  for (anchor in required) available <- available[abs(positions[available] - positions[anchor]) >= spacing]
  last <- -Inf
  count <- length(required)
  for (row in available[order(positions[available])]) {
    if (positions[row] - last >= spacing) {
      count <- count + 1L
      last <- positions[row]
    }
  }
  count
}

#' @keywords internal
#' @noRd
panel_constraints <- function(input, n_markers, spacing) {
  fixed <- input$truth$fixed
  chroms <- unique(fixed$CHROM[input$scoring])
  rows <- lapply(chroms, function(chrom) input$candidates[fixed$CHROM[input$candidates] == chrom])
  names(rows) <- chroms
  positions <- fixed$POS
  capacity <- vapply(rows, function(r) panel_capacity(r, intersect(r, input$mandatory), positions, spacing), numeric(1))
  quotas <- vapply(rows, function(r) max(1L, sum(input$mandatory %in% r)), numeric(1))
  if (any(capacity < quotas) || sum(quotas) > n_markers || sum(capacity) < n_markers) {
    panel_abort("Requested marker count is infeasible under mandatory, chromosome and spacing constraints.")
  }
  weights <- vapply(rows, function(r) max(1, diff(range(positions[r]))), numeric(1))
  while (sum(quotas) < n_markers) {
    active <- which(quotas < capacity)
    ideal <- (n_markers - sum(quotas)) * weights[active] / sum(weights[active])
    add <- pmin(floor(ideal), capacity[active] - quotas[active])
    quotas[active] <- quotas[active] + add
    spare <- n_markers - sum(quotas)
    ranked <- active[order(-(ideal - floor(ideal)), active)]
    for (i in ranked) {
      if (spare > 0 && quotas[i] < capacity[i]) {
        quotas[i] <- quotas[i] + 1L
        spare <- spare - 1L
      }
    }
  }
  list(rows = rows, quotas = quotas, positions = positions, spacing = spacing)
}

#' @keywords internal
#' @noRd
panel_feasible <- function(selected, constraints, complete = FALSE) {
  total <- 0L
  for (i in seq_along(constraints$rows)) {
    rows <- constraints$rows[[i]]
    chosen <- intersect(rows, selected)
    if (length(chosen) > constraints$quotas[i] ||
          panel_capacity(rows, chosen, constraints$positions, constraints$spacing) < constraints$quotas[i]) {
      return(FALSE)
    }
    total <- total + length(chosen)
  }
  !complete || total == sum(constraints$quotas)
}

#' @keywords internal
#' @noRd
panel_ld <- function(input, window, max_pairs) {
  fixed <- input$truth$fixed
  groups <- unique(input$samples$group)
  n_targets <- length(input$scoring)
  neighbors <- vector("list", length(input$candidates))
  pairs <- 0
  for (i in seq_along(input$candidates)) {
    row <- input$candidates[i]
    neighbors[[i]] <- which(
      fixed$CHROM[input$scoring] == fixed$CHROM[row] &
        abs(fixed$POS[input$scoring] - fixed$POS[row]) <= window
    )
    pairs <- pairs + length(neighbors[[i]]) * length(groups)
    if (pairs > max_pairs) panel_abort("max_ld_pairs exceeded; reduce the marker sets or LD window.")
  }
  edges <- vector("list", length(input$candidates))
  names(edges) <- as.character(input$candidates)
  for (i in seq_along(input$candidates)) {
    targets <- integer()
    values <- numeric()
    for (g in seq_along(groups)) {
      metadata <- input$samples
      cols <- match(
        metadata$sample_id[metadata$partition == "design" & metadata$group == groups[g]],
        input$truth$samples
      )
      a <- input$counts[input$candidates[i], cols]
      for (j in neighbors[[i]]) {
        b <- input$counts[input$scoring[j], cols]
        keep <- !is.na(a) & !is.na(b)
        r2 <- if (sum(keep) >= 20L && stats::var(a[keep]) > 0 && stats::var(b[keep]) > 0) {
          stats::cor(a[keep], b[keep])^2
        } else {
          0
        }
        targets <- c(targets, (g - 1L) * n_targets + j)
        values <- c(values, r2)
      }
    }
    edges[[i]] <- list(target = targets, value = values)
  }
  list(edges = edges, n_targets = n_targets * length(groups), n_pairs = pairs)
}

#' @keywords internal
#' @noRd
panel_coverage <- function(rows, ld) {
  coverage <- numeric(ld$n_targets)
  for (row in rows) {
    e <- ld$edges[[as.character(row)]]
    coverage[e$target] <- pmax(coverage[e$target], e$value)
  }
  coverage
}

#' @keywords internal
#' @noRd
panel_proxy <- function(rows, ld) {
  mean(panel_coverage(rows, ld))
}

#' @keywords internal
#' @noRd
panel_gain <- function(row, coverage, ld) {
  e <- ld$edges[[as.character(row)]]
  sum(pmax(coverage[e$target], e$value) - coverage[e$target]) / ld$n_targets
}

#' @keywords internal
#' @noRd
panel_seed <- function(input, constraints, ld = NULL) {
  selected <- input$mandatory
  while (length(selected) < sum(constraints$quotas)) {
    remaining <- setdiff(input$candidates, selected)
    feasible <- if (constraints$spacing == 0) {
      used <- tabulate(match(input$truth$fixed$CHROM[selected], names(constraints$rows)), length(constraints$rows))
      remaining[used[match(input$truth$fixed$CHROM[remaining], names(constraints$rows))] < constraints$quotas[
        match(input$truth$fixed$CHROM[remaining], names(constraints$rows))
      ]]
    } else {
      remaining[vapply(remaining, function(row) panel_feasible(c(selected, row), constraints), logical(1))]
    }
    if (!length(feasible)) panel_abort("No feasible panel completion.")
    score <- if (is.null(ld)) {
      distances <- rep(-Inf, length(feasible))
      for (i in seq_along(constraints$rows)) {
        available <- which(feasible %in% constraints$rows[[i]])
        if (!length(available)) next
        targets <- seq(min(constraints$positions[constraints$rows[[i]]]),
          max(constraints$positions[constraints$rows[[i]]]),
          length.out = constraints$quotas[i] + 2L
        )
        targets <- targets[-c(1L, length(targets))]
        occupied <- intersect(selected, constraints$rows[[i]])
        for (old in occupied) targets <- targets[-which.min(abs(targets - constraints$positions[old]))]
        distances[available] <- vapply(feasible[available], function(row) {
          -min(abs(targets - constraints$positions[row]))
        }, numeric(1))
      }
      distances
    } else {
      coverage <- panel_coverage(selected, ld)
      vapply(feasible, function(row) panel_gain(row, coverage, ld), numeric(1))
    }
    selected <- sort(c(selected, feasible[which.max(score)]))
  }
  selected
}

#' Best unvisited single swap, ranked without tuning or test genotypes.
#' @keywords internal
#' @noRd
panel_proposal <- function(input, constraints, ld, incumbent, visited) {
  best <- NULL
  best_score <- -Inf
  for (old in setdiff(incumbent, input$mandatory)) {
    retained <- setdiff(incumbent, old)
    coverage <- panel_coverage(retained, ld)
    base_score <- mean(coverage)
    chrom <- input$truth$fixed$CHROM[old]
    replacements <- setdiff(input$candidates[input$truth$fixed$CHROM[input$candidates] == chrom], incumbent)
    for (new in replacements) {
      rows <- sort(c(retained, new))
      key <- paste(rows, collapse = ",")
      if (key %in% visited || (constraints$spacing > 0 && !panel_feasible(rows, constraints, TRUE))) next
      score <- base_score + panel_gain(new, coverage, ld)
      if (score > best_score) {
        best <- rows
        best_score <- score
      }
    }
  }
  best
}
