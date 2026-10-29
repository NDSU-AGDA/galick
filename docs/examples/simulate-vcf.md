---
title: Simulate breeding populations from phased founder VCFs
description: Use a genetic map and phased founders to simulate linked descendants, observed SNP panels, missing calls and calling errors with Galick.
---
# Linked genotype simulation

This worked example is also available as an installed vignette. Small constructed
inputs demonstrate the API; they are not evidence of biological performance.


```r
library(galick)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

`simulate_vcf()` creates linked descendants from supplied phased founders and a
genetic map. Haldane meiosis models transmission; initial LD and population
structure come from founders, not invented species presets. The simulator uses
galick's existing dependencies. AlphaSimR is used only for independent development
validation and is not needed to run this vignette or the package.

## A small founder panel

The following panel is an instructional fixture, not empirical breeding data.
For a biological study supply correctly phased, representative real founders and
their map. Recombination cannot recover variation absent from those founders.

```r
founders <- tempfile(fileext = ".vcf")
map <- data.frame(chrom = "1", pos = seq(100, 1200, 100), cm = seq(0, 55, 5))
records <- vapply(seq_len(nrow(map)), function(i) {
  calls <- if (i <= 6) c("0|0", "1|1", "0|1", "1|0") else c("0|0", "1|1", "1|0", "0|1")
  paste(c("1", map$pos[i], paste0("rs", i), "A", "G", ".", "PASS", ".", "GT", calls), collapse = "\t")
}, character(1))
writeLines(c(
  "##fileformat=VCFv4.2", "##contig=<ID=1>",
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tf1\tf2\tf3\tf4",
  records
), founders)
```

Map units are cM, not Morgans or base pairs. Every founder position must occur
exactly once in the map, with nondecreasing cM within a chromosome. Zero-distance
markers cosegregate. Heterozygous phase blocks cannot be stitched arbitrarily.

## Outbred populations

```r
simulation <- simulate_vcf(
  founders, map, output_dir = tempfile("population_"),
  n_target = 20L, n_reference = 10L,
  population_size = 50L, generations = 3L, seed = 42L,
  missingness = list(rate = 0.1, sample_sd = 0.5, marker_sd = 0.5),
  error_rate = 0.01
)
simulation
simulation$paths
simulation$observation
head(simulation$samples)
```

The census size is not an effective population size. Parents are distinct
individuals; relatives can still mate and descendants need not be unrelated.
For several populations, provide `founder_data` with `sample_id` and `population`,
and matching named `n_target` and `n_reference` vectors. There is no migration.

Truth and reference calls are complete and phased, with disjoint sample IDs.
Observed target calls are unphased. Missingness applies only to retained panel
cells; errors flip one allele copy only in remaining called cells. The logit
intercept is calibrated to the requested expected missing rate. Realized rates
and observation events are reported separately.

## F2 and selfed lineages

```r
f2 <- simulate_vcf(
  founders, map, tempfile("f2_"), scenario = "f2", parents = c("f1", "f2"),
  n_target = 30L, n_reference = 0L, seed = 12L
)
selfed <- simulate_vcf(
  founders, map, tempfile("selfed_"), scenario = "selfed", parents = c("f1", "f2"),
  n_target = 30L, n_reference = 0L, generations = 6L, seed = 12L
)
selfed$summary
```

Both parents must be homozygous at every retained SNP. Founders are generation
0, F1 is generation 1, and F6 has five selfings after F1. Selfed lineages retain
residual heterozygosity; they are not forced into fully homozygous lines.

## An observed marker panel

```r
panel <- data.frame(chrom = "1", pos = c(100, 500, 900), ref = "A", alt = "G")
sparse <- simulate_vcf(
  founders, map, tempfile("array_"), n_target = 10L, n_reference = 10L,
  population_size = 25L, generations = 2L, observed_markers = panel,
  missingness = list(rate = 0), error_rate = 0, seed = 71L
)
validate_vcf(sparse$paths[["observed"]])$summary$n_records
```

Off-panel records are physically absent from observed VCFs but remain in truth
and reference VCFs. Changing the panel, missingness or error settings cannot
change genetic truth for the same founder/model/seed inputs.

For `evaluate_imputation()`, use the complete target truth and assign tune/test
partitions explicitly in the target sample metadata. Choose a deployment-appropriate
holdout; the simulator does not assert family independence or choose a split.

```r
samples <- as.data.frame(simulation$samples)
samples <- samples[samples$role == "target", ]
samples$partition <- rep(c("tune", "test"), each = 10)
study <- evaluate_imputation(
  simulation$paths[["truth"]], samples,
  candidates = list(
    base = list(ref = simulation$paths[["reference"]]),
    longer = list(ref = simulation$paths[["reference"]], extra_args = "iterations=20")
  ), baseline = "base", output_dir = tempfile("evaluation_")
)
```

## Tetraploid boundary

`inheritance = "autotetraploid_bivalent"` requires four-copy phased founders
and supports only `scenario = "outbred"`. Each chromosome randomly forms two
bivalents; one recombinant homolog is transmitted from each. Quadrivalents,
double reduction, preferential pairing, allopolyploids and mixed ploidy are
unsupported. Four-copy output remains valid general VCF but is refused by
`run_beagle()` and `evaluate_imputation()`; generating it does not validate an
imputation method for that biology.

## Diagnostics and limits

`simulation.json` records model version, RNG settings, hashes, map and founder
assignments. Tables report marker allele frequencies and heterozygosity by role
and population, pedigree, replay-checked homolog-switch events, and distance-binned
dosage LD. LD uses up to ten adjacent marker lags per chromosome; it is not an
unbiased survey of all locus pairs. Undefined correlations remain explicit.

The model has no crossover interference, mutation, selection, migration, sex
chromosomes, depth or likelihood simulation. Synthetic REF/ALT and assembly
labels are inherited from the supplied panel; no FASTA is checked. No universal
realism score is produced. Representative founders and a suitable map remain
necessary to make a simulation relevant to a real breeding population.

Bundles are staged and validated before publication. Existing nonempty output
directories are refused. `max_allele_cells` limits planned population arrays,
not total RAM or ledger size. Founder parsing, VCF strings, LD tables and event
logs add overhead. Use moderate workloads and measure resource use before scaling.
