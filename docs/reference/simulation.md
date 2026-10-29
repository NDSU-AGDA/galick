---
title: Simulate linked breeding genotypes and VCF files in R
description: Generate founder-conditioned diploid crosses, selfed lines and restricted autotetraploid genotypes using Haldane meiosis and genetic maps.
---
# Linked genotype simulation

`simulate_vcf()` models descendants conditional on a phased founder VCF and an
exact cM map. It does not invent founder LD, species histories or novel alleles.
Simulation needs no external software beyond the package's R dependencies.
For a runnable example, use the [linked simulation tutorial](../examples/simulate-vcf.md).
The following call is a template for your own founder VCF and map.

```r
result <- simulate_vcf(
  founder_vcf = "founders.vcf.gz", map = marker_map, # chrom, pos, cm
  output_dir = "simulation", scenario = "outbred",
  n_target = 100L, n_reference = 100L,
  population_size = 500L, generations = 3L, seed = 2026L
)
result$paths
result$observation
result$diagnostics
```

## Models

- Diploid outbred populations: distinct parents sampled within declared populations.
  Census size is not effective population size. No migration or selection is modeled.
- Biparental F2: two homozygous parents, followed by F1 and F2. Parents are generation 0.
- Selfed lineages: F1 followed by selfing; F6 has five selfings, with residual heterozygosity.
- Outbred autotetraploids: random bivalent pairing only, requiring four-copy phased
  founders. No quadrivalents, double reduction or preferential pairing. No allopolyploid
  or mixed-ploidy support. These outputs do not enable polyploid Beagle execution.

Haldane meiosis uses Poisson gametic crossovers without interference. Maps may
contain zero-distance intervals; crossover distances are converted from cM to
Morgans. Founders require complete, chromosome-consistent phase. Multiple
heterozygous phase sets cannot be stitched arbitrarily.

## Observation

An optional observed-marker panel removes off-panel records from observed VCFs.
Calibrated logit missingness supports sample/marker heterogeneity; calling errors
flip one allele copy only in retained called cells. No DS, GP, depth or likelihoods
are fabricated. Observation settings use a separate RNG stream and cannot alter truth.

## Outputs And Validation

`galick_simulation` inherits from `galick_result`. Output bundles contain truth,
observed and optional reference VCFs, sample roles/populations, pedigree, replay-checked
homolog switches, a PLINK-format genetic map, observation events and missingness tables,
marker diagnostics, binned dosage LD and a versioned `simulation.json`.

Outputs are staged and validated before publication. Nonempty destinations are
refused. The allele-array ceiling is not a total-RAM or ledger-size bound.
No reference FASTA is checked. LD diagnostics use at most ten adjacent marker lags
per chromosome; undefined correlations are counted, not replaced by zero.

Development checks compare inheritance summaries against AlphaSimR under
matching model assumptions. The [validation fixtures](https://github.com/NDSU-AGDA/galick/blob/main/tests/testthat/fixtures/simulation-alphasimr.json)
and [comparison script](https://github.com/NDSU-AGDA/galick/blob/main/tools/validate-simulation-alphasimr.R)
are available in the source repository. AlphaSimR is not required to run Galick.
