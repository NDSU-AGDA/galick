# Simulation

`simulate_vcf()` models descendants conditional on a phased founder VCF and an
exact cM map. It does not invent founder LD, species histories or novel alleles.
No new runtime dependency is required.

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

`beagler_simulation` inherits from `beagler_result`. Output bundles contain truth,
observed and optional reference VCFs, sample roles/populations, pedigree, replay-checked
homolog switches, a PLINK-format genetic map, observation events and missingness tables,
marker diagnostics, binned dosage LD and a versioned `simulation.json`.

Outputs are staged and validated before publication. Nonempty destinations are
refused. The allele-array ceiling is not a total-RAM or ledger-size bound.
No reference FASTA is checked. LD diagnostics use at most ten adjacent marker lags
per chromosome; undefined correlations are counted, not replaced by zero.

The installed `vignette("simulate-vcf", package = "beagler")` contains runnable
examples. Independent development validation uses AlphaSimR 2.1.0 with `v = 1`,
`p = 0`, `quadProb = 0` and one thread. Results and tolerances are archived in
`tests/testthat/fixtures/simulation-alphasimr.json`; regenerate them with
`tools/validate-simulation-alphasimr.R` in a development environment containing
AlphaSimR and pkgload. Neither is added as a simulation runtime dependency.

`tools/profile-simulation.R` provides a reproducible development workload: 1,000
SNPs on two chromosomes, 20 founders, census 200, three generations, 100 targets
and 50 references. On the development machine it took 11.7 seconds, with 266 MiB
maximum resident memory and a 2.4 MB output bundle. These are workload-specific
measurements, not a whole-genome performance guarantee.
