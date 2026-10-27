# Galick Roadmap

Galick supports genotype preparation and controlled imputation experiments for
plant and animal breeding. The source namespace and installed help define the
current API; proposed work below is not an availability promise.

## Implemented in 0.1.0

| Area | Public functions |
|---|---|
| Validation | `validate_genotype_table`, `validate_marker_map`, `validate_vcf`, `validate_reference_match` |
| Conversion | `csv_to_vcf`, `hapmap_to_vcf`, `vcf_to_dosage`, `vcf_to_hardcalls` |
| Export and merge | `export_genotypes`, `merge_vcf`, `detect_genotype_ploidy`, `detect_subgenome_labels` |
| Repair | `repair_vcf` |
| Beagle execution | `run_beagle`, `install_beagle`, `beagle_version` |
| Chromosome workflows | `split_vcf_by_chrom`, `run_beagle_by_chrom`, `merge_beagle_outputs` |
| PLINK bridge | `vcf_to_plink`, `plink_to_vcf` |
| Reporting | `imputation_report` |
| Benchmarking | `mask_genotypes`, `compare_imputed_vs_true`, `benchmark_imputation`, `evaluate_imputation` |
| Linked simulation | `simulate_vcf` |
| Marker panel design | `design_marker_panel` |

There are 28 public functions. Core result classes use the `galick_` prefix;
conditions use `galick_error_`. Run records use schema 2.0.0. Simulation,
evaluation and panel-design record schemas remain 1.0.0.

## Release Gates

- Full regression suite and available external-tool integration pass.
- Full `R CMD check`, including vignettes and the manual, has no errors,
  warnings or notes in the verification environment.
- README and website examples parse, match the public API and execute against
  bundled data where declared runnable.
- The strict documentation build succeeds, and diagrams/layouts are inspected
  at desktop and mobile sizes.
- CI results are reviewed on the supported operating systems before tagging.
- Representative breeding datasets support any biological-performance claims.

The GitHub repository is `DidierMurilloF/galick`. Pages deployment and release
tagging are explicit administrative actions, not consequences of a local check.

## Next Priorities

1. Complete real-data plant and animal case studies using independently held-out
   families, populations or generations appropriate to deployment.
2. Measure imputation recovery, missing predictions and group-specific results;
   distinguish them from downstream trait-prediction accuracy.
3. Profile larger regions and reference cohorts before claiming genome-scale
   performance or changing the current in-memory architecture.

## Boundaries and Deferred Work

- Beagle execution, evaluation and panel design remain diploid and biallelic.
- Simulation additionally supports random-bivalent autotetraploid outbred
  inheritance, not general polyploid meiosis or polyploid Beagle imputation.
- Founders and genetic maps determine simulation realism; there is no phenotype
  simulation, species preset, mutation model or selection model.
- Marker-panel design is a bounded heuristic and does not certify assay
  manufacturability, a global optimum or improved breeding-value prediction.
- GWAS, kinship matrices, mixed models, heritability estimation and genomic
  prediction are outside the current analytical scope.

See [the documentation audit](docs/about/release-audit.md) and
[migration guide](docs/getting-started/migration.md) for current contracts and
known limitations.
