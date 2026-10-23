# Reference

Function reference grouped by area. Every export listed here is
documented with a man page accessible from R via `?function_name`;
the pages here organize them by workflow stage and surface the most
important options.

## At a glance

| Area | What it does | Functions |
|---|---|---|
| [Validate](validate.md) | Inspect inputs before doing anything destructive. | `validate_genotype_table` · `validate_marker_map` · `validate_vcf` · `validate_reference_match` |
| [Convert (to VCF)](convert-to-vcf.md) | Turn breeder tables into VCF v4.2. | `csv_to_vcf` · `hapmap_to_vcf` |
| [Convert (from VCF)](convert-from-vcf.md) | Pull genotypes back into R as matrices. | `vcf_to_dosage` · `vcf_to_hardcalls` |
| [Repair](repair.md) | Fix common VCF issues that block Beagle. | `repair_vcf` |
| [Run Beagle](run-beagle.md) | Imputation and phasing. | `run_beagle` · `install_beagle` · `beagle_version` |
| [Chromosome workflows](chromosome-workflows.md) | Split / run / merge across chromosomes. | `split_vcf_by_chrom` · `run_beagle_by_chrom` · `merge_beagle_outputs` |
| [PLINK bridge](plink-bridge.md) | Convert between VCF and PLINK BED/BIM/FAM. | `vcf_to_plink` · `plink_to_vcf` |
| [Reports](reports.md) | Markdown summaries. | `imputation_report` |
| [Benchmarking](benchmarking.md) | Mask, impute, compare and evaluate configurations. | `mask_genotypes` · `compare_imputed_vs_true` · `benchmark_imputation` · `evaluate_imputation` |
| [Simulation](simulation.md) | Founder-conditioned linked breeding genotypes. | `simulate_vcf` |
| [Panel design](benchmarking.md#design_marker_panel) | Constrained marker selection with independent confirmation. | `design_marker_panel` |

## Object model

Every function returns one of two S3 classes:

- [`beagler_result`](../getting-started/result-classes.md#beagler_result) — for validators, repair, conversion, comparison.
- [`beagler_run`](../getting-started/result-classes.md#beagler_run) — for execution (Beagle, PLINK).

See [The result classes](../getting-started/result-classes.md) for
slot-by-slot detail.
