---
title: Worked Galick examples for genotype data in R
description: Follow source-aligned tutorials for VCF conversion, diagnosis, repair, Beagle imputation, linked simulation, accuracy evaluation and marker panels.
---
# Examples

Examples use bundled fixtures or explicitly constructed instructional data.
Paths to real inputs are identified as placeholders. Small artificial examples
demonstrate the API, not biological accuracy or whole-genome scalability.

| Workflow | External requirements |
|---|---|
| [CSV to VCF](csv-to-vcf.md) | None |
| [HapMap to VCF](hapmap-to-vcf.md) | None |
| [Validate and repair](validate-and-repair.md) | None |
| [Run Beagle](run-beagle.md) | Java and Beagle |
| [Chromosome workflows](chromosome-parallel.md) | Java and Beagle; bgzip for compressed merged output |
| [Accuracy benchmark](benchmark-accuracy.md) | Java and Beagle |
| [Export and merge](export-and-merge.md) | bgzip only for compressed output; PLINK for interoperability checks |
| [Linked simulation](simulate-vcf.md) | None |
| [Held-out evaluation](evaluate-imputation.md) | Java and Beagle for execution |
| [Marker panel design](design-marker-panel.md) | Java and Beagle for execution |

The last three tutorials are also installed as R vignettes. Run
`browseVignettes("galick")` to access the built versions.
