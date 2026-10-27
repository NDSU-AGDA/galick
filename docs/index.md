# Galick

**Genotype preparation, simulation and imputation evaluation in R.**

Prepare breeding genotype data for Beagle, inspect conversion losses and validate
imputation experiments. Galick adds linked-genotype simulation, held-out comparison
of Beagle configurations and constrained marker-panel design to the core VCF
workflow. It does not estimate trait effects or breeding values.

[Installation](getting-started/installation.md) ·
[Quick tour](getting-started/quick-tour.md) ·
[Function reference](reference/index.md) ·
[Migration](getting-started/migration.md)

## Start with bundled data

No external tools are needed for this example.

```r
library(galick)

csv <- system.file("extdata", "tiny.csv", package = "galick")
map <- system.file("extdata", "tiny_map.tsv", package = "galick")
vcf <- tempfile(fileext = ".vcf")

conversion <- csv_to_vcf(csv, map, vcf, encoding = "AHB")
stopifnot(conversion$ok)
stopifnot(validate_vcf(vcf, target = "beagle")$ok)

counts <- vcf_to_dosage(vcf)
counts[1:3, 1:4]
```

The fixture has five samples and 20 markers. It illustrates software behavior,
not expected accuracy in a breeding population. `vcf_to_dosage()` reads GT allele
counts, not expected DS dosage.

## Choose a workflow

| Task | Start here |
|---|---|
| Prepare CSV or HapMap | [Conversion](reference/convert-to-vcf.md) |
| Inspect or repair a VCF | [Validation](reference/validate.md), [repair](reference/repair.md) |
| Run and inspect Beagle | [Beagle execution](examples/run-beagle.md) |
| Export or merge genotype data | [Export and merge](examples/export-and-merge.md) |
| Generate linked experimental data | [Simulation](examples/simulate-vcf.md) |
| Compare Beagle configurations | [Held-out evaluation](examples/evaluate-imputation.md) |
| Choose a lower-density SNP panel | [Panel design](examples/design-marker-panel.md) |

## From data to a checked experiment

<div class="diagram" markdown>

```mermaid
flowchart TD
    I[CSV or HapMap] --> V[Convert and validate VCF]
    V --> C{Errors remain?}
    C -->|Yes| P[Inspect issues and request repair]
    P --> V
    C -->|No| B[Run Beagle]
    B --> O[Inspect status and export genotypes]
    F[Phased founders and cM map] --> S[Simulate linked descendants]
    S --> T[Truth, observed and reference VCFs]
    T --> E[Evaluate configurations or design a panel]
    E --> D[Design or tune]
    D --> H[Freeze choice]
    H --> X[Confirm on separate test samples]
    classDef default fill:#f1f7f5,stroke:#55776d,color:#1c302b
    classDef decision fill:#fff4df,stroke:#a57836,color:#493b24
    class C decision
```

</div>

## Understand the result

Validation, conversion and repair return `galick_result`; tool execution returns
`galick_run`. Simulation, evaluation and panel design have dedicated subclasses.
Inspect collected errors and warnings. A successful experiment can have no
recommended winner; `$ok` is not evidence of biological superiority.

See [Result classes](getting-started/result-classes.md) and the
[28-function reference](reference/index.md).

## Scope

Beagle execution, evaluation and panel design require diploid biallelic input.
The simulator also supports restricted random-bivalent autotetraploids, without
enabling polyploid imputation. Most VCF operations are in-memory. Trusted truth,
representative founders and biologically appropriate holdouts remain essential.

Version 0.1.0 is pre-release. Galick is MIT licensed; external tools retain their
own licenses. Source and issues are hosted at
[DidierMurilloF/galick](https://github.com/DidierMurilloF/galick).
