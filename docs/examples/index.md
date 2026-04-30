# Examples

End-to-end walkthroughs. All examples use the bundled fixtures from
`inst/extdata/` so you can run them as soon as you've installed
beagler.

<div class="grid cards" markdown>

-   :material-file-table:{ .lg .middle } [**CSV to VCF**](csv-to-vcf.md)

    ---

    Convert a wide breeder genotype table into VCF v4.2 with explicit
    REF/ALT mapping. Covers encoding detection, validation,
    failure modes.

-   :material-file-document:{ .lg .middle } [**HapMap to VCF**](hapmap-to-vcf.md)

    ---

    TASSEL HapMap conversion with IUPAC decoding and an optional
    reference build header.

-   :material-shield-check:{ .lg .middle } [**Validate and repair**](validate-and-repair.md)

    ---

    The `validate_vcf` → `repair_vcf` → `validate_vcf` workflow on
    a deliberately broken fixture. Walks through every repair step.

-   :material-rocket-launch:{ .lg .middle } [**Run Beagle**](run-beagle.md)

    ---

    Phase a tiny VCF, inspect the `beagler_run` object, and read
    the v1.0.0 JSON run record.

-   :material-source-branch:{ .lg .middle } [**Chromosome-parallel imputation**](chromosome-parallel.md)

    ---

    Split a multi-chromosome VCF, run Beagle on each chromosome in
    parallel, and merge the outputs.

-   :material-chart-line:{ .lg .middle } [**Benchmark accuracy**](benchmark-accuracy.md)

    ---

    Mask 10 % of genotypes, impute, and compute concordance rate /
    dosage R² / RMSE on the cells we hid.

</div>
