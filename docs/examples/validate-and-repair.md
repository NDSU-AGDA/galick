# Validate and repair

The `validate_vcf` → `repair_vcf` → `validate_vcf` workflow on a
deliberately broken fixture.

## A dirty starting point

`inst/extdata/tiny_dirty.vcf` ships nine records that exercise every
drop condition `repair_vcf()` knows about:

```r
library(galick)

dirty <- system.file("extdata", "tiny_dirty.vcf", package = "galick")
cat(readLines(dirty), sep = "\n")
```

```text
##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID    REF ALT QUAL FILTER INFO FORMAT s1     s2
chr1   200 rs1   A   T   .    PASS   .    GT     0/0    1/1
chr1   100 .     C   G   .    PASS   .    GT     0/1    1/1
chr1   300 rs3   G   A,C .    PASS   .    GT     0/0    0/1
chr1   100 rs4   T   A   .    PASS   .    GT     1/1    0/0
chr2   -5  rs5   A   T   .    PASS   .    GT     0/0    0/1
chr2   50  rs6   G   G   .    PASS   .    GT     0/0    0/1
chr2   100 rs1   C   A   .    PASS   .    GT     0/0    1/1
chr2   200 rs8   T   G   .    PASS   .    GT     0      1
chr2   300 rs9   A   C   .    PASS   .    GT:DP  0/1:30 1/1:25
```

## Validate first

Interactive sessions display the current validation stage and record-level
progress for genotype checks, including elapsed time and ETA. Use
`progress = FALSE` to suppress it or `progress = TRUE` to request it explicitly.
Reading/parsing uses an indeterminate indicator; it does not estimate a file-read
percentage or reduce the validator's memory requirements.

```r
v <- validate_vcf(dirty, target = "beagle")
v$ok

v$errors$code
v$warnings$code
```

## Repair

### Inspect before changing data

Error entries and affected records are different counts. A single issue entry
can refer to many records, and the same record can violate several rules.
`validate_vcf()` groups overlapping findings and counts unique record impact.

For a detailed, read-only investigation, request extended diagnostics:

```r
diagnosis <- validate_vcf(dirty, target = "beagle", diagnostics = "extended")
summary(diagnosis)

evidence <- diagnosis$summary$diagnostics
evidence$issue_impact
evidence$issue_groups
evidence$unsupported_alleles
evidence$allele_usage
evidence$genotype_counts
evidence$genotype_summary
evidence$quality_fields
evidence$quality_summary
evidence$quality_comparisons
evidence$actions
```

For example, `REF=A, ALT=G,-` defines allele indices `0=A`, `1=G`, and `2=-`.
A `2/2` call references the unsupported `-` allele twice, but counts as **one**
dependent genotype cell. A `0/2` or `2|.` call also depends on that allele.
Deleting `,-` alone leaves those calls unresolved. Even an unused allele
requires review of allele-numbered INFO/FORMAT annotations before removal.

Extended mode summarizes only unambiguous `Number=1` numeric FORMAT declarations
(`Integer` or `Float`) on interpretable problem records. Scalar fields with
arbitrary names are supported; undeclared, string and vector-valued fields are
not assigned numeric interpretations. Missing, non-finite and nonnumeric
values are counted separately and excluded from quantiles. Numeric summaries
are descriptive, not probabilities or universal quality thresholds.

`quality_summary` contains pooled exact quantiles by field and genotype category.
`quality_by_record` preserves record-level medians. `quality_comparisons` compares
unsupported-allele calls with fully called supported alleles **at the same
record**, reporting lower/equal/higher median counts and median differences.
Non-comparable records and the number of positive baselines used for ratios are
reported explicitly. Ratios use positive baseline medians only. Low intensity does not prove a
deletion or justify converting a called genotype into missing data.

Coverage and counts describe only the scanned prefix in a partial scan, not the
unexamined file. Record-scoped issue groups exclude file/header findings; those
remain visible in `issue_impact`. Genotype-cell impact counts directly flagged
calls rather than treating every call in a structurally affected row as invalid.
Absent GT fields and malformed record/FORMAT shapes limit inspection and are
reported in the coverage/impact summary.

Standard diagnostics do not collect quality-value arrays. Extended profiling
retains finite values on problem records for exact quantiles, so it uses more
memory. It reuses validation's parsed FORMAT/GT data; it does not reread the file.
The existing parser still loads the full file even for a limited scan.

### Apply a chosen repair

```r
clean <- tempfile(fileext = ".vcf")
result <- repair_vcf(dirty, clean)
result$ok
```

`$summary$step_summary` records what each step did:

```r
result$summary$step_summary
```

Step actions also surface as `$warnings`:

```r
result$warnings$code
```

## Re-validate

```r
validate_vcf(clean, target = "beagle")$ok
```

## What lands in the output

`repair_vcf()` writes a *minimal* VCF: just `CHROM`, `POS`, `ID`,
`REF`, `ALT`, and `FORMAT=GT`. The `GT:DP` cell from the dirty input
gets rewritten to GT-only:

```r
cat(readLines(clean)[6:9], sep = "\n")
```

## Step selection

You can opt out of steps. If the remaining steps can't make the
output pass Beagle validation, the function reports
`repair_vcf_remaining_*` errors and refuses to write:

```r
out <- tempfile(fileext = ".vcf")
result <- repair_vcf(dirty, out, steps = "drop_invalid_records")
result$ok
result$errors$code
file.exists(out)
```

## Chromosome-style normalization

If your VCF mixes `chr1` and `1` for the same chromosome, the default
`chrom_style = "auto"` keeps them separate. Pass `"numeric"` or
`"chr_prefixed"` to normalize labels in the output **and** use the
normalized form for duplicate detection:

```r
# chr1 + 1 with the same POS → seen as duplicates under "numeric"
repair_vcf("messy.vcf", "out.vcf", chrom_style = "numeric")
```

## File-level blockers (no repair)

Some problems can't be safely repaired. They show up as `$errors`
and the function writes nothing:

- `repair_vcf_missing_chrom_header` — no `#CHROM` line.
- `repair_vcf_no_samples` — Beagle needs samples.
- `repair_vcf_duplicated_samples` — no obvious resolution.

Fix these by hand (or upstream) before retrying.
