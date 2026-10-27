# CSV to VCF

Convert a wide breeder genotype table to VCF v4.2 with explicit
REF/ALT mapping from a marker map.

## Inputs

```r
library(galick)

csv_path <- system.file("extdata", "tiny.csv", package = "galick")
map_path <- system.file("extdata", "tiny_map.tsv", package = "galick")

read.csv(csv_path)[1:3, 1:6]

read.delim(map_path)[1:5, ]
```

## Convert

```r
out <- tempfile(fileext = ".vcf")
result <- csv_to_vcf(csv_path, map_path, out)

result$ok
result$summary$inferred_encoding
result$summary$n_samples
result$summary$n_markers
```

The output VCF passes Beagle-target validation:

```r
validate_vcf(out, target = "beagle")$ok
```

## Inspect the output

```r
cat(readLines(out)[1:8], sep = "\n")
```

## Failure modes

If validation fails, no VCF is written — the result carries `$errors`:

```r
bad_csv <- data.frame(
  sample_id = c("s1", "s1", "s2"),     # duplicate sample
  rs1 = c("A", "B", "H"),
  rs2 = c("H", "A", "B")
)
bad_map <- data.frame(
  marker_id = c("rs1", "rs2"),
  chrom = c("chr1", "chr1"),
  pos = c(100, 200),
  ref = c("A", "C"),
  alt = c("T", "G")
)

bad_out <- tempfile(fileext = ".vcf")
result <- csv_to_vcf(bad_csv, bad_map, bad_out, sample_col = "sample_id")

result$ok
result$errors$code
file.exists(bad_out)
```

## Encoding overrides

Auto-detection works for clean data, but you can pin the encoding:

```r
csv_to_vcf(csv_path, map_path, tempfile(fileext = ".vcf"),
           encoding = "AHB")
```

The detector picks the smallest alphabet that fits the observed
tokens, so `{A}` alone resolves to `AHB`, not `IUPAC`. See
[Validate](../reference/validate.md) for the full list.

## REF/ALT mapping conventions

For non-DNA encodings, `ahb_ref` controls which AHB token codes
homo-REF:

```r
# Default: A → 0/0, B → 1/1
csv_to_vcf(csv_path, map_path, tempfile(fileext = ".vcf"),
           encoding = "AHB", ahb_ref = "A")

# Swapped: A → 1/1, B → 0/0
csv_to_vcf(csv_path, map_path, tempfile(fileext = ".vcf"),
           encoding = "AHB", ahb_ref = "B")
```

Same idea for `aabb_ref_homo` with `AABB` encoding. Other encodings
infer REF/ALT from the marker map directly.
