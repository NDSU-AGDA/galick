# CSV to VCF

Convert a wide breeder genotype table to VCF v4.2 with explicit
REF/ALT mapping from a marker map.

## Inputs

```r
library(beagler)

csv_path <- system.file("extdata", "tiny.csv", package = "beagler")
map_path <- system.file("extdata", "tiny_map.tsv", package = "beagler")

read.csv(csv_path)[1:3, 1:6]
#>   sample_id rs1 rs2 rs3 rs4 rs5
#> 1        s1   B   A   B   H   B
#> 2        s2   A   B   H   B   A
#> 3        s3   H   B   A   B   H

read.delim(map_path)[1:5, ]
#>   marker_id chrom  pos ref alt
#> 1       rs1  chr1  100   A   T
#> 2       rs2  chr1 1100   C   G
#> 3       rs3  chr1 2100   G   A
#> 4       rs4  chr1 3100   T   C
#> 5       rs5  chr1 4100   A   T
```

## Convert

```r
out <- tempfile(fileext = ".vcf")
result <- csv_to_vcf(csv_path, map_path, out)

result$ok
#> [1] TRUE
result$summary$inferred_encoding
#> [1] "AHB"
result$summary$n_samples
#> [1] 5
result$summary$n_markers
#> [1] 20
```

The output VCF passes Beagle-target validation:

```r
validate_vcf(out, target = "beagle")$ok
#> [1] TRUE
```

## Inspect the output

```r
cat(readLines(out)[1:8], sep = "\n")
#> ##fileformat=VCFv4.2
#> ##fileDate=20260428
#> ##source=beagler/0.1.0
#> ##contig=<ID=chr1>
#> ##contig=<ID=chr2>
#> ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#> #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  s1  s2  s3  s4  s5
#> chr1  100  rs1  A  T  .  PASS  .  GT  1/1  0/0  0/1  1/1  1/1
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
#> [1] FALSE
result$errors$code
#> [1] "genotype_table_duplicated_samples"
file.exists(bad_out)
#> [1] FALSE
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
