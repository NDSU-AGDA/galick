# HapMap to VCF

TASSEL HapMap (`.hmp.txt`) is the most common format coming out of
GBS pipelines for plant breeding. It encodes a marker's REF/ALT
inline (`A/T`) and uses single-character IUPAC codes for genotypes.

## A peek at the format

```r
library(beagler)

hmp_path <- system.file("extdata", "tiny.hmp.txt", package = "beagler")
cat(readLines(hmp_path)[1:3], sep = "\n")
#> rs# alleles chrom pos strand assembly# center protLSID assayLSID panel QCcode s1 s2 s3 s4 s5
#> rs1 A/T chr1 100 + NA NA NA NA NA NA T A W T T
#> rs2 C/G chr1 1100 + NA NA NA NA NA NA G C S G C
```

Each cell is a single character:

- the **REF** letter (`A` for `A/T`) means homo-REF (`0/0`)
- the **ALT** letter (`T`) means homo-ALT (`1/1`)
- the matching **IUPAC ambiguity code** (`W` for `A/T`) means
  heterozygous (`0/1`)

| REF / ALT pair | Het code |
|---|---|
| A / G or G / A | `R` |
| C / T or T / C | `Y` |
| C / G or G / C | `S` |
| A / T or T / A | `W` |
| G / T or T / G | `K` |
| A / C or C / A | `M` |

## Convert

```r
out <- tempfile(fileext = ".vcf")
result <- hapmap_to_vcf(hmp_path, out)

result$ok
#> [1] TRUE
result$summary$n_samples
#> [1] 5
result$summary$n_markers
#> [1] 20
```

The output passes Beagle-target validation:

```r
validate_vcf(out, target = "beagle")$ok
#> [1] TRUE
```

## Reference build

If you know the assembly, pass it as `build`:

```r
out <- tempfile(fileext = ".vcf")
hapmap_to_vcf(hmp_path, out, build = "GRCh38")

grep("^##reference", readLines(out), value = TRUE)
#> [1] "##reference=GRCh38"
```

It lands as `##reference=` in the VCF header. Beagle ignores it; it's
there for downstream tools and audit trails.

## Failure modes

Same strict-fail discipline as `csv_to_vcf()`. The function rejects
multiallelic alleles, non-A/C/G/T tokens, and IUPAC codes that don't
match the declared REF/ALT for that marker:

```r
# Marker says alleles=A/T but a sample cell uses 'Y' (C/T het) — error
bad_hmp <- tempfile(fileext = ".hmp.txt")
writeLines(c(
  "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanel\tQCcode\ts1\ts2",
  "rs1\tA/T\tchr1\t100\t+\tNA\tNA\tNA\tNA\tNA\tNA\tY\tA"
), bad_hmp)

bad_out <- tempfile(fileext = ".vcf")
result <- hapmap_to_vcf(bad_hmp, bad_out)
result$ok
#> [1] FALSE
result$errors$code
#> [1] "hapmap_to_vcf_invalid_iupac_token"
file.exists(bad_out)
#> [1] FALSE
```

## Round-trip with CSV

The bundled fixtures `tiny.csv` + `tiny_map.tsv` and `tiny.hmp.txt`
encode the **same** genotypes (one in AHB, one in IUPAC). Both
conversions produce VCFs with identical body content:

```r
csv_out <- tempfile(fileext = ".vcf")
hmp_out <- tempfile(fileext = ".vcf")

csv_to_vcf(
  system.file("extdata", "tiny.csv", package = "beagler"),
  system.file("extdata", "tiny_map.tsv", package = "beagler"),
  csv_out
)
hapmap_to_vcf(hmp_path, hmp_out)

# Compare body lines (skip headers; they differ on date / source)
csv_body <- readLines(csv_out)
csv_body <- csv_body[!startsWith(csv_body, "#")]
hmp_body <- readLines(hmp_out)
hmp_body <- hmp_body[!startsWith(hmp_body, "#")]

identical(csv_body, hmp_body)
#> [1] TRUE
```
