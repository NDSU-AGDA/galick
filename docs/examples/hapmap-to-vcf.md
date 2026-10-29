---
title: Convert TASSEL HapMap calls to VCF with Galick
description: Follow a HapMap-to-VCF example using declared marker alleles, IUPAC genotype decoding and validation before publishing the converted file.
---
# HapMap to VCF

TASSEL HapMap (`.hmp.txt`) stores marker alleles inline (`A/T`) and can use
single-character IUPAC codes for genotypes. Galick interprets the first declared
allele as REF and the second as ALT; this alone does not establish reference-genome
orientation.

## A peek at the format

```r
library(galick)

hmp_path <- system.file("extdata", "tiny.hmp.txt", package = "galick")
cat(readLines(hmp_path)[1:3], sep = "\n")
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
stopifnot(result$ok)

result$ok
result$summary$n_samples
result$summary$n_markers
```

The output passes Beagle-target validation:

```r
validate_vcf(out, target = "beagle")$ok
```

## Reference build

If you know the assembly, pass it as `build`:

```r
out <- tempfile(fileext = ".vcf")
hapmap_to_vcf(hmp_path, out, build = "GRCh38")

grep("^##reference", readLines(out), value = TRUE)
```

It records `##reference=` in the VCF header. This labels the file; it does not
convert coordinates or verify alleles against that assembly. The example label
does not identify the assembly of these artificial fixtures.

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
stopifnot(!result$ok, !file.exists(bad_out))
result$ok
result$errors$code
file.exists(bad_out)
```

## Round-trip with CSV

The bundled fixtures `tiny.csv` + `tiny_map.tsv` and `tiny.hmp.txt`
encode the **same** genotypes (one in AHB, one in IUPAC). Both
conversions produce VCFs with identical body content:

```r
csv_out <- tempfile(fileext = ".vcf")
hmp_out <- tempfile(fileext = ".vcf")

csv_to_vcf(
  system.file("extdata", "tiny.csv", package = "galick"),
  system.file("extdata", "tiny_map.tsv", package = "galick"),
  csv_out
)
hapmap_to_vcf(hmp_path, hmp_out)

# Compare body lines (skip headers; they differ on date / source)
csv_body <- readLines(csv_out)
csv_body <- csv_body[!startsWith(csv_body, "#")]
hmp_body <- readLines(hmp_out)
hmp_body <- hmp_body[!startsWith(hmp_body, "#")]

identical(csv_body, hmp_body)
```
