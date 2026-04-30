# Convert from VCF

Pull genotype data out of a VCF as an R matrix for downstream
analysis (PCA, GRM construction, GWAS, ML pipelines, ...).

Both functions accept `.vcf` and `.vcf.gz` and return matrices with
**rows = samples, columns = markers**.

## `vcf_to_dosage()`

```r
vcf_to_dosage(
  vcf,
  samples = NULL,
  markers = NULL
) -> matrix
```

Returns a numeric matrix of ALT-allele counts. For diploid biallelic
sites:

| GT | Dosage |
|---|---|
| `0/0` | `0` |
| `0/1`, `1/0` | `1` |
| `1/1` | `2` |
| `./.` | `NA` |

For multiallelic GTs (e.g. `1/2`), the function counts non-zero
alleles, giving dosage `2` — the conservative interpretation when
ALT-allele identity isn't disambiguated by GT alone.

**Example**

```r
d <- vcf_to_dosage("imputed.vcf.gz")
dim(d)                             # samples × markers
d[1:3, 1:5]
#>    rs1 rs2 rs3 rs4 rs5
#> s1   2   0   1   2   1
#> s2   0   2   1   0   2
#> s3   1   1   2   1   0
```

**Filtering**

```r
# Just three samples, two markers
d <- vcf_to_dosage(vcf,
                   samples = c("s1", "s2", "s3"),
                   markers = c("rs1", "rs5"))
```

Marker columns use the VCF `ID` field. If `ID` is `.`, the function
falls back to `chrom_pos_ref_alt` so the matrix has unique column names.

## `vcf_to_hardcalls()`

```r
vcf_to_hardcalls(
  vcf,
  samples  = NULL,
  markers  = NULL,
  encoding = c("012", "AHB", "AABB")
) -> character matrix
```

Same shape and filtering as `vcf_to_dosage()`, but returns character
hard calls. Mapping is dosage-based:

| Dosage | `"012"` | `"AHB"` | `"AABB"` |
|---|---|---|---|
| 0 | `"0"` | `"A"` | `"AA"` |
| 1 | `"1"` | `"H"` | `"AB"` |
| 2 | `"2"` | `"B"` | `"BB"` |
| NA | `NA` | `NA` | `NA` |

```r
hard <- vcf_to_hardcalls("imputed.vcf.gz", encoding = "AHB")
hard[1:3, 1:5]
#>    rs1 rs2 rs3 rs4 rs5
#> s1   B   A   H   B   H
#> s2   A   B   H   A   B
#> s3   H   H   B   H   A
```
