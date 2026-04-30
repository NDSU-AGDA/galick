# M4 recoder tests. The M0 canary skip-helper is gone now that
# recode_to_vcf_gt() is a real implementation.

test_that("recode_to_vcf_gt maps 012 encoding", {
  expect_equal(beagler:::recode_to_vcf_gt(c(0L, 1L, 2L), encoding = "012"),
               c("0/0", "0/1", "1/1"))
})

test_that("recode_to_vcf_gt maps centered encoding", {
  expect_equal(beagler:::recode_to_vcf_gt(c(-1L, 0L, 1L), encoding = "centered"),
               c("0/0", "0/1", "1/1"))
})

test_that("recode_to_vcf_gt maps AHB with default ahb_ref = 'A'", {
  expect_equal(beagler:::recode_to_vcf_gt(c("A", "H", "B"), encoding = "AHB"),
               c("0/0", "0/1", "1/1"))
})

test_that("recode_to_vcf_gt swaps AHB when ahb_ref = 'B'", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("A", "H", "B"), encoding = "AHB", ahb_ref = "B"),
    c("1/1", "0/1", "0/0")
  )
})

test_that("recode_to_vcf_gt maps AABB with default aabb_ref_homo = 'AA'", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("AA", "AB", "BB"), encoding = "AABB"),
    c("0/0", "0/1", "1/1")
  )
})

test_that("recode_to_vcf_gt swaps AABB when aabb_ref_homo = 'BB'", {
  expect_equal(
    beagler:::recode_to_vcf_gt(
      c("AA", "AB", "BB"),
      encoding = "AABB",
      aabb_ref_homo = "BB"
    ),
    c("1/1", "0/1", "0/0")
  )
})

test_that("recode_to_vcf_gt VCF encoding is passthrough", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("0/0", "0/1", "1/0", "1/1"), encoding = "VCF"),
    c("0/0", "0/1", "0/1", "1/1")
  )
})

test_that("recode_to_vcf_gt maps REFHET", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("REF", "HET", "ALT"), encoding = "REFHET"),
    c("0/0", "0/1", "1/1")
  )
})

test_that("recode_to_vcf_gt maps IUPAC for ref=A, alt=G", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("A", "R", "G"), encoding = "IUPAC", ref = "A", alt = "G"),
    c("0/0", "0/1", "1/1")
  )
})

test_that("recode_to_vcf_gt maps IUPAC for ref=C, alt=T", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("C", "Y", "T"), encoding = "IUPAC", ref = "C", alt = "T"),
    c("0/0", "0/1", "1/1")
  )
})

test_that("recode_to_vcf_gt maps NA tokens to ./.", {
  expect_equal(beagler:::recode_to_vcf_gt(NA_integer_, encoding = "012"), "./.")
  expect_equal(beagler:::recode_to_vcf_gt(NA_character_, encoding = "AHB"), "./.")
})

test_that("recode_to_vcf_gt maps explicit missing codes to ./.", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c(".", "--", "N", "./.", ""), encoding = "AHB"),
    rep("./.", 5)
  )
})

test_that("recode_to_vcf_gt rejects out-of-alphabet tokens", {
  expect_error(
    beagler:::recode_to_vcf_gt(3L, encoding = "012"),
    class = "beagler_error_invalid_genotype_token"
  )
  expect_error(
    beagler:::recode_to_vcf_gt(-1L, encoding = "012"),
    class = "beagler_error_invalid_genotype_token"
  )
  expect_error(
    beagler:::recode_to_vcf_gt(2L, encoding = "centered"),
    class = "beagler_error_invalid_genotype_token"
  )
  expect_error(
    beagler:::recode_to_vcf_gt("Z", encoding = "AHB"),
    class = "beagler_error_invalid_genotype_token"
  )
  expect_error(
    beagler:::recode_to_vcf_gt("BAD", encoding = "REFHET"),
    class = "beagler_error_invalid_genotype_token"
  )
})

test_that("recode_to_vcf_gt rejects unknown encoding", {
  expect_error(
    beagler:::recode_to_vcf_gt("A", encoding = "made_up"),
    class = "beagler_error_invalid_genotype_token"
  )
})

test_that("recode_to_vcf_gt IUPAC requires distinct A/C/G/T ref/alt", {
  expect_error(
    beagler:::recode_to_vcf_gt("A", encoding = "IUPAC", ref = "X", alt = "G"),
    class = "beagler_error_invalid_genotype_token"
  )
  expect_error(
    beagler:::recode_to_vcf_gt("A", encoding = "IUPAC", ref = "A", alt = "A"),
    class = "beagler_error_invalid_genotype_token"
  )
  expect_error(
    beagler:::recode_to_vcf_gt("A", encoding = "IUPAC", ref = NULL, alt = "G"),
    class = "beagler_error_invalid_genotype_token"
  )
})

test_that("recode_to_vcf_gt IUPAC rejects token inconsistent with ref/alt", {
  # ref=A, alt=G; "Y" is C/T heterozygous, not consistent
  expect_error(
    beagler:::recode_to_vcf_gt("Y", encoding = "IUPAC", ref = "A", alt = "G"),
    class = "beagler_error_invalid_genotype_token"
  )
})

test_that("recode_to_vcf_gt normalizes whitespace and case", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c(" a ", "h", "b"), encoding = "AHB"),
    c("0/0", "0/1", "1/1")
  )
  expect_equal(
    beagler:::recode_to_vcf_gt(c("ref", " het ", "ALT"), encoding = "REFHET"),
    c("0/0", "0/1", "1/1")
  )
})

test_that("recode_to_vcf_gt preserves order across mixed missing/non-missing input", {
  expect_equal(
    beagler:::recode_to_vcf_gt(c("A", NA, "B", "H"), encoding = "AHB"),
    c("0/0", "./.", "1/1", "0/1")
  )
})

test_that("recode_to_vcf_gt handles empty input", {
  expect_equal(beagler:::recode_to_vcf_gt(character(), encoding = "012"), character())
})
