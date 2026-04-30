test_that("chrom_style returns chr_prefixed for chr-prefixed labels", {
  expect_equal(beagler:::chrom_style(c("chr1", "chr2", "chr3")), "chr_prefixed")
  expect_equal(beagler:::chrom_style(c("Chr01", "Chr02")), "chr_prefixed")
  expect_equal(beagler:::chrom_style(c("CHR1", "CHRX")), "chr_prefixed")
  expect_equal(beagler:::chrom_style(c("chrX", "chrY", "chrMT")), "chr_prefixed")
})

test_that("chrom_style returns numeric for plain integer-like labels", {
  expect_equal(beagler:::chrom_style(c("1", "2", "3")), "numeric")
  expect_equal(beagler:::chrom_style(c(1, 2, 3)), "numeric")
})

test_that("chrom_style returns letter for letter-only labels", {
  expect_equal(beagler:::chrom_style(c("X", "Y", "MT")), "letter")
})

test_that("chrom_style returns mixed when recognized styles are intermixed", {
  expect_equal(beagler:::chrom_style(c("chr1", "2")), "mixed")
  expect_equal(beagler:::chrom_style(c("1", "X")), "mixed")
  expect_equal(beagler:::chrom_style(c("chr1", "X")), "mixed")
})

test_that("chrom_style returns unknown for wheat-style labels", {
  expect_equal(beagler:::chrom_style(c("1A", "2B", "3D")), "unknown")
})

test_that("chrom_style returns unknown when any label is unrecognized", {
  expect_equal(beagler:::chrom_style(c("chr1", "1A")), "unknown")
  expect_equal(beagler:::chrom_style(c("1", "1A")), "unknown")
})

test_that("chrom_style handles edge cases", {
  expect_equal(beagler:::chrom_style(character()), "unknown")
  expect_equal(beagler:::chrom_style(c(NA_character_, NA_character_)), "unknown")
  expect_equal(beagler:::chrom_style(c("1", NA_character_)), "numeric")
})

test_that("normalize_chrom_label strips chr/Chr/CHR prefix to numeric form", {
  expect_equal(
    beagler:::normalize_chrom_label(c("chr1", "Chr2", "CHR3"), target = "numeric"),
    c("1", "2", "3")
  )
})

test_that("normalize_chrom_label adds lowercase chr prefix", {
  expect_equal(
    beagler:::normalize_chrom_label(c("1", "Chr2", "X"), target = "chr_prefixed"),
    c("chr1", "chr2", "chrX")
  )
})

test_that("normalize_chrom_label preserves NA", {
  expect_equal(
    beagler:::normalize_chrom_label(c("1", NA_character_), target = "numeric"),
    c("1", NA_character_)
  )
})

test_that("normalize_chrom_label as_is returns unchanged", {
  expect_equal(
    beagler:::normalize_chrom_label(c("chr1", "Chr2"), target = "as_is"),
    c("chr1", "Chr2")
  )
})
