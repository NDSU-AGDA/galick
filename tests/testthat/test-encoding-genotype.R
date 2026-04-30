test_that("detect_genotype_encoding identifies AHB", {
  expect_equal(beagler:::detect_genotype_encoding(c("A", "H", "B")), "AHB")
  expect_equal(beagler:::detect_genotype_encoding(c("A", "B")), "AHB")
})

test_that("detect_genotype_encoding identifies 012", {
  expect_equal(beagler:::detect_genotype_encoding(c("0", "1", "2")), "012")
  expect_equal(beagler:::detect_genotype_encoding(c("0", "2")), "012")
})

test_that("detect_genotype_encoding identifies AABB", {
  expect_equal(beagler:::detect_genotype_encoding(c("AA", "AB", "BB")), "AABB")
})

test_that("detect_genotype_encoding identifies VCF", {
  expect_equal(beagler:::detect_genotype_encoding(c("0/0", "0/1", "1/1")), "VCF")
  expect_equal(beagler:::detect_genotype_encoding(c("0/0", "1/0", "1/1")), "VCF")
})

test_that("detect_genotype_encoding identifies REFHET", {
  expect_equal(beagler:::detect_genotype_encoding(c("REF", "HET", "ALT")), "REFHET")
})

test_that("detect_genotype_encoding identifies IUPAC", {
  expect_equal(beagler:::detect_genotype_encoding(c("A", "C", "G", "T")), "IUPAC")
  expect_equal(beagler:::detect_genotype_encoding(c("R", "Y", "S")), "IUPAC")
})

test_that("detect_genotype_encoding tiebreaker picks smallest alphabet that fits", {
  # {A} fits both AHB and IUPAC; AHB is smaller (3 vs 11)
  expect_equal(beagler:::detect_genotype_encoding("A"), "AHB")
})

test_that("detect_genotype_encoding returns 'mixed' when tokens span encodings", {
  # A and 1 are individually recognized but no single encoding holds both
  expect_equal(beagler:::detect_genotype_encoding(c("A", "1")), "mixed")
  expect_equal(beagler:::detect_genotype_encoding(c("AA", "0/0")), "mixed")
  expect_equal(beagler:::detect_genotype_encoding(c("A", "1", "B")), "mixed")
})

test_that("detect_genotype_encoding returns 'unknown' for unrecognized tokens", {
  expect_equal(beagler:::detect_genotype_encoding(c("X", "Q", "Z")), "unknown")
  expect_equal(beagler:::detect_genotype_encoding(c("A", "X")), "unknown")
})

test_that("detect_genotype_encoding returns 'unknown' for empty after strip", {
  expect_equal(beagler:::detect_genotype_encoding(character()), "unknown")
  expect_equal(beagler:::detect_genotype_encoding(c(NA_character_, NA_character_)), "unknown")
  expect_equal(
    beagler:::detect_genotype_encoding(c("", ".", "--", "N"),
                                       missing_codes = c(NA, ".", "--", "", "N")),
    "unknown"
  )
})

test_that("detect_genotype_encoding strips whitespace and uppercases", {
  expect_equal(beagler:::detect_genotype_encoding(c(" a ", "h", "b")), "AHB")
  expect_equal(beagler:::detect_genotype_encoding(c("ref", "het", "alt")), "REFHET")
  expect_equal(beagler:::detect_genotype_encoding(c("aa", " ab ", "bb")), "AABB")
})

test_that("detect_genotype_encoding strips missing codes before detection", {
  # 'N' is in the default missing_codes, so AHB+N input is just AHB after strip
  expect_equal(
    beagler:::detect_genotype_encoding(c("A", "H", "B", "N")),
    "AHB"
  )
  # If N is removed from missing_codes, IUPAC would now contain N; AHB doesn't.
  # Observed = {A, H, B, N}. AHB doesn't fit (N missing). IUPAC doesn't fit
  # (H, B missing). Recognized? H,B in AHB; A,N in IUPAC. -> mixed.
  expect_equal(
    beagler:::detect_genotype_encoding(
      c("A", "H", "B", "N"),
      missing_codes = c(NA, ".", "--", "")
    ),
    "mixed"
  )
})

test_that("detect_genotype_encoding ignores duplicate occurrences", {
  expect_equal(beagler:::detect_genotype_encoding(rep("A", 100L)), "AHB")
})
