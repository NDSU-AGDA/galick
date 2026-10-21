test_that("VCF export preserves phase sets and additional FORMAT values", {
  input <- scientific_vcf("s", c(
    "1\t10\ta\tA\tG\t.\tPASS\t.\tGT:PS:DS\t0|1:10:0.8",
    "1\t20\tb\tA\tG\t.\tPASS\t.\tGT:PS:DS\t1|0:20:1.2"
  ), c(
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">",
    "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">"
  ))
  out <- withr::local_tempfile(fileext = ".vcf")
  export_genotypes(input, out, "vcf")
  expect_true(any(grepl("0|1:10:0.8", readLines(out), fixed = TRUE)))
  expect_true(any(grepl("1|0:20:1.2", readLines(out), fixed = TRUE)))
})

test_that("VCF export retains reference dictionaries, not guessed lengths", {
  input <- scientific_vcf(
    "s", "1\t10\ta\tA\tG\t.\tPASS\t.\tGT\t0/1",
    c("##contig=<ID=1,assembly=build1,md5=abc>", "##contig=<ID=2,length=500>")
  )
  out <- withr::local_tempfile(fileext = ".vcf")
  export_genotypes(input, out, "vcf")
  expect_true("##contig=<ID=1,assembly=build1,md5=abc>" %in% readLines(out))
  expect_true("##contig=<ID=2,length=500>" %in% readLines(out))
  expect_false(any(grepl("ID=1,length=10", readLines(out), fixed = TRUE)))
})

test_that("polyploid allele counts are invariant to unphased allele ordering", {
  input <- scientific_vcf(
    c("a", "b"),
    "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/0/1/1\t1/1/0/0"
  )
  out <- withr::local_tempfile(fileext = ".tsv")
  export_genotypes(input, out, "dosage")
  expect_equal(read.delim(out)$rs, c(2, 2))
  expect_error(export_genotypes(input, tempfile(), "plink"),
    class = "beagler_error_export_polyploid_refused"
  )
  expect_error(
    export_genotypes(input, tempfile(), "plink", downcode = "first_two"),
    class = "beagler_error_export_invalid_downcode"
  )
})

test_that("dosage export distinguishes expected DS from GT counts", {
  input <- scientific_vcf(
    "s",
    "1\t10\trs\tA\tG\t.\tPASS\t.\tGT:DS\t0|1:0.8",
    "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage\">"
  )
  out <- withr::local_tempfile(fileext = ".tsv")
  report <- export_genotypes(input, out, "dosage", field = "DS")
  expect_equal(read.delim(out)$rs, 0.8)
  expect_equal(report$summary$field, "DS")
})

test_that("HapMap exports round-trip through the public converter", {
  input <- scientific_vcf(
    c("a", "b", "c"),
    "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0|1\t1/1\t0/."
  )
  out <- withr::local_tempfile(fileext = ".hmp.txt")
  back <- withr::local_tempfile(fileext = ".vcf")
  report <- export_genotypes(input, out, "hapmap")
  expect_equal(unname(unlist(read.delim(out)[1, 12:14])), c("R", "G", "N"))
  expect_true(hapmap_to_vcf(out, back)$ok)
  expect_true("export_phase_dropped" %in% report$warnings$code)
  expect_true("export_partial_calls_missing" %in% report$warnings$code)
})

test_that("PED export masks partial calls and refuses multiallelic sites", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/.")
  prefix <- tempfile()
  withr::defer(unlink(paste0(prefix, c(".ped", ".map"))))
  report <- export_genotypes(input, prefix, "plink")
  pair <- unname(unlist(read.table(paste0(prefix, ".ped"))[1, 7:8]))
  expect_equal(pair, c(0, 0))
  expect_true("export_partial_calls_missing" %in% report$warnings$code)
  multi <- scientific_vcf("s", "1\t10\trs\tA\tG,T\t.\tPASS\t.\tGT\t1/2")
  expect_error(export_genotypes(multi, tempfile(), "plink"),
    class = "beagler_error_export_multiallelic"
  )
})

test_that("chromosome names report labels without inferring organism ploidy", {
  input <- scientific_vcf("s", c(
    "1A\t10\ta\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "1B\t10\tb\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "1D\t10\td\tA\tG\t.\tPASS\t.\tGT\t0/1"
  ))
  expect_equal(detect_subgenome_labels(input), c("A", "B", "D"))
  expect_equal(detect_genotype_ploidy(input), 2L)
})

test_that("sample merge aligns repeated coordinates by allele identity", {
  a <- scientific_vcf("a", c(
    "1\t10\tG_site\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "1\t10\tT_site\tA\tT\t.\tPASS\t.\tGT\t1/1"
  ))
  b <- scientific_vcf("b", c(
    "1\t10\tT_site\tA\tT\t.\tPASS\t.\tGT\t0/0",
    "1\t10\tG_site\tA\tG\t.\tPASS\t.\tGT\t1/1"
  ))
  out <- withr::local_tempfile(fileext = ".vcf")
  merge_vcf(c(a, b), out)
  parsed <- read_genotype_vcf(out)
  expect_equal(unname(parsed$gt[parsed$fixed$ID == "G_site", "b"]), "1/1")
  expect_equal(unname(parsed$gt[parsed$fixed$ID == "T_site", "b"]), "0/0")
})

test_that("sample merge recomputes counts and reports discarded annotations", {
  meta <- c(
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Count\">",
    "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">",
    "##FILTER=<ID=q10,Description=\"Low quality\">"
  )
  a <- scientific_vcf(
    "a", "1\t10\trs\tA\tG\t50\tPASS\tAC=1;AN=2;AF=0.5\tGT\t0/1", meta
  )
  b <- scientific_vcf(
    "b", "1\t10\trs\tA\tT\t5\tq10\tAC=2;AN=2;AF=1\tGT\t1/1", meta
  )
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- merge_vcf(c(a, b), out)
  parsed <- read_genotype_vcf(out)
  expect_equal(parsed$fixed$INFO, "AC=1,2;AN=4")
  expect_equal(parsed$fixed$FILTER, "q10")
  expect_equal(parsed$fixed$QUAL, ".")
  expect_true("AF" %in% result$summary$dropped_info_fields)
})

test_that("vertical merges preserve distinct overlapping SNPs and indels", {
  a <- scientific_vcf("s", "1\t10\tsnp\tA\tG\t.\tPASS\t.\tGT\t0/1")
  b <- scientific_vcf("s", "1\t10\tdel\tAT\tA\t.\tPASS\t.\tGT\t0/1")
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- merge_vcf(c(a, b), out, "by_marker")
  expect_equal(result$summary$sites_out, 2L)
  expect_equal(nrow(result$summary$sites_dropped), 0L)
})

test_that("merge refuses incompatible reference dictionaries", {
  a <- scientific_vcf(
    "a", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "##contig=<ID=1,length=100,assembly=build1>"
  )
  b <- scientific_vcf(
    "b", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "##contig=<ID=1,length=200,assembly=build2>"
  )
  expect_error(merge_vcf(c(a, b), tempfile(fileext = ".vcf")),
    class = "beagler_error_merge_header_conflict"
  )
})

test_that("vertical merge keeps independent phase blocks independent", {
  meta <- "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">"
  a <- scientific_vcf("s", "1\t10\ta\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:1", meta)
  b <- scientific_vcf("s", "1\t20\tb\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:1", meta)
  out <- withr::local_tempfile(fileext = ".vcf")
  merge_vcf(c(a, b), out, "by_marker")
  parsed <- read_genotype_vcf(out)
  expect_length(unique(parsed$format_fields$PS[, 1]), 2L)
})

test_that("reordered duplicate alleles align by unique marker IDs", {
  a <- scientific_vcf("a", c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "1\t10\ttwo\tA\tG\t.\tPASS\t.\tGT\t1/1"
  ))
  b <- scientific_vcf("b", c(
    "1\t10\ttwo\tA\tG\t.\tPASS\t.\tGT\t0/0",
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t1/1"
  ))
  out <- withr::local_tempfile(fileext = ".vcf")
  merge_vcf(c(a, b), out)
  expect_equal(unname(read_genotype_vcf(out)$gt[, "b"]), c("1/1", "0/0"))
})

test_that("ambiguous repeated-coordinate alignments fail without writing", {
  a <- scientific_vcf("a", c(
    "1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "1\t10\t.\tA\tT\t.\tPASS\t.\tGT\t1/1"
  ))
  b <- scientific_vcf("b", c(
    "1\t10\t.\tA\tG,T\t.\tPASS\t.\tGT\t0/0",
    "1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t1/1"
  ))
  out <- withr::local_tempfile(fileext = ".vcf")
  expect_error(merge_vcf(c(a, b), out),
    class = "beagler_error_merge_ambiguous_markers"
  )
  expect_false(file.exists(out))
})

test_that("empty merge reports dropped records and preserves old output", {
  a <- scientific_vcf("a", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t2/2")
  b <- scientific_vcf("b", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/0")
  out <- withr::local_tempfile(fileext = ".vcf")
  writeLines("existing output", out)
  result <- merge_vcf(c(a, b), out, overwrite = TRUE)
  expect_false(result$ok)
  expect_equal(result$errors$code, "merge_no_records_left")
  expect_equal(readLines(out), "existing output")
  expect_equal(result$summary$sites_dropped$record_index, c(1L, 1L))
  expect_setequal(
    result$summary$sites_dropped$input_path, normalizePath(c(a, b))
  )
  expect_equal(result$summary$sites_dropped$id, c("one", "one"))
})

test_that("coalescing unphased equivalent calls counts each sample once", {
  a <- scientific_vcf("s", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/1")
  b <- scientific_vcf("s", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t1/0")
  out <- withr::local_tempfile(fileext = ".vcf")
  merge_vcf(c(a, b), out, sample_collision = "coalesce")
  expect_equal(read_genotype_vcf(out)$fixed$INFO, "AC=1;AN=2")
  phased <- scientific_vcf("s", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0|1")
  expect_error(
    merge_vcf(c(a, phased), tempfile(fileext = ".vcf"),
              sample_collision = "coalesce"),
    class = "beagler_error_merge_phase_conflict"
  )
})

test_that("implicit phase blocks stay separate across files", {
  a <- scientific_vcf("s", c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0|1",
    "1\t20\ttwo\tA\tG\t.\tPASS\t.\tGT\t1|0"
  ))
  b <- scientific_vcf("s", "1\t30\tthree\tA\tG\t.\tPASS\t.\tGT\t0|1")
  out <- withr::local_tempfile(fileext = ".vcf")
  merge_vcf(c(a, b), out, "by_marker")
  ps <- read_genotype_vcf(out)$format_fields$PS[, 1]
  expect_equal(ps[1], ps[2])
  expect_false(ps[2] == ps[3])
})

test_that("DS requests reject missing or invalid dosage instead of using GT", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  out <- withr::local_tempfile(fileext = ".tsv")
  expect_error(export_genotypes(input, out, "dosage", field = "DS"),
    class = "beagler_error_export_missing_ds"
  )
  for (ds in c("bad", "Inf", "-0.1", "2.1", "0.1,0.9")) {
    input <- scientific_vcf("s", paste0(
      "1\t10\trs\tA\tG\t.\tPASS\t.\tGT:DS\t0/1:", ds
    ))
    expect_error(export_genotypes(input, out, "dosage", field = "DS"),
      class = "beagler_error_export_invalid_ds"
    )
  }
  expect_false(file.exists(out))
})

test_that("VCF overwrite removes stale indexes only after successful writing", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  dir <- withr::local_tempdir()
  out <- file.path(dir, "calls.vcf")
  writeLines("original", out)
  writeLines("old index", paste0(out, ".tbi"))
  expect_error(export_genotypes(input, out, "vcf"),
    class = "beagler_error_export_output_exists"
  )
  expect_equal(readLines(out), "original")
  expect_true(file.exists(paste0(out, ".tbi")))
  export_genotypes(input, out, "vcf", overwrite = TRUE)
  expect_false(file.exists(paste0(out, ".tbi")))
  expect_equal(length(list.files(dir)), 1L)
  expect_equal(nrow(read_genotype_vcf(out)$fixed), 1L)
})

test_that("compressed VCF export never silently substitutes ordinary gzip", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  out <- withr::local_tempfile(fileext = ".vcf.gz")
  writeLines("original", out)
  withr::local_envvar(PATH = "")
  expect_error(export_genotypes(input, out, "vcf", overwrite = TRUE),
    class = "beagler_error_bgzip_not_found"
  )
  expect_equal(readLines(out), "original")
})

test_that("empty HapMap export writes a header without inventing marker rows", {
  input <- scientific_vcf("s", character())
  out <- withr::local_tempfile(fileext = ".hmp.txt")
  result <- export_genotypes(input, out, "hapmap")
  expect_equal(result$summary$sites_out, 0L)
  expect_equal(nrow(read.delim(out)), 0L)
})

test_that("VCF export retains per-record FORMAT fields while sorting", {
  input <- scientific_vcf("s", c(
    "1\t20\ttwo\tA\tG\t.\tPASS\t.\tGT:DS:PS\t0|1:0.8:20",
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:10"
  ), c(
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">",
    "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage, expected\">"
  ))
  out <- withr::local_tempfile(fileext = ".vcf")
  export_genotypes(input, out, "vcf")
  parsed <- read_genotype_vcf(out)
  expect_equal(parsed$formats, c("GT:PS", "GT:DS:PS"))
  expect_equal(parsed$format_fields$PS[, 1], c("10", "20"))
  expect_equal(parsed$format_fields$DS[, 1], c(".", "0.8"))
})

test_that("merge rejects conflicting typed metadata", {
  a <- scientific_vcf(
    "a", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT:DS\t0/1:1",
    "##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Dosage, expected\">"
  )
  b <- scientific_vcf(
    "b", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT:DS\t0/1:1",
    "##FORMAT=<ID=DS,Number=R,Type=Float,Description=\"Dosage\">"
  )
  expect_error(merge_vcf(c(a, b), tempfile(fileext = ".vcf")),
    class = "beagler_error_merge_header_conflict"
  )
})

test_that("marker merge permits an empty input after reported filtering", {
  a <- scientific_vcf("s", "1\t10\tbad\tA\tG\t.\tPASS\t.\tGT\t2/2")
  b <- scientific_vcf("s", "1\t20\tgood\tA\tG\t.\tPASS\t.\tGT\t0/1")
  out <- withr::local_tempfile(fileext = ".vcf")
  report <- merge_vcf(c(a, b), out, "by_marker")
  expect_equal(read_genotype_vcf(out)$fixed$ID, "good")
  expect_equal(report$summary$sites_dropped$id, "bad")
})

test_that("malformed FORMAT values are not silently discarded", {
  for (record in c("GT:DS:DS\t0/1:1:2", "GT:DS\t0/1:1:2", "GT:DS\t0/1:")) {
    input <- scientific_vcf("s", paste0(
      "1\t10\trs\tA\tG\t.\tPASS\t.\t", record
    ))
    expect_error(export_genotypes(input, tempfile(fileext = ".vcf"), "vcf"),
      class = "beagler_error_genotype_vcf_invalid"
    )
  }
})
