write_export_fixture <- function(path) {
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr2,length=1000>",
    "##contig=<ID=chr1>",
    "##galick-test=preserve-me",
    paste(
      c(
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
        "INFO", "FORMAT", "s1", "s2"
      ),
      collapse = "\t"
    ),
    "chr1\t100000000\trs2\tA\tG\t9.5\tPASS\tAF=0.0000123456789\tGT\t0|.\t1|0",
    "chr2\t10\trs1\tC\tT\t.\tPASS\t.\tGT\t0/0/1/1\t1|0|.|1"
  ), path)
  path
}

test_that("subgenome labels and genotype ploidy remain separate", {
  path <- withr::local_tempfile(fileext = ".vcf")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    "chr1A\t1\ta\tA\tG\t.\tPASS\t.\tGT\t0/0/0/1/1/1",
    "chr1B\t1\tb\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "chr1D\t1\td\tA\tG\t.\tPASS\t.\tGT\t."
  ), path)

  expect_equal(detect_subgenome_labels(path), c("A", "B", "D"))
  genotype_ploidy <- detect_genotype_ploidy(path)
  expect_equal(genotype_ploidy, 6L)
})

test_that("VCF export preserves GT semantics, metadata, and valid positions", {
  input <- write_export_fixture(withr::local_tempfile(fileext = ".vcf"))
  output <- withr::local_tempfile(fileext = ".vcf")

  report <- export_genotypes(input, output, format = "vcf")
  expect_s3_class(report, "galick_result")
  expect_equal(report$summary$sites_in, 2L)
  expect_equal(report$summary$sites_out, 2L)
  expect_equal(report$summary$genotype_ploidy, 4L)

  lines <- readLines(output)
  expect_identical(lines[1], "##fileformat=VCFv4.2")
  expect_true("##galick-test=preserve-me" %in% lines)
  expect_true("##contig=<ID=chr1>" %in% lines)
  expect_false(any(grepl("1e\\+08", lines, ignore.case = TRUE)))

  parsed <- galick:::read_genotype_vcf(output)
  expect_equal(parsed$fixed$ID, c("rs1", "rs2"))
  expect_equal(
    unname(parsed$gt[parsed$fixed$ID == "rs1", ]),
    c("0/0/1/1", "1|0|.|1")
  )
  expect_equal(
    unname(parsed$gt[parsed$fixed$ID == "rs2", ]),
    c("0|.", "1|0")
  )
})

test_that("diploid formats refuse polyploid GT; counts retain all copies", {
  input <- write_export_fixture(withr::local_tempfile(fileext = ".vcf"))
  prefix <- withr::local_tempfile()
  dosage <- withr::local_tempfile(fileext = ".tsv")

  expect_error(
    export_genotypes(input, prefix, format = "plink"),
    class = "galick_error_export_polyploid_refused"
  )
  expect_error(
    export_genotypes(input, tempfile(fileext = ".hmp.txt"), format = "hapmap"),
    class = "galick_error_export_polyploid_refused"
  )

  report <- export_genotypes(input, dosage, format = "dosage")
  expect_true(file.exists(dosage))
  expect_equal(read.delim(dosage)$rs1, c(2, NA_real_))
  expect_equal(report$summary$field, "GT")
})

test_that("delimited exports preserve allele slots and phase", {
  input <- write_export_fixture(withr::local_tempfile(fileext = ".vcf"))
  csv <- withr::local_tempfile(fileext = ".csv")

  export_genotypes(input, csv, format = "csv")

  table <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
  expect_equal(table$s1, c("0|.", "0/0/1/1"))
  expect_equal(table$s2, c("1|0", "1|0|.|1"))
  expect_equal(table$POS, c(100000000, 10))
})

test_that("exported VCF parses in vcfR when available", {
  skip_if_not_installed("vcfR")
  input <- write_export_fixture(withr::local_tempfile(fileext = ".vcf"))
  output <- withr::local_tempfile(fileext = ".vcf")
  export_genotypes(input, output, format = "vcf")

  parsed <- vcfR::read.vcfR(output, verbose = FALSE)
  expect_equal(nrow(parsed@fix), 2L)
})

test_that("VCF round-trips haploid through hexaploid GT", {
  input <- withr::local_tempfile(fileext = ".vcf")
  vcf_output <- withr::local_tempfile(fileext = ".vcf")
  calls <- c("1", "0|.", "0/0/1/1", "0|1|.|0|1|1")
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1",
    paste0("chr1\t1\tp1\tA\tG\t.\tPASS\t.\tGT\t", calls[1]),
    paste0("chr1\t2\tp2\tC\tT\t.\tPASS\t.\tGT\t", calls[2]),
    paste0("chr1\t3\tp4\tG\tA\t.\tPASS\t.\tGT\t", calls[3]),
    paste0("chr1\t4\tp6\tT\tC\t.\tPASS\t.\tGT\t", calls[4])
  ), input)

  export_genotypes(input, vcf_output, format = "vcf")

  vcf_roundtrip <- galick:::read_genotype_vcf(vcf_output)
  expect_equal(unname(vcf_roundtrip$gt[, 1]), calls)
  expect_equal(detect_genotype_ploidy(vcf_output), 6L)
})
