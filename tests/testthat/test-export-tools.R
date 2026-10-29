test_that("PLINK imports PED exports with partially missing calls", {
  skip_if_no_plink()
  input <- scientific_vcf(c("a", "b", "c"), c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/.\t0/1\t1/1",
    "1\t20\ttwo\tC\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1"
  ))
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "exported")
  export_genotypes(input, prefix, "plink")
  out <- file.path(dir, "checked")
  result <- run_cli(find_plink(), c(
    "--file", prefix, "--make-bed", "--out", out,
    "--threads", "1", "--memory", "128"
  ))
  expect_equal(result$status, 0L, info = result$stderr)
  expect_true(all(file.exists(paste0(out, c(".bed", ".bim", ".fam")))))
  expect_equal(nrow(read.table(paste0(out, ".bim"))), 2L)
  expect_equal(nrow(read.table(paste0(out, ".fam"))), 3L)
})

test_that("BGZF exports are readable by tabix and bcftools", {
  skip_if(Sys.which("bgzip") == "", "bgzip not found on PATH")
  skip_if(Sys.which("tabix") == "", "tabix not found on PATH")
  input <- scientific_vcf("s", c(
    "1\t20\ttwo\tA\tG\t.\tPASS\t.\tGT:PS\t0|1:20",
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:10"
  ), c(
    "##contig=<ID=1,length=100>",
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">"
  ))
  out <- file.path(withr::local_tempdir(), "exported.vcf.gz")
  export_genotypes(input, out, "vcf")
  bytes <- readBin(out, "raw", n = 18L)
  expect_identical(bytes[13:16], as.raw(c(66, 67, 2, 0)))
  indexed <- run_cli(Sys.which("tabix"), c("-p", "vcf", out))
  expect_equal(indexed$status, 0L, info = indexed$stderr)
  query <- run_cli(Sys.which("tabix"), c(out, "1:10-10"))
  expect_equal(query$status, 0L, info = query$stderr)
  expect_match(
    query$stdout, "1\t10\tone\tA\tG\t.\tPASS\t.\tGT:PS\t1|0:10",
    fixed = TRUE
  )
  export_genotypes(input, out, "vcf", overwrite = TRUE)
  expect_false(file.exists(paste0(out, ".tbi")))
  skip_if_no_bcftools()
  checked <- run_cli(Sys.which("bcftools"), c("view", "-Ov", out))
  expect_equal(checked$status, 0L, info = checked$stderr)
  expect_match(checked$stdout, "GT:PS\t1|0:10", fixed = TRUE)
})

test_that("bcftools reads reconciled allele annotations and phase blocks", {
  skip_if_no_bcftools()
  meta <- c(
    "##contig=<ID=1,length=100>",
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">"
  )
  a <- scientific_vcf(
    "a", "1\t10\trs\tA\tG\t.\tPASS\tAF=0.5\tGT:PS:AD\t0|1:1:8,9", meta
  )
  b <- scientific_vcf(
    "b", "1\t10\trs\tA\tT\t.\tPASS\tAF=1\tGT:PS:AD\t1|1:1:0,7", meta
  )
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- merge_vcf(c(a, b), out)
  expect_equal(result$summary$dropped_format_fields, "AD")
  checked <- run_cli(Sys.which("bcftools"), c(
    "query", "-f", "%ALT\t%AC\t%AN[\t%GT:%PS]\n", out
  ))
  expect_equal(checked$status, 0L, info = checked$stderr)
  expect_equal(trimws(checked$stdout), "G,T\t1,2\t4\t0|1:1\t2|2:2")
})
