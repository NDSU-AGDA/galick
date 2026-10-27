test_that("comparison refuses allele swaps and ambiguous variant identities", {
  truth <- scientific_vcf(
    c("a", "b", "c"),
    "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
  )
  swapped <- scientific_vcf(
    c("a", "b", "c"),
    "1\t10\trs\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
  )
  result <- compare_imputed_vs_true(swapped, truth)
  expect_false(result$ok)
  expect_true("compare_imputed_variant_mismatch" %in% result$errors$code)
  expect_true(is.na(result$summary$concordance_rate))
})

test_that("expected dosage is scored explicitly and constant variance is defined", {
  truth <- scientific_vcf(
    c("a", "b", "c"),
    "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
  )
  predicted <- scientific_vcf(
    c("a", "b", "c"),
    "1\t10\trs\tA\tG\t.\tPASS\t.\tGT:DS\t0/0:1\t0/1:1\t1/1:1"
  )
  result <- compare_imputed_vs_true(predicted, truth, field = "DS")
  expect_true(result$ok)
  expect_equal(result$summary$field, "DS")
  expect_equal(result$summary$rmse, sqrt(2 / 3))
  expect_true(is.na(result$summary$dosage_r2))
  expect_equal(result$summary$n_undefined_r2, 1L)
  absent <- compare_imputed_vs_true(truth, truth, field = "DS")
  expect_false(absent$ok)
  expect_true("compare_imputed_missing_ds" %in% absent$errors$code)
})

test_that("missing calls and dropped variants remain in comparison denominators", {
  truth <- scientific_vcf(c("a", "b", "c"), c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "1\t20\ttwo\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1"
  ))
  predicted <- scientific_vcf(
    c("a", "b", "c"),
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t./."
  )
  result <- compare_imputed_vs_true(predicted, truth)
  expect_true(result$ok)
  expect_equal(result$summary$n_eligible, 6L)
  expect_equal(result$summary$n_compared, 2L)
  expect_equal(result$summary$n_missing_predictions, 4L)
  expect_equal(result$summary$completion_rate, 1 / 3)
  expect_equal(result$summary$concordance_rate, 1 / 3)
  expect_equal(result$summary$conditional_concordance_rate, 1)
  expect_equal(nrow(result$summary$per_marker), 2L)
})

test_that("singleton comparison returns an undefined correlation without error", {
  truth <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  result <- compare_imputed_vs_true(truth, truth)
  expect_true(result$ok)
  expect_equal(result$summary$concordance_rate, 1)
  expect_true(is.na(result$summary$dosage_r2))
})

test_that("a missing mask record is not interpreted as an intentionally masked call", {
  truth <- scientific_vcf("s", c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "1\t20\ttwo\tC\tT\t.\tPASS\t.\tGT\t0/1"
  ))
  mask <- scientific_vcf("s", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t./.")
  result <- compare_imputed_vs_true(truth, truth, mask)
  expect_false(result$ok)
  expect_true("compare_imputed_incomplete_mask" %in% result$errors$code)
})

test_that("VCF validation bounds GT indices by the record allele count", {
  invalid <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/2")
  expect_true("vcf_invalid_gt" %in% validate_vcf(invalid)$errors$code)
  expect_false(validate_vcf(invalid, target = "beagle")$ok)
  out <- withr::local_tempfile(fileext = ".vcf")
  repaired <- repair_vcf(invalid, out)
  expect_false(repaired$ok)
  expect_false(file.exists(out))
})

test_that("chromosome merging refuses failed jobs and missing files", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  good <- new_galick_run(status = 0L, output_paths = input)
  bad <- new_galick_run(status = 1L, output_paths = tempfile(fileext = ".vcf"))
  out <- withr::local_tempfile(fileext = ".vcf")
  expect_error(merge_beagle_outputs(list(good, bad), out),
    class = "galick_error_merge_beagle_invalid_input"
  )
  bad$status <- 0L
  expect_error(merge_beagle_outputs(list(good, bad), out),
    class = "galick_error_merge_beagle_invalid_input"
  )
  expect_false(file.exists(out))
})

test_that("chromosome merging aligns samples and retains all FORMAT values", {
  meta <- c(
    '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Dosage">',
    '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Probabilities">'
  )
  a <- scientific_vcf(
    c("a", "b"),
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT:DS:GP\t0/0:0.1:0.9,0.1,0\t1/1:1.9:0,0.1,0.9", meta
  )
  b <- scientific_vcf(
    c("b", "a"),
    "2\t10\ttwo\tC\tT\t.\tPASS\t.\tGT:DS:GP\t0/1:1:0,1,0\t0/0:0:1,0,0", meta
  )
  runs <- lapply(c(a, b), function(path) new_galick_run(status = 0L, output_paths = path))
  out <- withr::local_tempfile(fileext = ".vcf")
  merge_beagle_outputs(runs, out)
  parsed <- read_genotype_vcf(out)
  expect_equal(unname(parsed$gt[2, ]), c("0/0", "0/1"))
  expect_equal(unname(parsed$format_fields$DS[2, ]), c("0", "1"))
  expect_equal(unname(parsed$format_fields$GP[2, ]), c("1,0,0", "0,1,0"))
})

test_that("reference preflight detects assembly and duplicate-position ambiguity", {
  a <- scientific_vcf(
    "target", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "##reference=assembly1"
  )
  b <- scientific_vcf(
    "reference", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0|1",
    "##reference=assembly2"
  )
  expect_false(validate_reference_match(a, b)$ok)
  duplicates <- scientific_vcf("reference", c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0|1",
    "1\t10\ttwo\tA\tT\t.\tPASS\t.\tGT\t0|1"
  ))
  expect_false(validate_reference_match(a, duplicates)$ok)
})

test_that("reference readiness checks phased complete calls and sample overlap", {
  a <- scientific_vcf("target", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  unphased <- scientific_vcf("reference", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  missing <- scientific_vcf("reference", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0|.")
  expect_false(validate_reference_match(a, unphased, require_phased = TRUE)$ok)
  expect_false(validate_reference_match(a, missing, require_phased = TRUE)$ok)
  overlap <- scientific_vcf("target", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0|1")
  expect_false(validate_reference_match(a, overlap, require_phased = TRUE)$ok)
})

test_that("file hashing is portable and checksum failure cannot resume", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  expect_equal(file_sha256(input), unname(digest::digest(file = input, algo = "sha256")))
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "imputed")
  output <- paste0(prefix, ".vcf.gz")
  file.copy(input, output)
  jsonlite::write_json(
    list(
      status = 0L,
      args = list(job_fingerprint = "expected"),
      inputs = list(list(path = input, sha256 = file_sha256(input))),
      outputs = list(list(path = output, sha256 = "wrong"))
    ),
    file.path(dir, "run_record.json"),
    auto_unbox = TRUE
  )
  expect_null(try_resume_chrom_run(input, dir, prefix, "expected"))
})

test_that("job fingerprints cover inputs, parameters, executables and versions", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  dir <- withr::local_tempdir()
  jar <- file.path(dir, "beagle.jar")
  java <- file.path(dir, "java")
  map <- file.path(dir, "map.txt")
  writeLines("jar1", jar)
  writeLines("java1", java)
  writeLines("map1", map)
  local_mocked_bindings(
    find_java = function(...) java,
    find_beagle_jar = function(...) jar,
    detect_java_version = function(...) list(raw = "java 17"),
    detect_beagle_version = function(...) list(raw = "beagle 5.4")
  )
  base <- list(gt = input, out = file.path(dir, "imputed"), map = map, seed = 1L)
  fingerprint <- function(args) do.call(prepare_beagle_job, args)$fingerprint
  original <- fingerprint(base)
  for (change in list(
    list(seed = 2L), list(nthreads = 2L), list(memory_gb = 8L),
    list(mode = "phase"), list(extra_args = "iterations=20"),
    list(ref = input)
  )) {
    expect_false(identical(original, fingerprint(utils::modifyList(base, change))))
  }
  for (path in c(input, map, jar, java)) {
    old <- readLines(path)
    writeLines(c(old, "changed"), path)
    expect_false(identical(original, fingerprint(base)))
    writeLines(old, path)
  }
  expect_equal(original, fingerprint(base))
  expect_error(do.call(prepare_beagle_job, c(base, list(extra_args = "gt=elsewhere"))),
    class = "galick_error_run_beagle_invalid_input"
  )
})

test_that("resume requires matching configuration, checksums and valid output", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "imputed")
  output <- paste0(prefix, ".vcf.gz")
  write_output <- function(lines) {
    con <- gzfile(output, "wt")
    writeLines(lines, con)
    close(con)
  }
  write_output(readLines(input))
  write_record <- function(checksum = file_sha256(output)) {
    write_run_record_json(
      dir,
      command = c("java", "-jar", "beagle.jar"),
      versions = list(galick = "0.1.0", r = "4.1.0", os = "test"),
      args = list(job_fingerprint = "expected"),
      inputs = list(list(path = input, sha256 = file_sha256(input))),
      outputs = list(list(path = output, sha256 = checksum)),
      runtime_sec = 1, status = 0L, started_at = "2026-10-21T00:00:00Z",
      ended_at = "2026-10-21T00:00:01Z"
    )
  }
  write_record()
  expect_equal(try_resume_chrom_run(input, dir, prefix, "expected")$status, 0L)
  expect_null(try_resume_chrom_run(input, dir, prefix, "changed_configuration"))
  write_record(strrep("0", 64L))
  expect_null(try_resume_chrom_run(input, dir, prefix, "expected"))
  write_output("not a VCF")
  write_record()
  expect_null(try_resume_chrom_run(input, dir, prefix, "expected"))
  writeLines("{}", file.path(dir, "run_record.json"))
  expect_null(try_resume_chrom_run(input, dir, prefix, "expected"))
})

test_that("masking and benchmarking preserve eligibility and random seeds", {
  input <- scientific_vcf(c("a", "b", "c"), c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t.\t0/.\t0/1",
    "1\t20\ttwo\tC\tT\t.\tPASS\t.\tGT\t0/0\t1/1\t0/1"
  ))
  out <- withr::local_tempfile(fileext = ".vcf")
  result <- mask_genotypes(input, out, fraction = 0.5, seed = 1L)
  expect_equal(result$summary$n_eligible, 4L)
})

test_that("a failed overwrite cannot reuse or destroy the previous output", {
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  dir <- withr::local_tempdir()
  prefix <- file.path(dir, "imputed")
  output <- paste0(prefix, ".vcf.gz")
  con <- gzfile(output, "wt")
  writeLines(readLines(input), con)
  close(con)
  writeLines("previous index", paste0(output, ".tbi"))
  previous <- file_sha256(output)
  local_mocked_bindings(
    prepare_beagle_job = function(...) {
      list(
        java = "java", argv = c("-jar", "beagle.jar"),
        versions = list(galick = "0.1.0", r = "4.1.0", os = "test"),
        files = list(gt = list(path = input, sha256 = file_sha256(input))),
        fingerprint = "test"
      )
    },
    run_cli = function(cmd, args, stdout_log, stderr_log) {
      writeLines("", stdout_log)
      writeLines("", stderr_log)
      list(status = 0L, runtime_sec = 1)
    }
  )
  result <- run_beagle(input, prefix, run_dir = dir, overwrite = TRUE)
  expect_equal(result$status, 1L)
  expect_length(result$output_paths, 0L)
  expect_equal(file_sha256(output), previous)
  expect_equal(readLines(paste0(output, ".tbi")), "previous index")
  expect_equal(jsonlite::fromJSON(result$run_record)$status, 1L)
})

test_that("chromosome concatenation rejects sample and header conflicts", {
  one <- scientific_vcf(
    "a", "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "##reference=build_one"
  )
  two <- scientific_vcf(
    "b", "2\t10\ttwo\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "##reference=build_two"
  )
  runs <- lapply(c(one, two), function(path) new_galick_run(status = 0L, output_paths = path))
  out <- withr::local_tempfile(fileext = ".vcf")
  expect_error(merge_beagle_outputs(runs, out),
    class = "galick_error_merge_beagle_invalid_input"
  )
  expect_false(file.exists(out))
  two <- scientific_vcf(
    "a", "2\t10\ttwo\tA\tG\t.\tPASS\t.\tGT\t0/1",
    "##reference=build_two"
  )
  runs[[2]]$output_paths <- two
  expect_error(merge_beagle_outputs(runs, out))
  expect_false(file.exists(out))
  attr(runs, "expected_chromosomes") <- c("1", "2", "3")
  expect_error(merge_beagle_outputs(runs, out),
    class = "galick_error_merge_beagle_invalid_input"
  )
})

test_that("compressed chromosome output is indexable BGZF", {
  skip_if(
    Sys.which("bgzip") == "" || Sys.which("tabix") == "",
    "bgzip and tabix required"
  )
  input <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/1")
  run <- new_galick_run(status = 0L, output_paths = input)
  output <- withr::local_tempfile(fileext = ".vcf.gz")
  merge_beagle_outputs(list(run), output)
  expect_equal(run_cli(Sys.which("tabix"), c("-p", "vcf", output))$status, 0L)
  queried <- run_cli(Sys.which("tabix"), c(output, "1:10-10"))
  expect_equal(queried$status, 0L)
  expect_match(queried$stdout, "rs")
})

test_that("Beagle accepts a plain phased reference without an index", {
  skip_if_no_java()
  skip_if_no_beagle()
  target <- scientific_vcf(c("t1", "t2", "t3"), c(
    "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0",
    "1\t300\trs3\tA\tT\t.\tPASS\t.\tGT\t1/1\t0/0\t0/1"
  ), "##contig=<ID=1>")
  reference <- scientific_vcf(c("r1", "r2", "r3"), c(
    "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
    "1\t200\trs2\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0",
    "1\t300\trs3\tA\tT\t.\tPASS\t.\tGT\t1|1\t0|0\t0|1"
  ), "##contig=<ID=1>")
  dir <- withr::local_tempdir()
  result <- run_beagle(target, file.path(dir, "imputed"),
    ref = reference,
    seed = 17L, run_dir = dir
  )
  expect_equal(result$status, 0L)
  expect_true(validate_vcf(result$output_paths, target = "beagle")$ok)
})

test_that("extraction and masking share semantic allele bounds checks", {
  invalid <- scientific_vcf("s", "1\t10\trs\tA\tG\t.\tPASS\t.\tGT\t0/2")
  output <- withr::local_tempfile(fileext = ".vcf")
  expect_error(vcf_to_dosage(invalid), class = "galick_error_unresolvable_alleles")
  expect_error(vcf_to_hardcalls(invalid), class = "galick_error_unresolvable_alleles")
  expect_error(mask_genotypes(invalid, output), class = "galick_error_unresolvable_alleles")
  expect_false(file.exists(output))
})

test_that("benchmark forwards the masking seed to Beagle", {
  input <- scientific_vcf(c("a", "b", "c"), c(
    "1\t10\tone\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1",
    "1\t20\ttwo\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\t0/0"
  ))
  observed <- NULL
  local_mocked_bindings(run_beagle = function(gt, out, mode, seed, overwrite, ...) {
    observed <<- seed
    new_galick_run(status = 1L)
  })
  result <- benchmark_imputation(input, withr::local_tempdir(), fraction = 0.5, seed = 42L)
  expect_identical(observed, 42L)
  expect_equal(result$mask$summary$n_masked, 3L)
})
