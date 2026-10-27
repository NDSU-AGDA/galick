test_that("parse_java_version_output handles Java 8 and Java 9+ formats", {
  java8 <- 'java version "1.8.0_281"\nJava(TM) SE Runtime Environment'
  v <- galick:::parse_java_version_output(java8)
  expect_equal(v$major, 8L)
  expect_equal(v$raw, "1.8.0_281")

  java11 <- 'openjdk version "11.0.21" 2023-10-17\nOpenJDK Runtime Environment'
  v <- galick:::parse_java_version_output(java11)
  expect_equal(v$major, 11L)
  expect_equal(v$raw, "11.0.21")

  java17 <- 'openjdk version "17" 2021-09-14\nOpenJDK Runtime Environment'
  v <- galick:::parse_java_version_output(java17)
  expect_equal(v$major, 17L)
})

test_that("parse_java_version_output returns NULL on garbage input", {
  expect_null(galick:::parse_java_version_output(""))
  expect_null(galick:::parse_java_version_output("hello there, no version here"))
  expect_null(galick:::parse_java_version_output(character()))
})

test_that("parse_beagle_version_output extracts (version X.Y) banner", {
  banner <- "beagle.29Oct24.c8e.jar  (version 5.4)\nUsage: ..."
  v <- galick:::parse_beagle_version_output(banner)
  expect_equal(v$major, 5L)
  expect_equal(v$minor, 4L)
  expect_equal(v$raw, "5.4")
})

test_that("parse_beagle_version_output returns NULL when no version line is present", {
  expect_null(galick:::parse_beagle_version_output("nothing version-shaped here"))
  expect_null(galick:::parse_beagle_version_output(""))
})

test_that("find_java raises a classed error when java is unreachable", {
  withr::with_envvar(
    new = c(JAVA_HOME = "", PATH = "/nonexistent_xyz_path"),
    expect_error(
      galick:::find_java(),
      class = "galick_error_java_not_found"
    )
  )
})

test_that("find_java honours an explicit valid path", {
  jpath <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  resolved <- galick:::find_java(jpath)
  expect_match(resolved, "Rscript")
})

test_that("find_java raises when explicit path does not exist", {
  expect_error(
    galick:::find_java("/nonexistent/path/to/java"),
    class = "galick_error_java_not_found"
  )
})

test_that("find_beagle_jar honours BEAGLE_PATH when set and the file exists", {
  fake <- withr::local_tempfile(fileext = ".jar")
  file.create(fake)
  withr::with_envvar(
    new = c(BEAGLE_PATH = fake),
    {
      resolved <- galick:::find_beagle_jar()
      expect_match(resolved, "\\.jar$")
    }
  )
})

test_that("find_beagle_jar raises when nothing is found", {
  withr::with_envvar(
    new = c(BEAGLE_PATH = ""),
    {
      cache_path <- file.path(tools::R_user_dir("galick", "cache"), "beagle.jar")
      withr::local_options(list(galick.skip_cache_check = TRUE))
      if (!file.exists(cache_path)) {
        expect_error(
          galick:::find_beagle_jar(),
          class = "galick_error_beagle_not_found"
        )
      } else {
        skip("Beagle JAR present in cache; cannot test 'not found' path here.")
      }
    }
  )
})

test_that("find_beagle_jar raises when explicit path does not exist", {
  expect_error(
    galick:::find_beagle_jar("/nonexistent/beagle.jar"),
    class = "galick_error_beagle_not_found"
  )
})

test_that("detect_java_version returns version >= 8 on the live system java", {
  skip_if_no_java()
  java <- galick:::find_java()
  v <- galick:::detect_java_version(java)
  expect_gte(v$major, 8L)
})

test_that("find_bcftools and find_tabix return path or NULL without raising", {
  expect_true(is.null(galick:::find_bcftools()) || is.character(galick:::find_bcftools()))
  expect_true(is.null(galick:::find_tabix()) || is.character(galick:::find_tabix()))
})
