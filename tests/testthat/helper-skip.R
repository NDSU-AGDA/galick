skip_if_no_java <- function() {
  if (Sys.which("java") == "") {
    testthat::skip("java not found on PATH")
  }
}

skip_if_no_beagle <- function() {
  jar <- Sys.getenv("BEAGLE_PATH", unset = "")
  if (jar == "" || !file.exists(jar)) {
    testthat::skip("BEAGLE_PATH not set or JAR missing")
  }
}

skip_if_no_plink <- function() {
  if (Sys.which("plink") == "") {
    testthat::skip("plink not found on PATH")
  }
}

skip_if_no_bcftools <- function() {
  if (Sys.which("bcftools") == "") {
    testthat::skip("bcftools not found on PATH")
  }
}

skip_if_no_tabix <- function() {
  if (Sys.which("tabix") == "") {
    testthat::skip("tabix not found on PATH")
  }
}
