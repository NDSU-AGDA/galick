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
  resolved <- tryCatch(find_plink(), galick_error_plink_not_found = function(e) NULL)
  if (is.null(resolved)) {
    testthat::skip("genetics PLINK not found on PATH (PuTTY plink is not supported)")
  }
  invisible(resolved)
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
