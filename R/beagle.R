#' Run Beagle for Phasing and Imputation
#'
#' This function uses Beagle software to perform genotype phasing and/or imputation.
#' Beagle is a tool for phasing genotypes, estimating haplotypes, and imputing ungenotyped markers.
#'
#' This function accepts genotype data in VCF, CSV, or TXT format. If a CSV or TXT file is provided,
#' it will be converted to VCF format before running Beagle.
#'
#' @param gt The path to the input genotype file (VCF, CSV, or TXT format).
#' @param out The prefix for the output file.
#' @param n_threads (Optional) Specifies the number of computational threads that will be used (default: 1).
#' @param burnin (Optional) Maximum number of burnin iterations to estimate the initial haplotype frequency model for genotype phasing (default: 3).
#' @param iterations (Optional) Number of iterations used to estimate genotype phase (default: 12). Higher values increase computation time but improve accuracy.
#' @param phase_states (Optional) Number of model states used to estimate genotype phase (default: 280).
#' @param impute (Optional) Whether to impute markers present in the reference panel but absent in the target (default: \code{TRUE}).
#'
#' @details
#' **Beagle Phasing Parameters**
#' \itemize{
#'   \item \code{burnin}: Number of burnin iterations (default: 3).
#'   \item \code{iterations}: Number of iterations for phasing (default: 12).
#'   \item \code{phase_states}: Number of model states for phasing (default: 280).
#' }
#'
#' **Beagle Imputation Parameters**
#' \itemize{
#'   \item \code{impute}: If \code{TRUE}, markers from the reference panel missing in the target will be imputed (default: \code{TRUE}).
#' }
#'
#' @details
#' For detailed information about the Beagle algorithms, refer to the original publications:
#'
#' - **Phasing**: Browning BL et al. (2021). Fast two-stage phasing of large-scale sequence data. *Am J Hum Genet* 108(10):1880-1890. \doi{10.1016/j.ajhg.2021.08.005}
#' - **Imputation**: Browning BL et al. (2018). A one-penny imputed genome from next generation reference panels. *Am J Hum Genet* 103(3):338-348. \doi{10.1016/j.ajhg.2018.07.015}
#'
#'
#' @examples
#' \dontrun{
#' run_beagle(
#'   gt = "input.vcf",
#'   out = "output",
#'   n_threads = 4,
#'   burnin = 5,
#'   iterations = 15,
#'   phase_states = 300,
#'   impute = TRUE
#' )
#' }
#'
#'
#' @export
run_beagle <- function(gt, out,
                       n_threads = 1,
                       burnin = 3,
                       iterations = 12,
                       phase_states = 280,
                       impute = TRUE) {
 # Call the Rust function, passing NULL for arguments not specified
 run_beagle_internal(
  gt = gt,
  out = out,
  n_threads = if (missing(n_threads)) NULL else n_threads,
  burnin = if (missing(burnin)) NULL else burnin,
  iterations = if (missing(iterations)) NULL else iterations,
  phase_states = if (missing(phase_states)) NULL else phase_states,
  impute = if (missing(impute)) NULL else impute
 )

 # After Beagle has run, decompress the output VCF file but keep the original .gz file
 output_vcf_gz <- paste0(out, ".vcf.gz")
 output_vcf <- paste0(out, ".vcf")

 if (file.exists(output_vcf_gz)) {
  R.utils::gunzip(
   output_vcf_gz,
   destname = output_vcf,
   remove = FALSE,
   overwrite = TRUE
  )
  message("VCF file decompressed successfully to: ", output_vcf)
 } else {
  warning("The .vcf.gz file does not exist: ", output_vcf_gz)
 }
}

