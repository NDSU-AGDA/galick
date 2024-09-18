use extendr_api::prelude::*;
use std::process::Command;

/// Run Beagle with a specified file.
///
/// This function uses Beagle software to perform genotype phasing and/or imputation.
/// Beagle is a tool for phasing genotypes, estimating haplotypes, and imputing ungenotyped markers.
///
/// # Parameters
/// - `gt`: The path to the input genotype file (VCF format).
/// - `out`: The prefix for the output file.
/// - `n_threads`: (Optional) Specifies the number of computational threads that will be used (default: number of CPU cores).
/// - `burnin`: (Optional) Maximum number of burnin iterations to estimate the initial haplotype frequency model for genotype phasing (default: 3).
/// - `iterations`: (Optional) Number of iterations used to estimate genotype phase (default: 12). Higher values increase computation time but improve accuracy.
/// - `phase_states`: (Optional) Number of model states used to estimate genotype phase (default: 280).
/// - `impute`: (Optional) Whether to impute markers present in the reference panel but absent in the target (default: true).
///
/// # Beagle Phasing Parameters
/// - `burnin`: Number of burnin iterations (default: 3).
/// - `iterations`: Number of iterations for phasing (default: 12).
/// - `phase_states`: Number of model states for phasing (default: 280).
///
/// # Beagle Imputation Parameters
/// - `impute`: If true, markers from the reference panel missing in the target will be imputed (default: true).
///
/// # Examples
/// ```r
/// run_beagle("input.vcf", "output", Some(4), Some(5), Some(15), Some(300), Some(true))
/// ```
///
/// The above example will run Beagle using 4 threads, 5 burnin iterations, 15 phasing iterations, 300 phase states, and with imputation enabled.
///
/// # Errors
/// This function will panic if the `BEAGLE_PATH` environment variable is not defined or if the Beagle command fails to execute.
#[doc(hidden)]
#[extendr]
fn run_beagle_internal(
    gt: &str,
    out: &str,
    n_threads: Option<u32>,
    burnin: Option<u32>,
    iterations: Option<u32>,
    phase_states: Option<u32>,
    impute: Option<bool>,
) {
    let beagle_path = std::env::var("BEAGLE_PATH").expect("No variable BEAGLE_PATH defined");

    // Unwrap the options, using default values if None
    let n_threads = n_threads.unwrap_or(1); // Default to 1 thread
    let burnin = burnin.unwrap_or(3); // Default burnin is 3
    let iterations = iterations.unwrap_or(12); // Default iterations is 12
    let phase_states = phase_states.unwrap_or(280); // Default phase states is 280
    let impute = impute.unwrap_or(true); // Default impute is true

    let _ = Command::new("java")
        .arg("-jar")
        .arg(beagle_path)
        .arg(format!("gt={gt}"))
        .arg(format!("out={out}"))
        .arg(format!("nthreads={n_threads}"))
        .arg(format!("burnin={burnin}"))
        .arg(format!("iterations={iterations}"))
        .arg(format!("phase-states={phase_states}"))
        .arg(format!("impute={}", impute))
        .status()
        .expect("Failed to run beagle");
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod beagler;
    fn run_beagle_internal;
}
