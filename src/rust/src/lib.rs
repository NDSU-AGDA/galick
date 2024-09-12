use std::process::Command;

use extendr_api::prelude::*;

/// Run beagle with a specified file
/// @export
#[extendr]
fn run_beagle(gt: &str, out: &str, n_threads: u32) {
    let beagle_path = std::env::var("BEAGLE_PATH").expect("No variable BEAGLE_PATH defined");
    let _ = Command::new("java")
        .arg("-jar")
        .arg(beagle_path)
        .arg(format!("gt={gt}"))
        .arg(format!("out={out}"))
        .arg(format!("nthreads={n_threads}"))
        .status()
        .expect("Failed to run beagle");
    if !out.ends_with(".vcf.gz") {
        std::fs::rename(format!("{out}.vcf.gz"), out).expect("Unable to rename file");
    }
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod beagler;
    fn run_beagle;
}
