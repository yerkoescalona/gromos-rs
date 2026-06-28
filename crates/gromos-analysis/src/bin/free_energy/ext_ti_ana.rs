//! ext_ti_ana — Thermodynamic Integration analysis (gromos-rs-compatible).
//!
//! Reads one `.trg` file per λ window, computes ⟨dH/dλ⟩ ± ee() for each,
//! integrates with the trapezoidal rule, and reports ΔG.
//!
//! Usage (GROMOS style):
//!   ext_ti_ana @trg run_lam0.trg run_lam1.trg ... @lambda 0.0 0.25 0.5 0.75 1.0
//!
//!   @trg <file> [<file> ...]   One .trg file per λ window (in λ order).
//!   @lambda <val> [<val> ...]  Explicit λ values; must match number of @trg files.
//!                              If omitted, λ is read from the first frame of each file.
//!   @skip <n>                  Discard first n frames of each window (equilibration).
//!   @verb                      Verbose: print per-frame dH/dλ.

use gromos_core::stat::Stat;
use gromos_io::read_free_energy_trajectory;
use std::process;

fn print_usage() {
    eprintln!("ext_ti_ana — Thermodynamic Integration free-energy analysis");
    eprintln!();
    eprintln!("Usage: ext_ti_ana @trg <file> [<file>...] [@lambda <val>...] [@skip <n>]");
    eprintln!();
    eprintln!("  @trg <file> [...]   .trg files, one per λ window (in ascending λ order)");
    eprintln!("  @lambda <val> [...] Explicit λ values (default: read from file headers)");
    eprintln!("  @skip <n>           Skip first n frames per window (default: 0)");
    eprintln!("  @verb               Print per-window frame listing");
}

fn main() {
    let args: Vec<String> = gromos_io::gromos_args();

    if args.len() < 3 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 3 { 1 } else { 0 });
    }

    let mut trg_files: Vec<String> = Vec::new();
    let mut explicit_lambdas: Vec<f64> = Vec::new();
    let mut skip = 0usize;
    let mut verbose = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--trg" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    trg_files.push(args[i].clone());
                    i += 1;
                }
                continue;
            },
            "--lambda" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    match args[i].parse::<f64>() {
                        Ok(v) => explicit_lambdas.push(v),
                        Err(_) => break,
                    }
                    i += 1;
                }
                continue;
            },
            "--skip" => {
                i += 1;
                skip = args[i].parse().unwrap_or(0);
            },
            "--verb" => verbose = true,
            _ => {},
        }
        i += 1;
    }

    if trg_files.is_empty() {
        eprintln!("Error: no @trg files specified");
        print_usage();
        process::exit(1);
    }

    // Read each .trg file
    struct Window {
        lambda: f64,
        stat: Stat,
    }

    let mut windows: Vec<Window> = Vec::new();

    for (wi, path) in trg_files.iter().enumerate() {
        let frames = match read_free_energy_trajectory(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Error reading {path}: {e}");
                process::exit(1);
            },
        };

        if frames.is_empty() {
            eprintln!("Warning: {path} contains no FREEENERGY03 frames — skipping");
            continue;
        }

        let lambda = if let Some(&lam) = explicit_lambdas.get(wi) {
            lam
        } else {
            frames[0].lambda
        };

        let mut stat = Stat::new();
        let n_skip = skip.min(frames.len());
        for frame in &frames[n_skip..] {
            stat.add(frame.dhdl_total);
        }

        if verbose {
            println!(
                "# λ={lambda:.6}  file={path}  frames={} (skipped={n_skip})",
                frames.len()
            );
            for (fi, frame) in frames[n_skip..].iter().enumerate() {
                println!(
                    "  frame {:5}  t={:.4}ps  dH/dλ={:.6e}",
                    fi, frame.time, frame.dhdl_total
                );
            }
        }

        windows.push(Window { lambda, stat });
    }

    if windows.is_empty() {
        eprintln!("Error: no valid λ windows loaded");
        process::exit(1);
    }

    // Sort by λ
    windows.sort_by(|a, b| a.lambda.partial_cmp(&b.lambda).unwrap());

    // Print per-window summary
    println!("#");
    println!("# ext_ti_ana — Thermodynamic Integration");
    println!(
        "# {:>10}  {:>14}  {:>14}  {:>10}",
        "lambda", "<dH/dλ>", "ee(dH/dλ)", "N"
    );
    println!("#");

    let mut lambdas: Vec<f64> = Vec::new();
    let mut aves: Vec<f64> = Vec::new();

    for w in &windows {
        let n = w.stat.n();
        if n == 0 {
            eprintln!("Warning: window λ={:.6} has 0 samples", w.lambda);
            continue;
        }
        let ave = w.stat.ave();
        let ee = if n >= 4 { w.stat.ee() } else { w.stat.rmsd() };
        println!("  {:10.6}  {:14.6e}  {:14.6e}  {:10}", w.lambda, ave, ee, n);
        lambdas.push(w.lambda);
        aves.push(ave);
    }

    if lambdas.len() < 2 {
        println!("#");
        println!("# Need ≥2 λ windows for integration");
        process::exit(0);
    }

    // Trapezoidal integration: ΔG = ∫⟨dH/dλ⟩dλ
    let mut delta_g = 0.0f64;
    for k in 0..lambdas.len() - 1 {
        let dl = lambdas[k + 1] - lambdas[k];
        delta_g += 0.5 * (aves[k] + aves[k + 1]) * dl;
    }

    println!("#");
    println!("# TI integration (trapezoidal):  ΔG = {delta_g:.6} kJ/mol");
    println!(
        "#   λ range: {:.4} – {:.4}  ({} windows)",
        lambdas[0],
        lambdas[lambdas.len() - 1],
        lambdas.len()
    );
}
