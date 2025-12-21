//! bar - Bennett Acceptance Ratio for free energy calculations
//!
//! Usage: bar @forward <file> @backward <file> @temp <T>
//!
//! Calculates free energy difference using BAR method

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

const R: f64 = 0.008314462618; // Gas constant in kJ/(mol·K)

fn fermi(x: f64) -> f64 {
    1.0 / (1.0 + x.exp())
}

fn bar_iteration(dw_fwd: &[f64], dw_bwd: &[f64], beta: f64, c: f64) -> f64 {
    let n_fwd = dw_fwd.len() as f64;
    let n_bwd = dw_bwd.len() as f64;

    let sum_fwd: f64 = dw_fwd.iter().map(|&dw| fermi(beta * (dw - c))).sum();

    let sum_bwd: f64 = dw_bwd.iter().map(|&dw| fermi(-beta * (dw + c))).sum();

    c + (1.0 / beta) * ((n_fwd / n_bwd) * (sum_bwd / sum_fwd)).ln()
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: bar @forward <file> @backward <file> @temp <T>");
        process::exit(1);
    }

    let mut fwd_file = None;
    let mut bwd_file = None;
    let mut temperature = 300.0;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@forward" => {
                i += 1;
                fwd_file = Some(args[i].clone());
            },
            "@backward" => {
                i += 1;
                bwd_file = Some(args[i].clone());
            },
            "@temp" => {
                i += 1;
                temperature = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    // Read energy differences from files
    let mut dw_fwd = Vec::new();
    if let Some(path) = fwd_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(line_str) = line {
                    if let Ok(val) = line_str.trim().parse::<f64>() {
                        dw_fwd.push(val);
                    }
                }
            }
        }
    }

    let mut dw_bwd = Vec::new();
    if let Some(path) = bwd_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(line_str) = line {
                    if let Ok(val) = line_str.trim().parse::<f64>() {
                        dw_bwd.push(val);
                    }
                }
            }
        }
    }

    if dw_fwd.is_empty() || dw_bwd.is_empty() {
        eprintln!("Error: No data loaded");
        process::exit(1);
    }

    let beta = 1.0 / (R * temperature);

    eprintln!("# BAR Free Energy Calculation");
    eprintln!("# Temperature: {} K", temperature);
    eprintln!("# Forward samples: {}", dw_fwd.len());
    eprintln!("# Backward samples: {}", dw_bwd.len());

    // Initial guess: average of forward and backward estimates
    let avg_fwd: f64 = dw_fwd.iter().sum::<f64>() / dw_fwd.len() as f64;
    let avg_bwd: f64 = dw_bwd.iter().sum::<f64>() / dw_bwd.len() as f64;
    let mut dg = (avg_fwd - avg_bwd) / 2.0;

    println!("# Iteration    dG (kJ/mol)    Change");

    // Iterate BAR equation
    for iter in 0..100 {
        let dg_new = bar_iteration(&dw_fwd, &dw_bwd, beta, dg);
        let change = (dg_new - dg).abs();

        println!("{:10} {:15.6} {:15.6}", iter + 1, dg_new, change);

        if change < 1e-6 {
            eprintln!("# Converged after {} iterations", iter + 1);
            break;
        }

        dg = dg_new;
    }

    println!("\n# Final Free Energy Difference: {:.6} kJ/mol", dg);
}
