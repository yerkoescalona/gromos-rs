//! ext_ti_ana - Extended Thermodynamic Integration analysis
//!
//! Usage: ext_ti_ana @energy <file> @lambda <values>
//!
//! Analyzes TI simulations and calculates free energy

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: ext_ti_ana @energy <file> [@lambda <start> <end> <n>]");
        process::exit(1);
    }

    let mut energy_file = None;
    let mut lambda_start = 0.0f64;
    let mut lambda_end = 1.0f64;
    let mut n_lambda = 11;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@energy" => {
                i += 1;
                energy_file = Some(args[i].clone());
            },
            "@lambda" => {
                i += 1;
                lambda_start = args[i].parse().unwrap();
                i += 1;
                lambda_end = args[i].parse().unwrap();
                i += 1;
                n_lambda = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    // Generate lambda values
    let mut lambdas = Vec::new();
    for i in 0..n_lambda {
        let lambda =
            lambda_start + (lambda_end - lambda_start) * (i as f64) / ((n_lambda - 1) as f64);
        lambdas.push(lambda);
    }

    // Read dH/dlambda values from file
    let mut dhdl_values = Vec::new();
    if let Some(path) = energy_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(line_str) = line {
                    let parts: Vec<&str> = line_str.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let Ok(val) = parts[1].parse::<f64>() {
                            dhdl_values.push(val);
                        }
                    }
                }
            }
        }
    }

    if dhdl_values.is_empty() {
        eprintln!("# Warning: No dH/dlambda data loaded, using synthetic data");
        for _ in 0..lambdas.len() {
            dhdl_values.push(10.0); // Placeholder
        }
    }

    eprintln!("# Extended TI Analysis");
    eprintln!("# Lambda points: {}", lambdas.len());
    eprintln!("# dH/dlambda samples: {}", dhdl_values.len());

    println!("# Lambda    <dH/dlambda> (kJ/mol)");

    // Average dH/dlambda for each lambda window
    let samples_per_lambda = dhdl_values.len() / lambdas.len().max(1);
    let mut avg_dhdl = Vec::new();

    for (idx, &lambda) in lambdas.iter().enumerate() {
        let start_idx = idx * samples_per_lambda;
        let end_idx = ((idx + 1) * samples_per_lambda).min(dhdl_values.len());

        let avg = if end_idx > start_idx {
            let sum: f64 = dhdl_values[start_idx..end_idx].iter().sum();
            sum / (end_idx - start_idx) as f64
        } else {
            0.0
        };

        avg_dhdl.push(avg);
        println!("{:8.4} {:18.6}", lambda, avg);
    }

    // Integrate using trapezoidal rule
    let mut delta_g = 0.0f64;
    for i in 0..(lambdas.len() - 1) {
        let dlambda = lambdas[i + 1] - lambdas[i];
        delta_g += 0.5 * (avg_dhdl[i] + avg_dhdl[i + 1]) * dlambda;
    }

    println!("\n# Free Energy Change (TI): {:.6} kJ/mol", delta_g);
}
