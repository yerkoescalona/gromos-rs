//! ext_ti_merge — Merge TI data from overlapping lambda windows.
//!
//! Reads `ext_ti_ana` output files (one per sampled λ), each containing
//! ⟨dH/dλ⟩ values at one or more prediction λ values.  Combines overlapping
//! windows via linear interpolation to produce a unified ⟨dH/dλ⟩ vs λ curve,
//! then integrates with the trapezoidal rule to give ΔG.
//!
//! The interpolation formula (GROMOS ext_ti_merge.cc):
//!
//!   w₁ = (λ_P − λ₂) / (λ₁ − λ₂)
//!   w₂ = (λ_P − λ₁) / (λ₂ − λ₁)
//!   ⟨dH/dλ⟩_P = w₁·⟨dH/dλ⟩₁ + w₂·⟨dH/dλ⟩₂
//!
//! Usage:
//!   ext_ti_merge @files <out1> <out2> ... @slam <λ₁> <λ₂> ...
//!
//!   @files   ext_ti_ana output files, one per sampled λ window
//!   @slam    Sampled λ values matching the file order
//!   @plam    Prediction λ values (default: same as @slam)

use gromos_io::gromos_args;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn print_usage() {
    eprintln!("ext_ti_merge — merge TI data from overlapping lambda windows");
    eprintln!();
    eprintln!("Usage: ext_ti_merge @files <f1> [f2...] @slam <λ1> [λ2...] [@plam <λ1> [λ2...]]");
    eprintln!();
    eprintln!("  @files  ext_ti_ana output files (one per sampled λ)");
    eprintln!("  @slam   Sampled λ values (same order as @files)");
    eprintln!("  @plam   Prediction λ values (default: same as @slam)");
}

/// Parse ⟨dH/dλ⟩ from an ext_ti_ana output file.
/// Expected format: lines with  λ  ⟨dH/dλ⟩  [ee]  [N]  (skip comment lines).
fn parse_dhdl_file(path: &str) -> Vec<(f64, f64)> {
    let file = File::open(path).unwrap_or_else(|e| {
        eprintln!("Error opening '{path}': {e}");
        process::exit(1);
    });
    let mut result = Vec::new();
    for line in BufReader::new(file).lines().filter_map(|l| l.ok()) {
        let t = line.trim();
        if t.starts_with('#') || t.is_empty() {
            continue;
        }
        let cols: Vec<f64> = t
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        if cols.len() >= 2 {
            result.push((cols[0], cols[1])); // (λ, ⟨dH/dλ⟩)
        }
    }
    result
}

/// Linear interpolation weight: value at x interpolated between (x1,y1) and (x2,y2).
fn lerp(x: f64, x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let w1 = (x - x2) / (x1 - x2);
    let w2 = (x - x1) / (x2 - x1);
    w1 * y1 + w2 * y2
}

fn main() {
    let args = gromos_args();

    if args.len() < 2 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut files: Vec<String> = Vec::new();
    let mut slam: Vec<f64> = Vec::new();
    let mut plam: Vec<f64> = Vec::new();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--files" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    files.push(args[i].clone());
                    i += 1;
                }
                continue;
            },
            "--slam" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    if let Ok(v) = args[i].parse::<f64>() {
                        slam.push(v);
                    } else {
                        break;
                    }
                    i += 1;
                }
                continue;
            },
            "--plam" => {
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    if let Ok(v) = args[i].parse::<f64>() {
                        plam.push(v);
                    } else {
                        break;
                    }
                    i += 1;
                }
                continue;
            },
            other if other.starts_with("--") => {
                eprintln!("Unknown argument: {other}");
                process::exit(1);
            },
            _ => {},
        }
        i += 1;
    }

    if files.is_empty() {
        eprintln!("Error: @files required");
        process::exit(1);
    }
    if slam.is_empty() {
        eprintln!("Error: @slam required");
        process::exit(1);
    }
    if files.len() != slam.len() {
        eprintln!(
            "Error: @files ({}) and @slam ({}) must have the same count",
            files.len(),
            slam.len()
        );
        process::exit(1);
    }
    if plam.is_empty() {
        plam = slam.clone();
    }

    // Read all files: data[k] = list of (λ, ⟨dH/dλ⟩) for sampled λ slam[k]
    let data: Vec<Vec<(f64, f64)>> = files.iter().map(|f| parse_dhdl_file(f)).collect();

    // For each sampled λ, build a lookup: slam[k] -> dhdl (the value at that λ, from that file)
    let mut dhdl_at_slam: Vec<f64> = Vec::new();
    for (k, entries) in data.iter().enumerate() {
        // Find the entry in this file closest to slam[k]
        let target = slam[k];
        let best = entries.iter().min_by(|a, b| {
            (a.0 - target)
                .abs()
                .partial_cmp(&(b.0 - target).abs())
                .unwrap()
        });
        let dhdl = best.map(|e| e.1).unwrap_or_else(|| {
            eprintln!("Warning: no data near λ={target} in {}", files[k]);
            0.0
        });
        dhdl_at_slam.push(dhdl);
    }

    // For each prediction λ, interpolate from the two nearest sampled λ windows
    let mut pred: Vec<(f64, f64)> = Vec::new();
    for &lp in &plam {
        // Find surrounding sampled λ values
        let dhdl_p = if slam.len() == 1 {
            dhdl_at_slam[0]
        } else {
            // Find bracketing pair (λ₁, λ₂) with λ₁ ≤ lp ≤ λ₂
            let pos = slam.partition_point(|&s| s <= lp);
            let (k1, k2) = if pos == 0 {
                (0, 1)
            } else if pos >= slam.len() {
                (slam.len() - 2, slam.len() - 1)
            } else {
                (pos - 1, pos)
            };
            lerp(lp, slam[k1], dhdl_at_slam[k1], slam[k2], dhdl_at_slam[k2])
        };
        pred.push((lp, dhdl_p));
    }

    // Sort by λ and integrate
    pred.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    println!("# ext_ti_merge — merged TI ⟨dH/dλ⟩");
    println!("# sampled λ: {:?}", slam);
    println!("# {:>10}  {:>14}", "lambda", "<dH/dλ>(kJ/mol)");
    for (lp, dhdl) in &pred {
        println!("  {:10.6}  {:14.6e}", lp, dhdl);
    }

    // Trapezoidal ΔG
    if pred.len() >= 2 {
        let dg: f64 = pred
            .windows(2)
            .map(|w| 0.5 * (w[0].1 + w[1].1) * (w[1].0 - w[0].0))
            .sum();
        println!("#");
        println!(
            "# ΔG (trapezoidal) = {dg:.6} kJ/mol  ({} points)",
            pred.len()
        );
    }
}
