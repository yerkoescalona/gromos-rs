//! prep_noe - Prepare NOE restraint files
//!
//! Usage: prep_noe @pairs <file> @force <k>
//!
//! Prepares NOE restraint topology

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: prep_noe @pairs <file> [@force <k>] [@dist <r0>]");
        process::exit(1);
    }

    let mut pairs_file = None;
    let mut force_const = 1000.0f64; // kJ/(mol·nm²)
    let mut target_dist = 0.3f64; // nm

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@pairs" => {
                i += 1;
                pairs_file = Some(args[i].clone());
            },
            "@force" => {
                i += 1;
                force_const = args[i].parse().unwrap();
            },
            "@dist" => {
                i += 1;
                target_dist = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut pairs = Vec::new();
    if let Some(path) = pairs_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(line_str) = line {
                    let parts: Vec<&str> = line_str.split_whitespace().collect();
                    if parts.len() >= 2 {
                        if let (Ok(i), Ok(j)) =
                            (parts[0].parse::<usize>(), parts[1].parse::<usize>())
                        {
                            pairs.push((i, j));
                        }
                    }
                }
            }
        }
    }

    eprintln!("# Preparing NOE restraints");
    eprintln!("# Pairs: {}", pairs.len());
    eprintln!("# Force constant: {} kJ/(mol·nm²)", force_const);
    eprintln!("# Target distance: {} nm", target_dist);

    println!("TITLE");
    println!("  NOE distance restraints");
    println!("END");
    println!("DISTANCERES");
    println!("  # Number of restraints");
    println!("    {}", pairs.len());

    for (idx, &(i, j)) in pairs.iter().enumerate() {
        println!("  # Restraint {}: atoms {} - {}", idx + 1, i, j);
        println!(
            "    {:6} {:6} {:12.4} {:12.4}",
            i, j, target_dist, force_const
        );
    }

    println!("END");

    eprintln!("# NOE preparation complete");
}
