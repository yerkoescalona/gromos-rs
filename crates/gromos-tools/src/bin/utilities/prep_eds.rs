//! prep_eds - Prepare EDS simulation files
//!
//! Usage: prep_eds @topo <file> @states <n>
//!
//! Prepares topology and parameters for EDS simulations

use gromos::io::topology::{build_topology, read_topology_file};
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: prep_eds @topo <file> @states <n> [@emin <E>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut n_states = 2;
    let mut emin = -100.0f64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@states" => {
                i += 1;
                n_states = args[i].parse().unwrap();
            },
            "@emin" => {
                i += 1;
                emin = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);

    eprintln!("# Preparing EDS simulation");
    eprintln!("# Number of states: {}", n_states);
    eprintln!("# E_min: {} kJ/mol", emin);

    println!("TITLE");
    println!("  EDS simulation parameters");
    println!("END");
    println!("EDS");
    println!("  # Number of states");
    println!("    {}", n_states);
    println!("  # E_min (kJ/mol)");
    println!("    {:.4}", emin);
    println!("  # State energies (initial)");

    for i in 0..n_states {
        println!("    # State {}", i + 1);
        println!("      0.0");
    }

    println!("END");

    eprintln!("# EDS preparation complete");
}
