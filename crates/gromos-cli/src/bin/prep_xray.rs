//! prep_xray - Prepare X-ray refinement files
//!
//! Usage: prep_xray @topo <file> @resolution <d>
//!
//! Prepares X-ray restraint topology

use gromos::io::topology::{build_topology, read_topology_file};
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: prep_xray @topo <file> [@resolution <d>] [@weight <w>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut resolution = 2.0f64; // Å
    let mut weight = 1.0f64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@resolution" => {
                i += 1;
                resolution = args[i].parse().unwrap();
            },
            "@weight" => {
                i += 1;
                weight = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);

    let n_atoms = topo.mass.len();

    eprintln!("# Preparing X-ray refinement");
    eprintln!("# Atoms: {}", n_atoms);
    eprintln!("# Resolution: {} Å", resolution);
    eprintln!("# Weight: {}", weight);

    println!("TITLE");
    println!("  X-ray refinement parameters");
    println!("END");
    println!("XRAYRES");
    println!("  # Resolution (Å)");
    println!("    {:.4}", resolution);
    println!("  # Weight");
    println!("    {:.4}", weight);
    println!("  # Number of atoms");
    println!("    {}", n_atoms);
    println!("  # B-factors (initial)");

    for i in 0..n_atoms {
        println!("    {:6} {:8.2}", i + 1, 20.0); // Default B-factor
    }

    println!("END");

    eprintln!("# X-ray preparation complete");
}
