//! atominfo - List atom characteristics
//!
//! Usage: atominfo @topo <file> [@range <start-end>]
//!
//! Lists properties of atoms from topology

use gromos_rs::io::topology::{build_topology, read_topology_file};
use std::env;
use std::process;

fn print_usage() {
    eprintln!("atominfo - List atom characteristics");
    eprintln!();
    eprintln!("Usage: atominfo @topo <file> [@range <start-end>]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo     Molecular topology file");
    eprintln!("  @range    Atom range (e.g., 1-100), default: all");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  atominfo @topo system.top");
    eprintln!("  atominfo @topo system.top @range 1-100");
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 3 { 1 } else { 0 });
    }

    let mut topo_file = None;
    let mut range_start = 0;
    let mut range_end = usize::MAX;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @topo");
                    process::exit(1);
                }
                topo_file = Some(args[i].clone());
            },
            "@range" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: Missing value for @range");
                    process::exit(1);
                }
                let parts: Vec<&str> = args[i].split('-').collect();
                if parts.len() == 2 {
                    range_start = parts[0].parse().unwrap_or(1) - 1; // 0-indexed
                    range_end = parts[1].parse().unwrap_or(usize::MAX);
                }
            },
            _ => {
                eprintln!("Unknown argument: {}", args[i]);
                print_usage();
                process::exit(1);
            },
        }
        i += 1;
    }

    if topo_file.is_none() {
        eprintln!("Error: Missing @topo argument");
        print_usage();
        process::exit(1);
    }

    // Read topology
    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap_or_else(|e| {
        eprintln!("Error reading topology: {}", e);
        process::exit(1);
    });

    let topo = build_topology(topo_data);

    let num_atoms = topo.mass.len();
    range_end = range_end.min(num_atoms);

    eprintln!("# Total atoms: {}", num_atoms);
    eprintln!("# Showing atoms: {}-{}", range_start + 1, range_end);
    eprintln!();

    println!(
        "{:>8} {:>10} {:>12} {:>8}",
        "Atom#", "Mass", "Charge", "IAC"
    );
    println!("{}", "-".repeat(50));

    for i in range_start..range_end {
        println!(
            "{:8} {:10.4} {:12.6} {:8}",
            i + 1,
            topo.mass[i],
            topo.charge[i],
            topo.iac[i]
        );
    }

    eprintln!();
    eprintln!("# Listed {} atoms", range_end - range_start);
}
