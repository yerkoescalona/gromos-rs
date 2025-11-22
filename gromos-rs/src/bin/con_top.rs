//! con_top - Concatenate topology files
//!
//! Usage: con_top @topo <file1> @topo <file2> ...
//!
//! Combines multiple topology files into one

use gromos_rs::io::topology::{build_topology, read_topology_file};
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: con_top @topo <file1> @topo <file2> ...");
        process::exit(1);
    }

    let mut topo_files = Vec::new();

    let mut i = 1;
    while i < args.len() {
        if args[i] == "@topo" {
            i += 1;
            if i < args.len() {
                topo_files.push(args[i].clone());
            }
        }
        i += 1;
    }

    if topo_files.len() < 2 {
        eprintln!("Error: Need at least 2 topology files");
        process::exit(1);
    }

    let mut all_topos = Vec::new();

    for topo_file in &topo_files {
        let topo_data = read_topology_file(topo_file).unwrap();
        let topo = build_topology(topo_data);
        all_topos.push(topo);
        eprintln!(
            "Read: {} ({} atoms)",
            topo_file,
            all_topos.last().unwrap().mass.len()
        );
    }

    let total_atoms: usize = all_topos.iter().map(|t| t.mass.len()).sum();
    let total_bonds: usize = all_topos.iter().map(|t| t.solute.bonds.len()).sum();

    println!("TITLE");
    println!("Combined topology from {} files", topo_files.len());
    println!("END");
    println!();
    println!("# Total atoms: {}", total_atoms);
    println!("# Total bonds: {}", total_bonds);
    println!();

    let mut atom_offset = 0;
    for (idx, topo) in all_topos.iter().enumerate() {
        println!("# From: {}", topo_files[idx]);
        for (i, &mass) in topo.mass.iter().enumerate() {
            println!(
                "# Atom {:6}: Mass {:8.4}, Charge {:10.6}",
                atom_offset + i + 1,
                mass,
                topo.charge[i]
            );
        }
        atom_offset += topo.mass.len();
    }

    eprintln!();
    eprintln!(
        "Combined {} atoms from {} topologies",
        total_atoms,
        topo_files.len()
    );
}
