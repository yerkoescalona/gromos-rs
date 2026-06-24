//! atominfo — Report information about selected atoms (gromos++ compatible).
//!
//! Usage:
//!   atominfo @topo <topology> @atomspec <spec> [@sort]
//!
//! Reads a GROMOS topology, resolves the AtomSpecifier, and prints the TITLE +
//! ATOMS block in gromos++ format for every selected atom.
//!
//! Output columns (ATOMS block):
//!   mol:atom  gromos_nr  residue_nr  residue_name  atom_name  iac(1-based)  charge

use gromos_core::{
    selection::AtomSelection,
    topology::Topology,
};
use gromos_io::{gromos_args, topology::build_topology, topology::read_topology_file};
use std::process;

fn print_usage() {
    eprintln!("atominfo — report information about selected atoms");
    eprintln!();
    eprintln!("Usage: atominfo @topo <topology> @atomspec <spec> [@sort]");
    eprintln!();
    eprintln!("  @topo     <file>    GROMOS topology (.top/.topo)");
    eprintln!("  @atomspec <spec>    gromos++ AtomSpecifier string");
    eprintln!("  @sort               sort selected atoms by global index (default: yes)");
}

fn main() {
    let args = gromos_args();

    if args.len() < 2 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut topo_file = None;
    let mut atomspec = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--topo"     => { i += 1; topo_file   = Some(args[i].clone()); }
            "--atomspec" => { i += 1; atomspec     = Some(args[i].clone()); }
            "--sort"     => {} // always sort
            other if other.starts_with("--") => {
                eprintln!("Unknown argument: {other}");
                process::exit(1);
            }
            _ => {}
        }
        i += 1;
    }

    let topo_path = topo_file.unwrap_or_else(|| { eprintln!("Error: @topo required"); process::exit(1); });
    let spec      = atomspec.unwrap_or_else(|| { eprintln!("Error: @atomspec required"); process::exit(1); });

    // Load topology (no solvation needed for atominfo)
    let topo_data = read_topology_file(&topo_path).unwrap_or_else(|e| {
        eprintln!("Error reading topology: {e}");
        process::exit(1);
    });
    let topo = build_topology(topo_data);

    // Resolve selection
    let selection = AtomSelection::from_string(&spec, &topo).unwrap_or_else(|e| {
        eprintln!("Error in atom specifier '{spec}': {e}");
        process::exit(1);
    });

    // Print TITLE block
    println!("TITLE");
    println!("\tatominfo");
    println!("\t{spec}");
    println!("END");

    // Print ATOMS block
    println!("ATOMS");
    for &global_idx in selection.indices() {
        let mol_nr   = topo.molecule_nr(global_idx).map(|m| m + 1).unwrap_or(1); // 1-based
        let local_nr = if topo.molecules.is_empty() {
            global_idx + 1
        } else {
            let start = topo.molecules.get(mol_nr - 1).map(|r| r.start).unwrap_or(0);
            global_idx - start + 1
        };
        let gromos_nr   = global_idx + 1;
        let res_nr      = topo.residue_nr(global_idx).unwrap_or(0);
        let res_name    = topo.residue_name(global_idx).unwrap_or("???");
        let atom_name   = topo.atom_name(global_idx).unwrap_or("???");
        let iac         = topo.iac.get(global_idx).map(|&v| v + 1).unwrap_or(0); // 1-based
        let charge      = topo.charge.get(global_idx).copied().unwrap_or(0.0);

        println!(
            "{:>11}  {:>9}  {:>9}  {:>10}  {:>9}  {:>9}  {:>9}",
            format!("{mol_nr}:{local_nr}"),
            gromos_nr,
            res_nr,
            res_name,
            atom_name,
            iac,
            charge,
        );
    }
    println!("END");
}
