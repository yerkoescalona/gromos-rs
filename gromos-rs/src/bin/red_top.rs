//! red_top - Reduce topology to selected atoms
//!
//! Usage: red_top @topo <file> @natoms <n>
//!
//! Creates reduced topology with first N atoms

use gromos_rs::io::topology::{build_topology, read_topology_file};
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: red_top @topo <file> @natoms <n>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut natoms = 0;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@natoms" => {
                i += 1;
                natoms = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);

    let total_atoms = topo.mass.len();
    let keep_atoms = natoms.min(total_atoms);

    eprintln!("Original topology: {} atoms", total_atoms);
    eprintln!("Reduced topology: {} atoms", keep_atoms);

    println!("TITLE");
    println!("Reduced topology: first {} atoms", keep_atoms);
    println!("END");
    println!();
    println!("# Atoms:");

    for i in 0..keep_atoms {
        println!(
            "# {:6}: Mass {:8.4}, Charge {:10.6}, IAC {:4}",
            i + 1,
            topo.mass[i],
            topo.charge[i],
            topo.iac[i]
        );
    }

    println!();
    println!("# Total: {} atoms", keep_atoms);
}
