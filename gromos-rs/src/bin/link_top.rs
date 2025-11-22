//! link_top - Link multiple topology files
//!
//! Usage: link_top @topo1 <file> @topo2 <file>
//!
//! Combines topologies into a single file

use gromos_rs::io::topology::read_topology_file;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: link_top @topo1 <file> @topo2 <file> [@topo3 <file> ...]");
        process::exit(1);
    }

    let mut topo_files = Vec::new();

    let mut i = 1;
    while i < args.len() {
        if args[i].starts_with("@topo") {
            i += 1;
            if i < args.len() {
                topo_files.push(args[i].clone());
            }
        }
        i += 1;
    }

    if topo_files.is_empty() {
        eprintln!("Error: No topology files specified");
        process::exit(1);
    }

    eprintln!("# Linking {} topology files", topo_files.len());

    let mut all_masses = Vec::new();
    let mut all_charges = Vec::new();
    let mut all_iacs = Vec::new();

    for (idx, topo_file) in topo_files.iter().enumerate() {
        eprintln!("# Reading topology {}: {}", idx + 1, topo_file);

        match read_topology_file(topo_file) {
            Ok(topo_data) => {
                all_masses.extend(topo_data.masses);
                all_charges.extend(topo_data.charges);
                all_iacs.extend(topo_data.iac);
            },
            Err(e) => {
                eprintln!("Warning: Could not read {}: {:?}", topo_file, e);
            },
        }
    }

    let total_atoms = all_masses.len();
    eprintln!("# Total atoms in linked topology: {}", total_atoms);

    println!("TITLE");
    println!("  Linked topology");
    println!("END");
    println!("TOPVERSION");
    println!("    2.0");
    println!("END");
    println!("ATOMTYPENAME");
    println!("    C");
    println!("END");
    println!("RESNAME");
    println!("    MOL");
    println!("END");
    println!("SOLUTEATOM");
    println!("  # NRP: number of solute atoms");
    println!("    {}", total_atoms);

    for i in 0..total_atoms {
        println!("  # {:6} RES  ATOM  IAC   MASS    CHARGE", i + 1);
        println!(
            "    {:6}   1     1  {:3}  {:8.4}  {:10.6}",
            i + 1,
            all_iacs[i],
            all_masses[i],
            all_charges[i]
        );
    }

    println!("END");
    println!("SOLUTEMOLECULES");
    println!("  # NSPM: number of separate molecules");
    println!("    1");
    println!("  # NSP[1...NSPM]: atom sequence number of last atom");
    println!("    {}", total_atoms);
    println!("END");
}
