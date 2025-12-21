//! addvirt_top - Add virtual atoms to topology
//!
//! Usage: addvirt_top @topo <file> @nvirt <n>
//!
//! Adds virtual (massless) atoms to topology

use gromos::io::topology::read_topology_file;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: addvirt_top @topo <file> @nvirt <n>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut n_virtual = 10;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@nvirt" => {
                i += 1;
                n_virtual = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();

    let original_atoms = topo_data.masses.len();
    let total_atoms = original_atoms + n_virtual;

    eprintln!("# Adding {} virtual atoms to topology", n_virtual);
    eprintln!("# Original atoms: {}", original_atoms);
    eprintln!("# Total atoms: {}", total_atoms);

    println!("TITLE");
    println!("  Topology with virtual atoms");
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

    // Print original atoms
    for i in 0..original_atoms {
        println!("  # {:6} RES  ATOM  IAC   MASS    CHARGE", i + 1);
        println!(
            "    {:6}   1     1  {:3}  {:8.4}  {:10.6}",
            i + 1,
            topo_data.iac[i],
            topo_data.masses[i],
            topo_data.charges[i]
        );
    }

    // Print virtual atoms (mass = 0)
    for i in 0..n_virtual {
        let atom_idx = original_atoms + i + 1;
        println!(
            "  # {:6} RES  ATOM  IAC   MASS    CHARGE (virtual)",
            atom_idx
        );
        println!("    {:6}   1     1    1    0.0000      0.000000", atom_idx);
    }

    println!("END");
    println!("SOLUTEMOLECULES");
    println!("  # NSPM: number of separate molecules");
    println!("    1");
    println!("  # NSP[1...NSPM]: atom sequence number of last atom");
    println!("    {}", total_atoms);
    println!("END");
}
