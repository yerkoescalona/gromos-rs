//! top2build - Convert topology to build file
//!
//! Usage: top2build @topo <file>
//!
//! Converts GROMOS topology to building block format

use gromos_rs::io::topology::read_topology_file;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: top2build @topo <file>");
        process::exit(1);
    }

    let mut topo_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();

    eprintln!("# Converting topology to build file");
    eprintln!("# Atoms: {}", topo_data.masses.len());

    println!("TITLE");
    println!("  Building blocks from topology");
    println!("END");
    println!("BUILDBLOCK");
    println!("  # Number of blocks: 1");
    println!("  # Block 1");
    println!("  MOL");
    println!("  # Number of atoms");
    println!("    {}", topo_data.masses.len());

    for (i, mass) in topo_data.masses.iter().enumerate() {
        let charge = topo_data.charges.get(i).copied().unwrap_or(0.0);
        let iac = topo_data.iac.get(i).copied().unwrap_or(1);

        println!("  # Atom {}", i + 1);
        println!("    {:6} {:10.4} {:12.6} {:4}", i + 1, mass, charge, iac);
    }

    println!("END");

    eprintln!("# Conversion complete");
}
