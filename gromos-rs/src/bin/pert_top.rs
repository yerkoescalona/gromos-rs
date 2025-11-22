//! pert_top - Create perturbation topology
//!
//! Usage: pert_top @topo1 <file> @topo2 <file>
//!
//! Creates perturbation topology for free energy calculations

use gromos_rs::io::topology::read_topology_file;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: pert_top @topo1 <file> @topo2 <file>");
        process::exit(1);
    }

    let mut topo1_file = None;
    let mut topo2_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo1" => {
                i += 1;
                topo1_file = Some(args[i].clone());
            },
            "@topo2" => {
                i += 1;
                topo2_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let topo1_data = read_topology_file(&topo1_file.unwrap()).unwrap();
    let topo2_data = read_topology_file(&topo2_file.unwrap()).unwrap();

    let n1 = topo1_data.masses.len();
    let n2 = topo2_data.masses.len();
    let n_atoms = n1.max(n2);

    eprintln!("# Creating perturbation topology");
    eprintln!("# State A atoms: {}", n1);
    eprintln!("# State B atoms: {}", n2);

    println!("TITLE");
    println!("  Perturbation topology");
    println!("END");
    println!("PERTATOM");
    println!("  # Number of perturbed atoms");
    println!("    {}", n_atoms);

    for i in 0..n_atoms {
        let mass_a = if i < n1 { topo1_data.masses[i] } else { 0.0 };
        let charge_a = if i < n1 { topo1_data.charges[i] } else { 0.0 };
        let iac_a = if i < n1 { topo1_data.iac[i] } else { 0 };

        let mass_b = if i < n2 { topo2_data.masses[i] } else { 0.0 };
        let charge_b = if i < n2 { topo2_data.charges[i] } else { 0.0 };
        let iac_b = if i < n2 { topo2_data.iac[i] } else { 0 };

        println!("  # Atom {}", i + 1);
        println!(
            "    # State A: mass={:.4} charge={:.6} IAC={}",
            mass_a, charge_a, iac_a
        );
        println!(
            "    # State B: mass={:.4} charge={:.6} IAC={}",
            mass_b, charge_b, iac_b
        );
        println!(
            "    {:6} {:10.4} {:12.6} {:4} {:10.4} {:12.6} {:4}",
            i + 1,
            mass_a,
            charge_a,
            iac_a,
            mass_b,
            charge_b,
            iac_b
        );
    }

    println!("END");

    eprintln!("# Perturbation topology created");
}
