//! make_top - Build GROMOS topology from scratch
//!
//! Usage: make_top @natoms <n> @mass <m> @charge <q>
//!
//! Creates a simple topology file

use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: make_top @natoms <n> [@mass <m>] [@charge <q>]");
        process::exit(1);
    }

    let mut natoms = 100;
    let mut mass = 12.0f64;
    let mut charge = 0.0f64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@natoms" => {
                i += 1;
                natoms = args[i].parse().unwrap();
            },
            "@mass" => {
                i += 1;
                mass = args[i].parse().unwrap();
            },
            "@charge" => {
                i += 1;
                charge = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    eprintln!("# Creating topology with {} atoms", natoms);
    eprintln!("# Mass: {}, Charge: {}", mass, charge);

    println!("TITLE");
    println!("  Generated topology");
    println!("END");
    println!("PHYSICALCONSTANTS");
    println!("  # FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)");
    println!("    138.9354");
    println!("  # HBAR: Planck's constant HBAR = H/(2* PI)");
    println!("    0.0635078");
    println!("  # SPDL: Speed of light (nm/ps)");
    println!("    299.7925");
    println!("  # BOLTZ: Boltzmann's constant");
    println!("    0.00831441");
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
    println!("    {}", natoms);

    for i in 1..=natoms {
        println!("  # {:6} RES  ATOM  IAC   MASS    CHARGE", i);
        println!("    {:6}   1     1    1  {:8.4}  {:10.6}", i, mass, charge);
    }

    println!("END");
    println!("SOLUTEMOLECULES");
    println!("  # NSPM: number of separate molecules");
    println!("    1");
    println!("  # NSP[1...NSPM]: atom sequence number of last atom");
    println!("    {}", natoms);
    println!("END");
}
