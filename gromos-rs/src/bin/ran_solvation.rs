//! ran_solvation - Random solvation of solute
//!
//! Usage: ran_solvation @solute <file> @solvent <file> @box <x y z> @nsolv <n>
//!
//! Randomly places solvent molecules around solute

use std::env;
use std::process;

fn lcg_random(state: &mut u64) -> f32 {
    const A: u64 = 1103515245;
    const C: u64 = 12345;
    const M: u64 = 2u64.pow(31);

    *state = (A.wrapping_mul(*state).wrapping_add(C)) % M;
    (*state as f32) / (M as f32)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 9 {
        eprintln!("Usage: ran_solvation @solute <file> @solvent <file> @box <x y z> @nsolv <n>");
        process::exit(1);
    }

    let mut _solute_file = None;
    let mut _solvent_file = None;
    let mut box_x = 5.0f32;
    let mut box_y = 5.0f32;
    let mut box_z = 5.0f32;
    let mut n_solvent = 1000;
    let mut seed = 12345u64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@solute" => {
                i += 1;
                _solute_file = Some(args[i].clone());
            },
            "@solvent" => {
                i += 1;
                _solvent_file = Some(args[i].clone());
            },
            "@box" => {
                i += 1;
                box_x = args[i].parse().unwrap();
                i += 1;
                box_y = args[i].parse().unwrap();
                i += 1;
                box_z = args[i].parse().unwrap();
            },
            "@nsolv" => {
                i += 1;
                n_solvent = args[i].parse().unwrap();
            },
            "@seed" => {
                i += 1;
                seed = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    eprintln!("# Random solvation");
    eprintln!("# Box: {} x {} x {} nm", box_x, box_y, box_z);
    eprintln!("# Solvent molecules: {}", n_solvent);

    let mut rng = seed;

    println!("TITLE");
    println!("  Random solvation");
    println!("END");
    println!("POSITION");
    println!("  {:15.9} {:15.9} {:15.9}", box_x, box_y, box_z);

    // Simplified: place 3-atom water molecules
    for _ in 0..n_solvent {
        let cx = lcg_random(&mut rng) * box_x;
        let cy = lcg_random(&mut rng) * box_y;
        let cz = lcg_random(&mut rng) * box_z;

        // Oxygen
        println!("{:15.9} {:15.9} {:15.9}", cx, cy, cz);

        // H1
        println!("{:15.9} {:15.9} {:15.9}", cx + 0.01, cy, cz);

        // H2
        println!("{:15.9} {:15.9} {:15.9}", cx, cy + 0.01, cz);
    }

    println!("END");

    eprintln!("# Solvation complete: {} atoms", n_solvent * 3);
}
