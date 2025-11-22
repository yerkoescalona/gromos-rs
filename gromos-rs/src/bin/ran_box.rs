//! ran_box - Generate random box configuration
//!
//! Usage: ran_box @topo <file> @box <x y z> @seed <n>
//!
//! Creates random positions within box dimensions

use gromos_rs::io::topology::{build_topology, read_topology_file};
use std::env;
use std::process;

fn lcg_random(state: &mut u64) -> f32 {
    // Linear congruential generator
    const A: u64 = 1103515245;
    const C: u64 = 12345;
    const M: u64 = 2u64.pow(31);

    *state = (A.wrapping_mul(*state).wrapping_add(C)) % M;
    (*state as f32) / (M as f32)
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: ran_box @topo <file> @box <x y z> [@seed <n>]");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut box_x = 3.0f32;
    let mut box_y = 3.0f32;
    let mut box_z = 3.0f32;
    let mut seed = 12345u64;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@box" => {
                i += 1;
                box_x = args[i].parse().unwrap();
                i += 1;
                box_y = args[i].parse().unwrap();
                i += 1;
                box_z = args[i].parse().unwrap();
            },
            "@seed" => {
                i += 1;
                seed = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let num_atoms = topo.mass.len();

    eprintln!("# Generating random positions for {} atoms", num_atoms);
    eprintln!("# Box: {} x {} x {} nm", box_x, box_y, box_z);
    eprintln!("# Seed: {}", seed);

    let mut rng_state = seed;

    println!("TITLE");
    println!("  Random box configuration");
    println!("END");
    println!("POSITION");
    println!("  {:15.9} {:15.9} {:15.9}", box_x, box_y, box_z);

    for _ in 0..num_atoms {
        let x = lcg_random(&mut rng_state) * box_x;
        let y = lcg_random(&mut rng_state) * box_y;
        let z = lcg_random(&mut rng_state) * box_z;
        println!("{:15.9} {:15.9} {:15.9}", x, y, z);
    }

    println!("END");
}
