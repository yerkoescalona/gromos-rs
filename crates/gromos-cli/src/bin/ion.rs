//! ion - Ion analysis and distribution
//!
//! Usage: ion @topo <file> @traj <file> @ion <atom>
//!
//! Analyzes ion positions and coordination

use gromos::io::topology::{build_topology, read_topology_file};
use gromos::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: ion @topo <file> @traj <file> @ion <atom_idx>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut ion_idx = 0;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@ion" => {
                i += 1;
                ion_idx = args[i].parse::<usize>().unwrap() - 1; // Convert to 0-indexed
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    let num_atoms = topo.mass.len();

    if ion_idx >= num_atoms {
        eprintln!(
            "Error: Ion index {} out of range (max {})",
            ion_idx + 1,
            num_atoms
        );
        process::exit(1);
    }

    eprintln!("# Ion analysis");
    eprintln!("# Ion atom: {}", ion_idx + 1);
    eprintln!("# Charge: {}", topo.charge[ion_idx]);

    println!("# Time (ps)    Ion_X (nm)    Ion_Y (nm)    Ion_Z (nm)    Coordination");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let ion_pos = frame.positions[ion_idx];

                // Count coordination (atoms within 0.35 nm)
                let cutoff = 0.35f32;
                let cutoff_sq = cutoff * cutoff;
                let mut coordination = 0;

                for j in 0..num_atoms {
                    if j != ion_idx {
                        let dx = frame.positions[j].x - ion_pos.x;
                        let dy = frame.positions[j].y - ion_pos.y;
                        let dz = frame.positions[j].z - ion_pos.z;
                        let dist_sq = dx * dx + dy * dy + dz * dz;

                        if dist_sq <= cutoff_sq {
                            coordination += 1;
                        }
                    }
                }

                println!(
                    "{:12.4} {:13.6} {:13.6} {:13.6} {:14}",
                    frame.time, ion_pos.x, ion_pos.y, ion_pos.z, coordination
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
