//! ener - Recalculate interaction energies from trajectory
//!
//! Usage: ener @topo <file> @traj <file>
//!
//! Calculates nonbonded interaction energies

use gromos_rs::io::topology::{build_topology, read_topology_file};
use gromos_rs::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

const COULOMB_CONST: f32 = 138.9354859; // kJ/(mol·nm·e²)

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 5 {
        eprintln!("Usage: ener @topo <file> @traj <file>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;

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
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    let num_atoms = topo.mass.len();

    eprintln!("# Energy recalculation");
    eprintln!("# Atoms: {}", num_atoms);

    println!("# Time (ps)    E_elec (kJ/mol)    E_vdw (kJ/mol)    E_tot (kJ/mol)");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let mut e_elec = 0.0f32;
                let mut e_vdw = 0.0f32;

                for i in 0..num_atoms {
                    for j in (i + 1)..num_atoms {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist_sq = dx * dx + dy * dy + dz * dz;
                        let dist = dist_sq.sqrt();

                        if dist > 0.01 {
                            // Electrostatic energy
                            let qi = topo.charge[i] as f32;
                            let qj = topo.charge[j] as f32;
                            e_elec += COULOMB_CONST * qi * qj / dist;

                            // Simple LJ-like VDW (placeholder)
                            let sigma = 0.3f32; // nm
                            let epsilon = 1.0f32; // kJ/mol
                            let sr = sigma / dist;
                            let sr6 = sr.powi(6);
                            let sr12 = sr6 * sr6;
                            e_vdw += 4.0 * epsilon * (sr12 - sr6);
                        }
                    }
                }

                let e_tot = e_elec + e_vdw;

                println!(
                    "{:12.4} {:18.4} {:17.4} {:17.4}",
                    frame.time, e_elec, e_vdw, e_tot
                );
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
