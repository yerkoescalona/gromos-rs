//! copy_box - Duplicate simulation box in space
//!
//! Usage: copy_box @topo <file> @traj <file> @replicate <nx ny nz>
//!
//! Creates replicated box by copying molecules in x, y, z directions

use gromos::io::topology::{build_topology, read_topology_file};
use gromos::io::trajectory::TrajectoryReader;
use gromos::math::Vec3;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: copy_box @topo <file> @traj <file> @replicate <nx ny nz>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut nx = 1;
    let mut ny = 1;
    let mut nz = 1;

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
            "@replicate" => {
                i += 1;
                nx = args[i].parse().unwrap();
                i += 1;
                ny = args[i].parse().unwrap();
                i += 1;
                nz = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    let num_atoms = topo.mass.len();
    eprintln!("# Replicating box {} x {} x {} times", nx, ny, nz);
    eprintln!("# Original atoms: {}", num_atoms);
    eprintln!(
        "# Total atoms after replication: {}",
        num_atoms * nx * ny * nz
    );

    println!("TITLE");
    println!("  Replicated simulation box");
    println!("END");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let box_dims = frame.box_dims;
                let new_box = Vec3::new(
                    box_dims.x * nx as f32,
                    box_dims.y * ny as f32,
                    box_dims.z * nz as f32,
                );

                println!("POSITION");
                println!("  {:15.9} {:15.9} {:15.9}", new_box.x, new_box.y, new_box.z);

                for ix in 0..nx {
                    for iy in 0..ny {
                        for iz in 0..nz {
                            let offset = Vec3::new(
                                ix as f32 * box_dims.x,
                                iy as f32 * box_dims.y,
                                iz as f32 * box_dims.z,
                            );

                            for pos in &frame.positions {
                                let new_pos = *pos + offset;
                                println!(
                                    "{:15.9} {:15.9} {:15.9}",
                                    new_pos.x, new_pos.y, new_pos.z
                                );
                            }
                        }
                    }
                }
                println!("END");
                break; // Only first frame
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
