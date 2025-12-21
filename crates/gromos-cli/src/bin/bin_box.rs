//! bin_box - Bin coordinates into grid cells
//!
//! Usage: bin_box @topo <file> @traj <file> @bins <nx ny nz>
//!
//! Divides simulation box into bins and counts atoms per bin

use gromos::io::topology::{build_topology, read_topology_file};
use gromos::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 9 {
        eprintln!("Usage: bin_box @topo <file> @traj <file> @bins <nx ny nz>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut nx = 10;
    let mut ny = 10;
    let mut nz = 10;

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
            "@bins" => {
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
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    eprintln!("# Binning into {} x {} x {} grid", nx, ny, nz);

    println!("# Time (ps)    Bin_X  Bin_Y  Bin_Z  Count");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let box_dims = frame.box_dims;
                let dx = box_dims.x / nx as f32;
                let dy = box_dims.y / ny as f32;
                let dz = box_dims.z / nz as f32;

                // Count atoms in each bin
                let mut bins = vec![vec![vec![0; nz]; ny]; nx];

                for pos in &frame.positions {
                    let ix = ((pos.x / dx).floor() as usize).min(nx - 1);
                    let iy = ((pos.y / dy).floor() as usize).min(ny - 1);
                    let iz = ((pos.z / dz).floor() as usize).min(nz - 1);
                    bins[ix][iy][iz] += 1;
                }

                // Print non-empty bins
                for ix in 0..nx {
                    for iy in 0..ny {
                        for iz in 0..nz {
                            if bins[ix][iy][iz] > 0 {
                                println!(
                                    "{:12.4} {:6} {:6} {:6} {:6}",
                                    frame.time, ix, iy, iz, bins[ix][iy][iz]
                                );
                            }
                        }
                    }
                }
                break; // Only first frame
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
