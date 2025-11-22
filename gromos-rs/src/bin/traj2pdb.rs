//! traj2pdb - Convert trajectory to PDB format
//!
//! Usage: traj2pdb @traj <file> @frame <n>
//!
//! Converts trajectory frame to PDB format

use gromos_rs::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: traj2pdb @traj <file> [@frame <n>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut target_frame = 1;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@frame" => {
                i += 1;
                target_frame = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    eprintln!("# Converting trajectory to PDB");
    eprintln!("# Target frame: {}", target_frame);

    let mut frame_idx = 0;

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                frame_idx += 1;

                if frame_idx == target_frame {
                    println!("REMARK   Frame {} at time {:.2} ps", frame_idx, frame.time);
                    println!(
                        "CRYST1{:9.3}{:9.3}{:9.3}  90.00  90.00  90.00 P 1           1",
                        frame.box_dims.x * 10.0, // Convert nm to Å
                        frame.box_dims.y * 10.0,
                        frame.box_dims.z * 10.0
                    );

                    for (idx, pos) in frame.positions.iter().enumerate() {
                        println!("ATOM  {:5}  CA  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C",
                                 idx + 1,
                                 (idx / 100) + 1,
                                 pos.x * 10.0,  // Convert nm to Å
                                 pos.y * 10.0,
                                 pos.z * 10.0);
                    }

                    println!("END");
                    break;
                }
            },
            Ok(None) => {
                if frame_idx < target_frame {
                    eprintln!(
                        "Warning: Only {} frames available, requested frame {}",
                        frame_idx, target_frame
                    );
                }
                break;
            },
            Err(_) => break,
        }
    }

    eprintln!("# Conversion complete");
}
