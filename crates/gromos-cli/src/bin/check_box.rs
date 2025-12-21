//! check_box - Check simulation box dimensions
//!
//! Usage: check_box @traj <file>
//!
//! Validates box dimensions throughout trajectory

use gromos::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: check_box @traj <file>");
        process::exit(1);
    }

    let mut traj_file = None;

    let mut i = 1;
    while i < args.len() {
        if args[i] == "@traj" {
            i += 1;
            traj_file = Some(args[i].clone());
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    println!("# Frame    Time (ps)    Box_X (nm)    Box_Y (nm)    Box_Z (nm)    Volume (nm³)");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let box_dims = frame.box_dims;
                let volume = box_dims.x * box_dims.y * box_dims.z;

                println!(
                    "{:8} {:12.4} {:12.6} {:12.6} {:12.6} {:12.4}",
                    frame_idx + 1,
                    frame.time,
                    box_dims.x,
                    box_dims.y,
                    box_dims.z,
                    volume
                );
                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Checked {} frames", frame_idx);
}
