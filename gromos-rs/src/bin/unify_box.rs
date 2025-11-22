//! unify_box - Unify box shapes in trajectory
//!
//! Usage: unify_box @topo <file> @traj <file> @box <x y z>
//!
//! Sets all frames to use the same box dimensions

use gromos_rs::io::topology::{build_topology, read_topology_file};
use gromos_rs::io::trajectory::TrajectoryReader;
use gromos_rs::math::Vec3;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 9 {
        eprintln!("Usage: unify_box @topo <file> @traj <file> @box <x y z>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut box_x = 0.0;
    let mut box_y = 0.0;
    let mut box_z = 0.0;

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
            "@box" => {
                i += 1;
                box_x = args[i].parse().unwrap();
                i += 1;
                box_y = args[i].parse().unwrap();
                i += 1;
                box_z = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    let unified_box = Vec3::new(box_x, box_y, box_z);

    println!("TITLE");
    println!(
        "Unified box dimensions: {:.3} {:.3} {:.3}",
        box_x, box_y, box_z
    );
    println!("END");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                println!("TIMESTEP");
                println!("{:15} {:15.4}", frame.step, frame.time);
                println!("END");
                println!("POSITIONRED");

                for (i, pos) in frame.positions.iter().enumerate() {
                    println!("{:6} {:15.9} {:15.9} {:15.9}", i + 1, pos.x, pos.y, pos.z);
                }

                println!("END");
                println!("BOX");
                println!(
                    " {:15.9} {:15.9} {:15.9}",
                    unified_box.x, unified_box.y, unified_box.z
                );
                println!("END");
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
