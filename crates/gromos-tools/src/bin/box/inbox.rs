//! inbox - Put atoms into the box
//!
//! Usage: inbox @topo <file> @traj <file> @pbc <type>
//!
//! Ensures all atoms are inside the periodic box

use gromos::io::topology::{build_topology, read_topology_file};
use gromos::io::trajectory::TrajectoryReader;
use gromos::math::Vec3;
use std::env;
use std::process;

fn put_in_box(mut pos: Vec3, box_dims: Vec3) -> Vec3 {
    if box_dims.x > 0.0 {
        while pos.x < 0.0 {
            pos.x += box_dims.x;
        }
        while pos.x >= box_dims.x {
            pos.x -= box_dims.x;
        }
    }
    if box_dims.y > 0.0 {
        while pos.y < 0.0 {
            pos.y += box_dims.y;
        }
        while pos.y >= box_dims.y {
            pos.y -= box_dims.y;
        }
    }
    if box_dims.z > 0.0 {
        while pos.z < 0.0 {
            pos.z += box_dims.z;
        }
        while pos.z >= box_dims.z {
            pos.z -= box_dims.z;
        }
    }
    pos
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: inbox @topo <file> @traj <file> @pbc <type>");
        process::exit(1);
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut pbc_type = String::new();

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
            "@pbc" => {
                i += 1;
                pbc_type = args[i].clone();
            },
            _ => {},
        }
        i += 1;
    }

    let topo_data = read_topology_file(&topo_file.unwrap()).unwrap();
    let _topo = build_topology(topo_data);
    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    println!("TITLE");
    println!("Atoms put into box");
    println!("END");

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let box_dims = frame.box_dims;

                println!("TIMESTEP");
                println!("{:15} {:15.4}", frame.step, frame.time);
                println!("END");
                println!("POSITIONRED");

                for (i, pos) in frame.positions.iter().enumerate() {
                    let inbox_pos = if pbc_type != "v" {
                        put_in_box(*pos, box_dims)
                    } else {
                        *pos
                    };

                    println!(
                        "{:6} {:15.9} {:15.9} {:15.9}",
                        i + 1,
                        inbox_pos.x,
                        inbox_pos.y,
                        inbox_pos.z
                    );
                }

                println!("END");
                println!("BOX");
                println!(
                    " {:15.9} {:15.9} {:15.9}",
                    box_dims.x, box_dims.y, box_dims.z
                );
                println!("END");
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }
}
