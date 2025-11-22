//! gch - Generate constraint history
//!
//! Usage: gch @traj <file> @constraints <file>
//!
//! Generates history of constraint violations

use gromos_rs::io::trajectory::TrajectoryReader;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: gch @traj <file> [@constraints <file>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut const_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@constraints" => {
                i += 1;
                const_file = Some(args[i].clone());
            },
            _ => {},
        }
        i += 1;
    }

    let mut constraints: Vec<(usize, usize, f32)> = Vec::new();

    if let Some(path) = const_file {
        if let Ok(file) = File::open(&path) {
            let reader = BufReader::new(file);
            for line in reader.lines() {
                if let Ok(line_str) = line {
                    let parts: Vec<&str> = line_str.split_whitespace().collect();
                    if parts.len() >= 3 {
                        if let (Ok(i), Ok(j), Ok(d)) = (
                            parts[0].parse::<usize>(),
                            parts[1].parse::<usize>(),
                            parts[2].parse::<f32>(),
                        ) {
                            constraints.push((i - 1, j - 1, d));
                        }
                    }
                }
            }
        }
    }

    if constraints.is_empty() {
        eprintln!("# Warning: No constraints loaded, using default (atoms 1-2, 0.1 nm)");
        constraints.push((0, 1, 0.1));
    }

    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    eprintln!("# Constraint history");
    eprintln!("# Constraints: {}", constraints.len());

    println!("# Frame    Time (ps)    Constraint    Distance (nm)    Violation (nm)");

    let mut frame_idx = 0;
    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                for (cidx, &(i, j, target)) in constraints.iter().enumerate() {
                    if i < frame.positions.len() && j < frame.positions.len() {
                        let dx = frame.positions[i].x - frame.positions[j].x;
                        let dy = frame.positions[i].y - frame.positions[j].y;
                        let dz = frame.positions[i].z - frame.positions[j].z;
                        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                        let viol = (dist - target).abs();

                        println!(
                            "{:8} {:12.4} {:12} {:16.6} {:16.6}",
                            frame_idx + 1,
                            frame.time,
                            cidx + 1,
                            dist,
                            viol
                        );
                    }
                }
                frame_idx += 1;
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    eprintln!("# Processed {} frames", frame_idx);
}
