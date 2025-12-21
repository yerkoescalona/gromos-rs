//! tser - Time series analysis
//!
//! Usage: tser @traj <file> @property <name>
//!
//! Analyzes time series properties from trajectories

use gromos::io::trajectory::TrajectoryReader;
use std::env;
use std::process;

fn calc_mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn calc_std(values: &[f64], mean: f64) -> f64 {
    let variance = values.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / values.len() as f64;
    variance.sqrt()
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: tser @traj <file> [@property <volume|density>]");
        process::exit(1);
    }

    let mut traj_file = None;
    let mut property = "volume";

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "@property" => {
                i += 1;
                property = &args[i];
            },
            _ => {},
        }
        i += 1;
    }

    let mut traj = TrajectoryReader::new(&traj_file.unwrap()).unwrap();

    eprintln!("# Time series analysis");
    eprintln!("# Property: {}", property);

    let mut values = Vec::new();
    let mut times = Vec::new();

    loop {
        match traj.read_frame() {
            Ok(Some(frame)) => {
                let value = match property {
                    "density" => {
                        let volume = frame.box_dims.x * frame.box_dims.y * frame.box_dims.z;
                        if volume > 0.0 {
                            (frame.positions.len() as f64) / (volume as f64)
                        } else {
                            0.0
                        }
                    },
                    _ => {
                        // volume
                        (frame.box_dims.x * frame.box_dims.y * frame.box_dims.z) as f64
                    },
                };

                values.push(value);
                times.push(frame.time as f64);
            },
            Ok(None) => break,
            Err(_) => break,
        }
    }

    if values.is_empty() {
        eprintln!("Error: No data points");
        process::exit(1);
    }

    let mean = calc_mean(&values);
    let std = calc_std(&values, mean);
    let min = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    eprintln!("# Analyzed {} frames", values.len());

    println!("# Time Series Statistics");
    println!("# Property: {}", property);
    println!("# Points: {}", values.len());
    println!("# Mean: {:.6}", mean);
    println!("# Std Dev: {:.6}", std);
    println!("# Min: {:.6}", min);
    println!("# Max: {:.6}", max);
    println!();
    println!("# Time (ps)    Value");

    for (time, value) in times.iter().zip(values.iter()) {
        println!("{:12.4} {:15.6}", time, value);
    }
}
