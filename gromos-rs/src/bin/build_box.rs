//! build_box - Build simulation box with molecules
//!
//! Usage: build_box @solute <file> @box <x y z> @density <rho>
//!
//! Builds simulation box with specified density

use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 7 {
        eprintln!("Usage: build_box @solute <file> @box <x y z> [@density <rho>]");
        process::exit(1);
    }

    let mut _solute_file = None;
    let mut box_x = 5.0f32;
    let mut box_y = 5.0f32;
    let mut box_z = 5.0f32;
    let mut density = 1000.0f64; // kg/m³

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@solute" => {
                i += 1;
                _solute_file = Some(args[i].clone());
            },
            "@box" => {
                i += 1;
                box_x = args[i].parse().unwrap();
                i += 1;
                box_y = args[i].parse().unwrap();
                i += 1;
                box_z = args[i].parse().unwrap();
            },
            "@density" => {
                i += 1;
                density = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    let volume = (box_x * box_y * box_z) as f64; // nm³
    let volume_m3 = volume * 1e-27; // Convert to m³
    let total_mass_kg = density * volume_m3;

    // Estimate number of molecules (assuming water, 18 g/mol)
    let molar_mass = 18.0; // g/mol
    let na = 6.022e23;
    let n_molecules = (total_mass_kg * 1000.0 / molar_mass * na) as usize;

    eprintln!("# Building simulation box");
    eprintln!("# Box: {} x {} x {} nm", box_x, box_y, box_z);
    eprintln!("# Volume: {:.4} nm³", volume);
    eprintln!("# Target density: {} kg/m³", density);
    eprintln!("# Estimated molecules: {}", n_molecules);

    println!("TITLE");
    println!("  Built simulation box");
    println!("END");
    println!("POSITION");
    println!("  {:15.9} {:15.9} {:15.9}", box_x, box_y, box_z);

    // Simple grid placement
    let n_side = (n_molecules as f64).cbrt().ceil() as usize;
    let spacing_x = box_x / n_side as f32;
    let spacing_y = box_y / n_side as f32;
    let spacing_z = box_z / n_side as f32;

    let mut count = 0;
    for ix in 0..n_side {
        for iy in 0..n_side {
            for iz in 0..n_side {
                if count >= n_molecules {
                    break;
                }

                let x = (ix as f32 + 0.5) * spacing_x;
                let y = (iy as f32 + 0.5) * spacing_y;
                let z = (iz as f32 + 0.5) * spacing_z;

                // Place 3-atom molecule (e.g., water)
                println!("{:15.9} {:15.9} {:15.9}", x, y, z);
                println!("{:15.9} {:15.9} {:15.9}", x + 0.01, y, z);
                println!("{:15.9} {:15.9} {:15.9}", x, y + 0.01, z);

                count += 1;
            }
        }
    }

    println!("END");

    eprintln!("# Placed {} molecules ({} atoms)", count, count * 3);
}
