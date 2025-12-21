//! nhoparam - Calculate Nose-Hoover thermostat parameters
//!
//! Usage: nhoparam @temp <T> @tau <time>
//!
//! Calculates optimal Nose-Hoover parameters

use std::env;
use std::process;

const K_B: f64 = 0.00831446; // kJ/(mol·K)

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        eprintln!("Usage: nhoparam @temp <T> [@tau <ps>] [@ndof <n>]");
        process::exit(1);
    }

    let mut temperature = 300.0;
    let mut tau = 0.1; // ps
    let mut n_dof = 3000;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@temp" => {
                i += 1;
                temperature = args[i].parse().unwrap();
            },
            "@tau" => {
                i += 1;
                tau = args[i].parse().unwrap();
            },
            "@ndof" => {
                i += 1;
                n_dof = args[i].parse().unwrap();
            },
            _ => {},
        }
        i += 1;
    }

    eprintln!("# Nose-Hoover thermostat parameters");
    eprintln!("# Target temperature: {} K", temperature);
    eprintln!("# Relaxation time: {} ps", tau);
    eprintln!("# Degrees of freedom: {}", n_dof);

    // Calculate coupling parameter Q
    let omega = 2.0 * std::f64::consts::PI / tau; // Angular frequency
    let q_mass = n_dof as f64 * K_B * temperature / (omega * omega);

    println!("# Nose-Hoover Parameters");
    println!("# Temperature (K): {}", temperature);
    println!("# Tau (ps): {}", tau);
    println!("# N_DOF: {}", n_dof);
    println!("# Q (thermostat mass): {:.6e}", q_mass);
    println!("# Omega (1/ps): {:.6}", omega);
    println!();
    println!("# GROMOS input block:");
    println!("NOSHOOVER");
    println!("  # Target temperature (K)");
    println!("    {:.2}", temperature);
    println!("  # Coupling time tau (ps)");
    println!("    {:.4}", tau);
    println!("END");
}
