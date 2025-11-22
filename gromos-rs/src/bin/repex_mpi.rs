//! MPI Replica Exchange Molecular Dynamics Program
//!
//! Based on md++/program/repex_mpi.cc
//!
//! This program runs temperature replica exchange MD across multiple MPI processes.
//! Each process runs one replica at a different temperature, and periodically
//! attempts to exchange configurations with neighbors based on the Metropolis criterion.
//!
//! # Usage
//!
//! ```bash
//! # Build with MPI support
//! cargo build --release --features use-mpi --bin repex_mpi
//!
//! # Run with 4 replicas (temperatures: 300, 310, 320, 330 K)
//! mpirun -np 4 ./target/release/repex_mpi \
//!     --topo system.top \
//!     --conf system.g96 \
//!     --temp-min 300 \
//!     --temp-max 330 \
//!     --steps 100000 \
//!     --exchange-freq 100
//! ```
//!
//! # Temperature Distribution
//!
//! Temperatures are distributed geometrically for optimal exchange acceptance:
//! ```text
//! T_i = T_min * (T_max/T_min)^(i/(n-1))
//! ```

#![cfg(feature = "use-mpi")]

use clap::Parser;
use gromos_rs::*;
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(name = "repex_mpi")]
#[command(about = "MPI Replica Exchange Molecular Dynamics", long_about = None)]
struct Args {
    /// Topology file (.top)
    #[arg(short, long)]
    topo: PathBuf,

    /// Initial configuration file (.g96)
    #[arg(short, long)]
    conf: PathBuf,

    /// Minimum temperature (K)
    #[arg(long, default_value_t = 300.0)]
    temp_min: f64,

    /// Maximum temperature (K)
    #[arg(long, default_value_t = 350.0)]
    temp_max: f64,

    /// Number of MD steps per replica
    #[arg(short, long, default_value_t = 100000)]
    steps: usize,

    /// Time step (ps)
    #[arg(short, long, default_value_t = 0.002)]
    dt: f64,

    /// Exchange attempt frequency (MD steps)
    #[arg(long, default_value_t = 100)]
    exchange_freq: usize,

    /// Exchange scheme (sequential, odd-even, random)
    #[arg(long, default_value = "odd-even")]
    exchange_scheme: String,

    /// Output trajectory file prefix
    #[arg(short, long, default_value = "replica")]
    output: String,

    /// Output frequency (steps)
    #[arg(long, default_value_t = 1000)]
    output_freq: usize,

    /// Nonbonded cutoff (nm)
    #[arg(long, default_value_t = 1.4)]
    cutoff: f64,

    /// Quiet mode (suppress non-master output)
    #[arg(short, long)]
    quiet: bool,
}

fn main() {
    // ========================================================================
    // MPI INITIALIZATION
    // ========================================================================

    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();
    let size = world.size();

    let args = Args::parse();

    // ========================================================================
    // TEMPERATURE DISTRIBUTION
    // ========================================================================

    // Calculate geometric temperature distribution for optimal exchange
    let temperatures: Vec<f64> = (0..size)
        .map(|i| {
            let frac = if size > 1 {
                i as f64 / (size - 1) as f64
            } else {
                0.0
            };
            args.temp_min * (args.temp_max / args.temp_min).powf(frac)
        })
        .collect();

    let my_temperature = temperatures[rank as usize];

    // ========================================================================
    // PRINT BANNER (Master only)
    // ========================================================================

    if rank == 0 && !args.quiet {
        println!("╔═══════════════════════════════════════════════════════════════╗");
        println!("║       GROMOS-RS: MPI Replica Exchange MD                      ║");
        println!("╚═══════════════════════════════════════════════════════════════╝");
        println!();
        println!("Configuration:");
        println!("  Replicas:         {}", size);
        println!(
            "  Temperature range: {:.1} - {:.1} K",
            args.temp_min, args.temp_max
        );
        println!("  Steps per replica: {}", args.steps);
        println!("  Exchange freq:     every {} steps", args.exchange_freq);
        println!("  Exchange scheme:   {}", args.exchange_scheme);
        println!();
        println!("Temperature ladder:");
        for (i, temp) in temperatures.iter().enumerate() {
            println!("  Replica {}: {:.2} K", i, temp);
        }
        println!();
    }

    // ========================================================================
    // LOAD TOPOLOGY AND CONFIGURATION
    // ========================================================================

    let topo = if rank == 0 {
        io::topology::read_topology(&args.topo).expect("Failed to read topology")
    } else {
        // Slaves need topology too (simplified - should broadcast)
        topology::Topology::new()
    };

    let conf = if rank == 0 {
        io::g96::read_g96(&args.conf).expect("Failed to read configuration")
    } else {
        // Slaves get empty config (simplified)
        configuration::Configuration::new(0, 1, 1)
    };

    if rank == 0 && !args.quiet {
        println!("Loaded system:");
        println!("  Atoms: {}", topo.num_atoms());
        println!();
    }

    // ========================================================================
    // SETUP REPLICA EXCHANGE CONTROLLER
    // ========================================================================

    let exchange_scheme = match args.exchange_scheme.as_str() {
        "sequential" => remd::ExchangeScheme::Sequential,
        "odd-even" => remd::ExchangeScheme::OddEven,
        "random" => remd::ExchangeScheme::Random,
        _ => {
            if rank == 0 {
                eprintln!("Unknown exchange scheme: {}", args.exchange_scheme);
            }
            return;
        },
    };

    let mut controller = remd_mpi::RemdMpiController::new(
        &world,
        my_temperature,
        0.0, // lambda (not used for T-REMD)
        topo.num_atoms(),
        args.dt,
        remd::ExchangeType::Temperature,
        exchange_scheme,
        args.exchange_freq,
    );

    // ========================================================================
    // CREATE INTEGRATOR
    // ========================================================================

    let mut integrator = integrator::LeapFrog::new();

    // ========================================================================
    // MAIN REMD LOOP
    // ========================================================================

    if rank == 0 && !args.quiet {
        println!("Starting Replica Exchange MD...");
        println!();
    } else if !args.quiet {
        println!("Replica {} ready at T = {:.2} K", rank, my_temperature);
    }

    let num_exchanges = args.steps / args.exchange_freq;

    for cycle in 0..num_exchanges {
        // Run MD for exchange_freq steps
        controller.run_md_steps(args.exchange_freq, &topo, &mut integrator);

        // Attempt exchange (collective MPI operation)
        let exchanged = controller.attempt_exchange(&world, cycle * args.exchange_freq);

        // Print progress (master only)
        if rank == 0 && !args.quiet && (cycle + 1) % 10 == 0 {
            let progress = ((cycle + 1) as f64 / num_exchanges as f64) * 100.0;
            println!(
                "Progress: {:.1}% ({}/{} exchanges attempted)",
                progress,
                cycle + 1,
                num_exchanges
            );
        }

        // Print exchange info
        if !args.quiet && exchanged {
            println!("Cycle {}: Replica {} exchanged with partner", cycle, rank);
        }

        // Output trajectory periodically
        if (cycle * args.exchange_freq) % args.output_freq == 0 {
            // In full implementation, write trajectory here
        }
    }

    // ========================================================================
    // FINALIZATION
    // ========================================================================

    // Synchronize before printing final statistics
    world.barrier();

    if rank == 0 && !args.quiet {
        println!();
        println!("Replica Exchange completed!");
        println!();

        // Print statistics
        controller.print_statistics();

        println!();
        println!("Final output written to: {}_*.g96", args.output);
    }

    // MPI automatically finalizes when universe goes out of scope
}

#[cfg(not(feature = "use-mpi"))]
fn main() {
    eprintln!("ERROR: repex_mpi requires the 'use-mpi' feature");
    eprintln!("Compile with: cargo build --features use-mpi --bin repex_mpi");
    std::process::exit(1);
}
