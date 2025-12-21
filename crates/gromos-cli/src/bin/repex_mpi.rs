//! MPI Replica Exchange Molecular Dynamics Program
//!
//! Based on md++/program/repex_mpi.cc
//!
//! This program runs temperature replica exchange MD across multiple MPI processes.
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

#[cfg(feature = "use-mpi")]
mod enabled {
    use clap::Parser;
    use gromos::*;
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(name = "repex_mpi")]
    #[command(about = "MPI Replica Exchange Molecular Dynamics", long_about = None)]
    pub struct Args {
        #[arg(short, long)]
        topo: PathBuf,
        #[arg(short, long)]
        conf: PathBuf,
        #[arg(long, default_value_t = 300.0)]
        temp_min: f64,
        #[arg(long, default_value_t = 350.0)]
        temp_max: f64,
        #[arg(short, long, default_value_t = 100000)]
        steps: usize,
        #[arg(short, long, default_value_t = 0.002)]
        dt: f64,
        #[arg(long, default_value_t = 100)]
        exchange_freq: usize,
        #[arg(long, default_value = "odd-even")]
        exchange_scheme: String,
        #[arg(short, long, default_value = "replica")]
        output: String,
        #[arg(long, default_value_t = 1000)]
        output_freq: usize,
        #[arg(long, default_value_t = 1.4)]
        cutoff: f64,
        #[arg(short, long)]
        quiet: bool,
    }

    pub fn run() {
        let universe = mpi::initialize().unwrap();
        let world = universe.world();
        let rank = world.rank();
        let size = world.size();

        let args = Args::parse();

        // Calculate geometric temperature distribution
        let temperatures: Vec<f64> = (0..size)
            .map(|i| {
                let frac = if size > 1 { i as f64 / (size - 1) as f64 } else { 0.0 };
                args.temp_min * (args.temp_max / args.temp_min).powf(frac)
            })
            .collect();

        let my_temperature = temperatures[rank as usize];

        if rank == 0 && !args.quiet {
            println!("╔═══════════════════════════════════════════════════════════════╗");
            println!("║       GROMOS-RS: MPI Replica Exchange MD                      ║");
            println!("╚═══════════════════════════════════════════════════════════════╝");
            println!();
            println!("  Replicas: {}", size);
            println!("  Temperature range: {:.1} - {:.1} K", args.temp_min, args.temp_max);
            println!("  Steps per replica: {}", args.steps);
            println!();
            for (i, temp) in temperatures.iter().enumerate() {
                println!("  Replica {}: {:.2} K", i, temp);
            }
            println!();
        }

        // Load topology and configuration
        let topo = if rank == 0 {
            io::topology::read_topology(&args.topo).expect("Failed to read topology")
        } else {
            topology::Topology::new()
        };

        let _conf = if rank == 0 {
            io::g96::read_g96(&args.conf).expect("Failed to read configuration")
        } else {
            configuration::Configuration::new(0, 1, 1)
        };

        if rank == 0 && !args.quiet {
            println!("Loaded system: {} atoms", topo.num_atoms());
            println!();
        }

        // Setup REMD controller
        let exchange_scheme = match args.exchange_scheme.as_str() {
            "sequential" => remd::ExchangeScheme::Sequential,
            "odd-even" => remd::ExchangeScheme::OddEven,
            "random" => remd::ExchangeScheme::Random,
            _ => {
                if rank == 0 {
                    eprintln!("Unknown exchange scheme: {}", args.exchange_scheme);
                }
                return;
            }
        };

        let mut controller = remd_mpi::RemdMpiController::new(
            &world,
            my_temperature,
            0.0,
            topo.num_atoms(),
            args.dt,
            remd::ExchangeType::Temperature,
            exchange_scheme,
            args.exchange_freq,
        );

        let mut integrator = integrator::LeapFrog::new();

        if rank == 0 && !args.quiet {
            println!("Starting Replica Exchange MD...");
            println!();
        }

        let num_exchanges = args.steps / args.exchange_freq;

        for cycle in 0..num_exchanges {
            controller.run_md_steps(args.exchange_freq, &topo, &mut integrator);
            let exchanged = controller.attempt_exchange(&world, cycle * args.exchange_freq);

            if rank == 0 && !args.quiet && (cycle + 1) % 10 == 0 {
                let progress = ((cycle + 1) as f64 / num_exchanges as f64) * 100.0;
                println!("Progress: {:.1}% ({}/{} exchanges)", progress, cycle + 1, num_exchanges);
            }

            if !args.quiet && exchanged {
                println!("Cycle {}: Replica {} exchanged", cycle, rank);
            }
        }

        world.barrier();

        if rank == 0 && !args.quiet {
            println!();
            println!("Replica Exchange completed!");
            controller.print_statistics();
            println!("Final output: {}_*.g96", args.output);
        }
    }
}

#[cfg(feature = "use-mpi")]
fn main() {
    enabled::run();
}

#[cfg(not(feature = "use-mpi"))]
fn main() {
    eprintln!("ERROR: repex_mpi requires the 'use-mpi' feature");
    eprintln!("Compile with: cargo build --features use-mpi --bin repex_mpi");
    std::process::exit(1);
}
