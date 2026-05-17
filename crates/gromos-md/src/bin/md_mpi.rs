//! MPI-parallel molecular dynamics simulation
//!
//! This is the MPI version of the md program, based on md++/program/md_mpi.cc
//!
//! # Usage
//!
//! ```bash
//! # Compile with MPI support
//! cargo build --release --features use-mpi --bin md_mpi
//!
//! # Run with 4 MPI processes
//! mpirun -np 4 ./target/release/md_mpi --topo system.top --conf system.g96 --steps 1000
//! ```

// This binary requires the use-mpi feature
#[cfg(feature = "use-mpi")]
mod enabled {
    use clap::Parser;
    use gromos::*;
    use std::fs::File;
    use std::io::Write;
    use std::path::PathBuf;

    #[derive(Parser, Debug)]
    #[command(name = "md_mpi")]
    #[command(about = "MPI-parallel molecular dynamics simulation", long_about = None)]
    pub struct Args {
        /// Topology file (.top)
        #[arg(short, long)]
        topo: PathBuf,

        /// Configuration file (.g96)
        #[arg(short, long)]
        conf: PathBuf,

        /// Number of MD steps
        #[arg(short, long, default_value_t = 1000)]
        steps: usize,

        /// Time step (ps)
        #[arg(short, long, default_value_t = 0.002)]
        dt: f64,

        /// Output trajectory file
        #[arg(short, long, default_value = "traj.g96")]
        output: PathBuf,

        /// Output frequency (steps)
        #[arg(long, default_value_t = 100)]
        output_freq: usize,

        /// Nonbonded cutoff (nm)
        #[arg(long, default_value_t = 1.4)]
        cutoff: f64,

        /// Enable parallel innerloop (Rayon within each MPI process)
        #[arg(long)]
        parallel_innerloop: bool,

        /// Quiet mode (suppress slave output)
        #[arg(short, long)]
        quiet: bool,
    }

    pub fn run() {
        // Initialize MPI
        let universe = mpi::initialize().unwrap();
        let world = universe.world();
        let rank = world.rank();
        let size = world.size();

        // Create MPI control
        let mpi_control = mpi::MpiControl::new(&world);

        // Setup output streams
        let mut output_file: Option<File> = None;
        if !mpi_control.is_master() {
            // Slaves write to their own log files
            let filename = format!("slave_{}.log", rank);
            output_file = Some(File::create(&filename).expect("Failed to create slave log file"));
        }

        // Parse command line arguments
        let args = Args::parse();

        // Print banner
        if mpi_control.is_master() {
            println!("=======================================================");
            println!("  GROMOS-RS: MPI Molecular Dynamics");
            println!("=======================================================");
            println!("  MPI processes: {}", size);
            println!("  Master rank: {}", mpi_control.master_id);
            println!("  Slave ranks: {}", mpi_control.num_slaves());
            println!("=======================================================");
            println!();
        } else if !args.quiet {
            if let Some(ref mut f) = output_file {
                writeln!(f, "Slave process {} initialized", rank).ok();
            }
        }

        // Load topology and configuration (master only)
        let mut topo = topology::Topology::new();
        let mut conf = configuration::Configuration::new(0, 1, 1);

        if mpi_control.is_master() {
            println!("Loading topology from {:?}...", args.topo);
            topo = io::topology::read_topology(&args.topo).expect("Failed to read topology");

            println!("Loading configuration from {:?}...", args.conf);
            conf = io::g96::read_g96(&args.conf).expect("Failed to read configuration");

            println!("  Atoms: {}", topo.num_atoms());
            println!("  Box: {:?}", conf.current().box_matrix);
            println!();
        }

        // Setup force field parameters
        let _crf_params = interaction::nonbonded::CRFParameters {
            crf_cut: args.cutoff,
            crf_2cut3i: 0.5 / (2.0 * args.cutoff.powi(3)),
            crf_cut3i: (1.0 - 0.5 / 2.0) / args.cutoff,
        };

        // Create integrator
        let mut integrator = if args.parallel_innerloop {
            integrator::LeapFrog::new().with_parallel()
        } else {
            integrator::LeapFrog::new()
        };

        if mpi_control.is_master() {
            println!(
                "Starting MD simulation ({} steps, dt={} ps)...",
                args.steps, args.dt
            );
            println!();
        }

        // Main MD loop
        for step in 0..args.steps {
            // Integration (master only)
            if mpi_control.is_master() {
                integrator.step(args.dt, &topo, &mut conf);
            }

            // TODO: MPI parallel force calculation
            // This requires the full mpi module to be available

            // Output (master only)
            if mpi_control.is_master() && step % args.output_freq == 0 {
                let energy_kinetic = conf.calculate_kinetic_energy(&topo);

                println!(
                    "Step {:6} | T = {:8.3} K | E_kin = {:12.6} kJ/mol",
                    step,
                    conf.calculate_temperature(&topo),
                    energy_kinetic
                );
            }
        }

        // Finalization
        if mpi_control.is_master() {
            println!();
            println!("Simulation completed successfully!");
            println!("Writing final configuration to {:?}...", args.output);

            io::g96::write_g96(&args.output, &conf).expect("Failed to write output");
        }
    }
}

#[cfg(feature = "use-mpi")]
fn main() {
    enabled::run();
}

#[cfg(not(feature = "use-mpi"))]
fn main() {
    eprintln!("ERROR: md_mpi requires the 'use-mpi' feature");
    eprintln!("Compile with: cargo build --features use-mpi --bin md_mpi");
    std::process::exit(1);
}
