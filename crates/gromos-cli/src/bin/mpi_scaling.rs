//! MPI Scaling Benchmark Tool
//!
//! This program measures the parallel efficiency and scaling behavior of the
//! MPI implementation across different numbers of processes.
//!
//! # Usage
//!
//! ```bash
//! # Build
//! cargo build --release --features use-mpi --bin mpi_scaling
//!
//! # Run scaling tests with 1, 2, 4, 8 processes
//! for np in 1 2 4 8; do
//!     echo "=== Testing with $np processes ==="
//!     mpirun -np $np ./target/release/mpi_scaling --atoms 1000 --steps 100
//! done
//! ```

// This binary requires the use-mpi feature
#[cfg(feature = "use-mpi")]
mod enabled {
    use clap::Parser;
    use gromos::*;
    use std::time::Instant;

    #[derive(Parser, Debug)]
    #[command(name = "mpi_scaling")]
    #[command(about = "MPI scaling benchmark for GROMOS-RS", long_about = None)]
    pub struct Args {
        /// Number of atoms in test system
        #[arg(short, long, default_value_t = 1000)]
        atoms: usize,

        /// Number of MD steps
        #[arg(short, long, default_value_t = 100)]
        steps: usize,

        /// Nonbonded cutoff (nm)
        #[arg(long, default_value_t = 1.4)]
        cutoff: f64,

        /// Box size (nm)
        #[arg(long, default_value_t = 5.0)]
        box_size: f32,

        /// Print detailed timing for each step
        #[arg(long)]
        detailed: bool,

        /// Output CSV format
        #[arg(long)]
        csv: bool,

        /// Warmup steps (not counted in timing)
        #[arg(long, default_value_t = 10)]
        warmup: usize,
    }

    fn create_test_system(
        n_atoms: usize,
        box_size: f32,
    ) -> (topology::Topology, configuration::Configuration) {
        use math::Vec3;

        let mut topo = topology::Topology::new();
        topo.num_solute_atoms = n_atoms;
        topo.num_atoms_cache = Some(n_atoms);

        // Simple atom properties
        topo.mass = vec![16.0; n_atoms]; // Like oxygen
        topo.inverse_mass = vec![1.0 / 16.0; n_atoms];
        topo.charge = vec![0.0; n_atoms]; // Neutral for simplicity
        topo.iac = vec![1; n_atoms]; // All same type

        // LJ parameters (water-like)
        let lj_param = topology::LJParameters {
            c6: 2.617e-3,  // kJ mol^-1 nm^6
            c12: 2.634e-6, // kJ mol^-1 nm^12
        };
        topo.lj_params = vec![vec![lj_param; 10]; 10];

        // Create configuration with random positions
        let mut conf = configuration::Configuration::new(n_atoms, 1, 1);

        // Random positions in box
        use rand::Rng;
        let mut rng = rand::thread_rng();
        for i in 0..n_atoms {
            conf.current_mut().pos[i] = Vec3::new(
                rng.gen::<f32>() * box_size,
                rng.gen::<f32>() * box_size,
                rng.gen::<f32>() * box_size,
            );
            conf.current_mut().vel[i] = Vec3::new(
                rng.gen::<f32>() - 0.5,
                rng.gen::<f32>() - 0.5,
                rng.gen::<f32>() - 0.5,
            );
        }

        // Set box
        conf.current_mut().box_matrix = [
            [box_size, 0.0, 0.0],
            [0.0, box_size, 0.0],
            [0.0, 0.0, box_size],
        ];

        (topo, conf)
    }

    fn generate_pairlist(
        n_atoms: usize,
        cutoff: f64,
        positions: &[math::Vec3],
        box_size: f32,
    ) -> Vec<(u32, u32)> {
        let mut pairlist = Vec::new();
        let cutoff_sq = (cutoff * cutoff) as f32;

        // Simple all-pairs with cutoff check
        for i in 0..n_atoms {
            for j in (i + 1)..n_atoms {
                let mut r = positions[j] - positions[i];

                // Apply minimum image convention
                if r.x > box_size / 2.0 {
                    r.x -= box_size;
                }
                if r.x < -box_size / 2.0 {
                    r.x += box_size;
                }
                if r.y > box_size / 2.0 {
                    r.y -= box_size;
                }
                if r.y < -box_size / 2.0 {
                    r.y += box_size;
                }
                if r.z > box_size / 2.0 {
                    r.z -= box_size;
                }
                if r.z < -box_size / 2.0 {
                    r.z += box_size;
                }

                if r.length_squared() < cutoff_sq {
                    pairlist.push((i as u32, j as u32));
                }
            }
        }

        pairlist
    }

    pub fn run() {
        // Initialize MPI
        let universe = mpi::initialize().unwrap();
        let world = universe.world();
        let rank = world.rank();
        let size = world.size();

        let mpi_control = mpi::MpiControl::new(&world);

        let args = Args::parse();

        // Master prints header
        if mpi_control.is_master() {
            if args.csv {
                println!("ranks,atoms,steps,time_s,time_per_step_ms,pairs");
            } else {
                println!("╔════════════════════════════════════════╗");
                println!("║      MPI Scaling Benchmark Tool       ║");
                println!("╚════════════════════════════════════════╝");
                println!();
                println!("Configuration:");
                println!("  MPI ranks:    {}", size);
                println!("  Atoms:        {}", args.atoms);
                println!("  Steps:        {}", args.steps);
                println!("  Warmup:       {}", args.warmup);
                println!("  Cutoff:       {} nm", args.cutoff);
                println!("  Box size:     {} nm", args.box_size);
                println!();
            }
        }

        // Create test system
        let (topo, mut conf) = create_test_system(args.atoms, args.box_size);

        // Setup force calculator with MPI
        let mut force_calc = mpi::MpiForceCalculator::new(&world, args.atoms, args.cutoff);

        // Generate initial pairlist
        let positions = &conf.current().pos;
        let pairlist = generate_pairlist(args.atoms, args.cutoff, positions, args.box_size);
        let n_pairs = pairlist.len();

        // Warmup
        for _ in 0..args.warmup {
            let _ = force_calc.calculate_forces(&conf.current().pos, &topo);
        }
        mpi_control.barrier();

        // Timing run
        let start = Instant::now();
        let mut total_energy = 0.0;

        for step in 0..args.steps {
            let (forces, energy) = force_calc.calculate_forces(&conf.current().pos, &topo);
            total_energy = energy;

            // Apply forces (simple integration)
            let dt: f32 = 0.001;
            for i in 0..args.atoms {
                let inv_mass = topo.inverse_mass[i] as f32;
                conf.current_mut().vel[i] += forces[i] * (inv_mass * dt);
                conf.current_mut().pos[i] += conf.current().vel[i] * dt;
            }

            if args.detailed && mpi_control.is_master() && step % 10 == 0 {
                println!("  Step {:4}: E = {:12.3} kJ/mol", step, total_energy);
            }
        }

        mpi_control.barrier();
        let elapsed = start.elapsed();

        // Print results
        if mpi_control.is_master() {
            let time_s = elapsed.as_secs_f64();
            let time_per_step_ms = time_s * 1000.0 / args.steps as f64;

            if args.csv {
                println!(
                    "{},{},{},{:.6},{:.4},{}",
                    size, args.atoms, args.steps, time_s, time_per_step_ms, n_pairs
                );
            } else {
                println!("Results:");
                println!("  Total time:      {:.3} s", time_s);
                println!("  Time per step:   {:.4} ms", time_per_step_ms);
                println!("  Pairs evaluated: {}", n_pairs);
                println!("  Final energy:    {:.3} kJ/mol", total_energy);
                println!();
                println!("To calculate speedup and efficiency:");
                println!("  1. Run with 1 process and note the time");
                println!("  2. Run with N processes");
                println!("  3. Speedup = T(1) / T(N)");
                println!("  4. Efficiency = Speedup / N");
                println!();
            }
        }
    }
}

#[cfg(feature = "use-mpi")]
fn main() {
    enabled::run();
}

#[cfg(not(feature = "use-mpi"))]
fn main() {
    eprintln!("ERROR: mpi_scaling requires the 'use-mpi' feature");
    eprintln!("Compile with: cargo build --features use-mpi --bin mpi_scaling");
    std::process::exit(1);
}
