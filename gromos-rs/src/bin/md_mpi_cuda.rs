//! MPI+CUDA Hybrid Molecular Dynamics
//!
//! This program runs MD simulations using both MPI and CUDA:
//! - MPI for multi-node parallelism
//! - CUDA for GPU acceleration on each node
//!
//! # Architecture
//!
//! ```text
//! Node 0                    Node 1
//! ┌─────────────────────┐  ┌─────────────────────┐
//! │ MPI Rank 0 ──> GPU 0│  │ MPI Rank 2 ──> GPU 0│
//! │ MPI Rank 1 ──> GPU 1│  │ MPI Rank 3 ──> GPU 1│
//! └─────────────────────┘  └─────────────────────┘
//!          ↕ MPI Allreduce ↕
//! ```
//!
//! # Usage
//!
//! ```bash
//! # Build with MPI and CUDA support
//! cargo build --release --features use-mpi,use-cuda --bin md_mpi_cuda
//!
//! # Run on 4 MPI ranks (2 nodes × 2 GPUs/node)
//! mpirun -np 4 \
//!     --map-by ppr:2:node \
//!     ./target/release/md_mpi_cuda \
//!     --topo system.top \
//!     --conf system.g96 \
//!     --steps 10000
//! ```
//!
//! Based on md++/program/md_mpi.cc with CUDA support

use clap::Parser;
use gromos_rs::*;
use mpi::traits::*;

#[derive(Parser, Debug)]
#[command(name = "md_mpi_cuda")]
#[command(about = "MPI+CUDA Hybrid Molecular Dynamics", long_about = None)]
struct Args {
    /// Topology file
    #[arg(long)]
    topo: String,

    /// Configuration file
    #[arg(long)]
    conf: String,

    /// Number of MD steps
    #[arg(long, default_value_t = 1000)]
    steps: usize,

    /// Timestep (ps)
    #[arg(long, default_value_t = 0.002)]
    dt: f64,

    /// Temperature (K)
    #[arg(long, default_value_t = 300.0)]
    temp: f64,

    /// Nonbonded cutoff (nm)
    #[arg(long, default_value_t = 1.4)]
    cutoff: f64,

    /// Output trajectory file
    #[arg(long, default_value = "traj.g96")]
    output: String,

    /// GPU device ID (default: auto-assign based on MPI rank)
    #[arg(long)]
    gpu_id: Option<usize>,

    /// Enable GPU-aware MPI (direct GPU-to-GPU communication)
    #[arg(long, default_value_t = false)]
    gpu_aware_mpi: bool,
}

fn main() {
    // Initialize MPI
    let universe = mpi::initialize().unwrap();
    let world = universe.world();
    let rank = world.rank();
    let size = world.size();

    let args = Args::parse();

    // Master prints header
    if rank == 0 {
        println!("╔════════════════════════════════════════╗");
        println!("║   GROMOS-RS: MPI+CUDA Hybrid MD       ║");
        println!("╚════════════════════════════════════════╝");
        println!();
        println!("MPI Configuration:");
        println!("  Total processes:  {}", size);
        println!(
            "  GPU-aware MPI:    {}",
            if args.gpu_aware_mpi { "Yes" } else { "No" }
        );
        println!();
        println!("Simulation Parameters:");
        println!("  Steps:            {}", args.steps);
        println!("  Timestep:         {} ps", args.dt);
        println!("  Temperature:      {} K", args.temp);
        println!("  Cutoff:           {} nm", args.cutoff);
        println!();
    }

    // Load topology and configuration (all ranks)
    // In real implementation, only master reads and broadcasts
    let topo = load_topology(&args.topo);
    let mut conf = load_configuration(&args.conf);
    let n_atoms = conf.positions.len();

    // Initialize MPI+CUDA force calculator
    let mut gpu_calc = match gpu::MpiGpuForceCalculator::new(&world, n_atoms, args.cutoff) {
        Ok(calc) => calc,
        Err(e) => {
            eprintln!("Rank {}: Failed to initialize GPU: {:?}", rank, e);
            std::process::exit(1);
        },
    };

    // Print GPU assignment
    world.barrier();
    println!(
        "Rank {}/{}: Using GPU {}",
        rank,
        size,
        gpu_calc.gpu_device_id()
    );
    world.barrier();

    if rank == 0 {
        println!("\nStarting MD simulation...\n");
    }

    // Initialize integrator (Leap-Frog)
    let mut integrator = integrator::LeapFrog::new(args.dt);

    // MD loop
    let mut total_energy = 0.0;
    let start_time = std::time::Instant::now();

    for step in 0..args.steps {
        // Calculate forces on GPU
        let (forces, potential_energy) =
            match gpu_calc.calculate_forces_mpi(&world, &conf.positions) {
                Ok(result) => result,
                Err(e) => {
                    eprintln!("Rank {}: GPU force calculation failed: {:?}", rank, e);
                    std::process::exit(1);
                },
            };

        // Integrate equations of motion (on CPU)
        // In full implementation, this could also be on GPU
        integrator.step(&topo, &mut conf, &forces);

        total_energy = potential_energy;

        // Print progress (master only)
        if rank == 0 && step % 100 == 0 {
            let elapsed = start_time.elapsed().as_secs_f64();
            let ns_per_day = (step as f64 * args.dt * 1e-3) / elapsed * 86400.0;
            println!(
                "Step {:6} | Energy: {:12.3} kJ/mol | Perf: {:.2} ns/day",
                step, total_energy, ns_per_day
            );
        }
    }

    // Final synchronization
    world.barrier();

    if rank == 0 {
        let elapsed = start_time.elapsed().as_secs_f64();
        let total_ns = args.steps as f64 * args.dt * 1e-3;
        let ns_per_day = total_ns / elapsed * 86400.0;

        println!("\n╔════════════════════════════════════════╗");
        println!("║         Simulation Complete!          ║");
        println!("╚════════════════════════════════════════╝");
        println!();
        println!("Performance:");
        println!("  Total time:       {:.2} s", elapsed);
        println!("  Simulated time:   {:.3} ns", total_ns);
        println!("  Performance:      {:.2} ns/day", ns_per_day);
        println!("  Final energy:     {:.3} kJ/mol", total_energy);
        println!();
    }
}

// Placeholder functions (in real implementation, these would be proper I/O)
fn load_topology(_path: &str) -> topology::Topology {
    topology::Topology::new()
}

fn load_configuration(_path: &str) -> configuration::Configuration {
    configuration::Configuration::new(100, 1, 1)
}
