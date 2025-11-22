//! Integration tests for MPI functionality
//!
//! These tests verify that the MPI implementation correctly distributes work
//! and produces the same results as the serial version.
//!
//! # Running the tests
//!
//! ```bash
//! # Run with 1 process (serial verification)
//! cargo test --features use-mpi --test mpi_tests
//!
//! # Run with multiple MPI processes
//! mpirun -np 4 cargo test --features use-mpi --test mpi_tests -- --test-threads=1
//! ```

#![cfg(feature = "use-mpi")]

use gromos_rs::interaction::nonbonded::*;
use gromos_rs::mpi::{MpiControl, MpiNonbondedMaster};
use gromos_rs::*;

#[test]
fn test_mpi_control_creation() {
    // Test serial MPI control
    let control = MpiControl::serial();
    assert!(!control.enabled);
    assert_eq!(control.rank, 0);
    assert_eq!(control.size, 1);
    assert!(control.is_master());
    assert!(!control.is_slave());
    assert_eq!(control.num_slaves(), 0);
}

#[test]
fn test_mpi_pair_range_division() {
    // Test that pair ranges are divided correctly among processes
    let test_cases = vec![
        (100, 4),  // 100 pairs, 4 processes
        (103, 4),  // 103 pairs, 4 processes (not evenly divisible)
        (1000, 8), // 1000 pairs, 8 processes
        (997, 7),  // 997 pairs, 7 processes
    ];

    for (total_pairs, num_processes) in test_cases {
        let mut ranges = Vec::new();

        for rank in 0..num_processes {
            let control = MpiControl {
                enabled: true,
                rank,
                size: num_processes,
                master_id: 0,
                simulation_id: 0,
            };

            let master = MpiNonbondedMaster::new(control);
            let (start, end) = master.calculate_pair_range(total_pairs);

            ranges.push((start, end));

            // Verify that each range is valid
            assert!(
                start < end,
                "Range must be non-empty: rank {} got [{}, {})",
                rank,
                start,
                end
            );
            assert!(end <= total_pairs, "Range must not exceed total pairs");
        }

        // Verify that ranges are contiguous and cover all pairs
        for i in 0..ranges.len() - 1 {
            assert_eq!(
                ranges[i].1,
                ranges[i + 1].0,
                "Ranges must be contiguous: ranks {} and {} have gap",
                i,
                i + 1
            );
        }

        // Verify that all pairs are covered
        assert_eq!(ranges[0].0, 0, "First range must start at 0");
        assert_eq!(
            ranges[num_processes as usize - 1].1,
            total_pairs,
            "Last range must end at total_pairs"
        );

        // Verify load balance (max difference should be at most 1)
        let sizes: Vec<usize> = ranges.iter().map(|(start, end)| end - start).collect();
        let min_size = *sizes.iter().min().unwrap();
        let max_size = *sizes.iter().max().unwrap();
        assert!(
            max_size - min_size <= 1,
            "Load imbalance too large for {} pairs on {} processes: min={}, max={}",
            total_pairs,
            num_processes,
            min_size,
            max_size
        );
    }
}

#[test]
fn test_force_storage_operations() {
    let n_atoms = 100;
    let mut storage = ForceStorage::new(n_atoms);

    // Test initialization
    assert_eq!(storage.forces.len(), n_atoms);
    assert_eq!(storage.e_lj, 0.0);
    assert_eq!(storage.e_crf, 0.0);

    // Add some forces and energies
    storage.forces[0] = Vec3::new(1.0, 2.0, 3.0);
    storage.e_lj = 100.0;
    storage.e_crf = 50.0;
    storage.virial[0][0] = 10.0;

    // Test clear
    storage.clear();
    assert_eq!(storage.forces[0], Vec3::ZERO);
    assert_eq!(storage.e_lj, 0.0);
    assert_eq!(storage.e_crf, 0.0);
    assert_eq!(storage.virial[0][0], 0.0);
}

/// Test that serial and parallel innerloop produce the same results
#[test]
fn test_serial_vs_parallel_innerloop() {
    use math::{Periodicity, Vec3};

    // Create a small test system
    let n_atoms = 10;
    let positions: Vec<Vec3> = (0..n_atoms)
        .map(|i| {
            Vec3::new(
                (i as f32) * 0.3,
                ((i * 2) as f32) * 0.2,
                ((i * 3) as f32) * 0.1,
            )
        })
        .collect();

    let charges: Vec<f32> = vec![0.5; n_atoms];
    let iac: Vec<u32> = vec![1; n_atoms];

    // Simple pairlist (all pairs)
    let mut pairlist = Vec::new();
    for i in 0..n_atoms {
        for j in (i + 1)..n_atoms {
            pairlist.push((i as u32, j as u32));
        }
    }

    // Create LJ parameters (simple uniform parameters)
    let mut lj_params = vec![
        vec![
            LJParameters {
                c6: 1.0e-3,
                c12: 1.0e-6
            };
            10
        ];
        10
    ];

    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.5 / (2.0 * 1.4_f64.powi(3)),
        crf_cut3i: (1.0 - 0.5 / 2.0) / 1.4,
    };

    let periodicity = Periodicity::rectangular([5.0, 5.0, 5.0]);

    // Calculate with serial innerloop
    let mut serial_storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        &positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &periodicity,
        &mut serial_storage,
    );

    // Calculate with parallel innerloop
    let parallel_storage = lj_crf_innerloop_parallel(
        &positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &periodicity,
        n_atoms,
    );

    // Compare energies (should be identical)
    let energy_tolerance = 1e-10;
    assert!(
        (serial_storage.e_lj - parallel_storage.e_lj).abs() < energy_tolerance,
        "LJ energy mismatch: serial={}, parallel={}",
        serial_storage.e_lj,
        parallel_storage.e_lj
    );
    assert!(
        (serial_storage.e_crf - parallel_storage.e_crf).abs() < energy_tolerance,
        "CRF energy mismatch: serial={}, parallel={}",
        serial_storage.e_crf,
        parallel_storage.e_crf
    );

    // Compare forces
    let force_tolerance = 1e-6;
    for i in 0..n_atoms {
        let diff = (serial_storage.forces[i] - parallel_storage.forces[i]).length();
        assert!(
            diff < force_tolerance,
            "Force mismatch at atom {}: diff={}, serial={:?}, parallel={:?}",
            i,
            diff,
            serial_storage.forces[i],
            parallel_storage.forces[i]
        );
    }

    // Compare virial
    for a in 0..3 {
        for b in 0..3 {
            let diff = (serial_storage.virial[a][b] - parallel_storage.virial[a][b]).abs();
            assert!(
                diff < energy_tolerance,
                "Virial[{}][{}] mismatch: serial={}, parallel={}",
                a,
                b,
                serial_storage.virial[a][b],
                parallel_storage.virial[a][b]
            );
        }
    }
}

#[test]
fn test_lj_crf_interaction_accuracy() {
    use math::Vec3;

    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.5 / (2.0 * 1.4_f64.powi(3)),
        crf_cut3i: (1.0 - 0.5 / 2.0) / 1.4,
    };

    // Test at a few distances
    let test_distances = vec![0.3, 0.5, 0.8, 1.0, 1.2];

    for dist in test_distances {
        let r = Vec3::new(dist, 0.0, 0.0);
        let c6 = 1.0e-3;
        let c12 = 1.0e-6;
        let q_prod = 0.25; // +0.5 * +0.5

        let (f_mag, e_lj, e_crf) = lj_crf_interaction(r, c6, c12, q_prod, &crf);

        // Verify finite values
        assert!(f_mag.is_finite(), "Force should be finite at r={}", dist);
        assert!(e_lj.is_finite(), "LJ energy should be finite at r={}", dist);
        assert!(
            e_crf.is_finite(),
            "CRF energy should be finite at r={}",
            dist
        );

        // Verify that force decreases with distance (at reasonable distances)
        if dist > 0.4 {
            let r2 = Vec3::new(dist * 1.1, 0.0, 0.0);
            let (f_mag2, _, _) = lj_crf_interaction(r2, c6, c12, q_prod, &crf);
            assert!(
                f_mag.abs() > f_mag2.abs(),
                "Force magnitude should decrease with distance: r={}: f={}, r={}: f={}",
                dist,
                f_mag,
                dist * 1.1,
                f_mag2
            );
        }
    }
}

/// Test that pair range calculation handles edge cases
#[test]
fn test_pair_range_edge_cases() {
    // Test with more processes than pairs
    let control = MpiControl {
        enabled: true,
        rank: 5,
        size: 10,
        master_id: 0,
        simulation_id: 0,
    };
    let master = MpiNonbondedMaster::new(control);
    let (start, end) = master.calculate_pair_range(5);

    // Some processes will get 0 pairs
    assert!(end >= start, "Range must be valid even if empty");
    assert!(end <= 5, "Range must not exceed total pairs");

    // Test with 1 pair
    let control = MpiControl {
        enabled: true,
        rank: 0,
        size: 4,
        master_id: 0,
        simulation_id: 0,
    };
    let master = MpiNonbondedMaster::new(control);
    let (start, end) = master.calculate_pair_range(1);
    assert_eq!(start, 0);
    assert_eq!(end, 1);
}

/// Test virial calculation
#[test]
fn test_virial_calculation() {
    use math::{Periodicity, Vec3};

    // Create a simple two-particle system
    let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];

    let charges = vec![0.5, -0.5];
    let iac = vec![1, 1];
    let pairlist = vec![(0, 1)];

    let lj_params = vec![
        vec![
            LJParameters {
                c6: 1.0e-3,
                c12: 1.0e-6
            };
            2
        ];
        2
    ];

    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.5 / (2.0 * 1.4_f64.powi(3)),
        crf_cut3i: (1.0 - 0.5 / 2.0) / 1.4,
    };

    let periodicity = Periodicity::rectangular([5.0, 5.0, 5.0]);

    let mut storage = ForceStorage::new(2);
    lj_crf_innerloop(
        &positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &periodicity,
        &mut storage,
    );

    // Virial should be non-zero (we have interactions)
    let virial_trace = storage.virial[0][0] + storage.virial[1][1] + storage.virial[2][2];
    assert!(
        virial_trace.abs() > 1e-10,
        "Virial trace should be non-zero for interacting particles"
    );

    // For particles along x-axis, virial[0][0] should dominate
    assert!(
        storage.virial[0][0].abs() > storage.virial[1][1].abs(),
        "Virial[0][0] should dominate for x-aligned particles"
    );
}

/// Benchmark-style test to measure relative performance
#[test]
fn test_performance_serial_vs_parallel() {
    use math::{Periodicity, Vec3};
    use std::time::Instant;

    // Create a larger system for meaningful timing
    let n_atoms = 500;
    let positions: Vec<Vec3> = (0..n_atoms)
        .map(|i| {
            let angle = (i as f32) * 0.1;
            Vec3::new(angle.cos() * 5.0, angle.sin() * 5.0, (i as f32) * 0.01)
        })
        .collect();

    let charges: Vec<f32> = vec![0.5; n_atoms];
    let iac: Vec<u32> = vec![1; n_atoms];

    // Create subset of pairs (not all pairs - that would be too slow)
    let mut pairlist = Vec::new();
    for i in 0..n_atoms.min(100) {
        for j in (i + 1)..(i + 20).min(n_atoms) {
            pairlist.push((i as u32, j as u32));
        }
    }

    let lj_params = vec![
        vec![
            LJParameters {
                c6: 1.0e-3,
                c12: 1.0e-6
            };
            10
        ];
        10
    ];
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.5 / (2.0 * 1.4_f64.powi(3)),
        crf_cut3i: (1.0 - 0.5 / 2.0) / 1.4,
    };
    let periodicity = Periodicity::rectangular([10.0, 10.0, 10.0]);

    // Time serial version
    let start = Instant::now();
    let mut serial_storage = ForceStorage::new(n_atoms);
    for _ in 0..5 {
        serial_storage.clear();
        lj_crf_innerloop(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            &mut serial_storage,
        );
    }
    let serial_time = start.elapsed();

    // Time parallel version
    let start = Instant::now();
    for _ in 0..5 {
        let _parallel_storage = lj_crf_innerloop_parallel(
            &positions,
            &charges,
            &iac,
            &pairlist,
            &lj_params,
            &crf,
            &periodicity,
            n_atoms,
        );
    }
    let parallel_time = start.elapsed();

    println!("\n=== Performance Comparison ===");
    println!("System: {} atoms, {} pairs", n_atoms, pairlist.len());
    println!("Serial time:   {:?}", serial_time);
    println!("Parallel time: {:?}", parallel_time);
    println!(
        "Speedup: {:.2}x",
        serial_time.as_secs_f64() / parallel_time.as_secs_f64()
    );

    // Note: Parallel might be slower for small systems due to overhead
    // This is just a sanity check that both complete
    assert!(serial_time.as_secs() < 10, "Serial version taking too long");
    assert!(
        parallel_time.as_secs() < 10,
        "Parallel version taking too long"
    );
}

#[cfg(test)]
mod mpi_integration_tests {
    use super::*;

    /// Test that requires actual MPI execution
    /// Run with: mpirun -np 4 cargo test --features use-mpi test_mpi_actual_run
    #[test]
    #[ignore] // Ignore by default, run explicitly with MPI
    fn test_mpi_actual_run() {
        // This test is meant to be run with mpirun
        if !mpi::is_initialized() {
            println!("MPI not initialized - test requires mpirun");
            return;
        }

        let universe = mpi::initialize().unwrap();
        let world = universe.world();
        let control = MpiControl::new(&world);

        println!("Process {} of {}: Testing MPI", control.rank, control.size);

        // All processes should agree on the size
        assert!(control.size > 0);
        assert!(control.rank < control.size);

        // Only one master
        if control.is_master() {
            assert_eq!(control.rank, control.master_id);
            println!("Master process confirmed");
        } else {
            assert_ne!(control.rank, control.master_id);
            println!("Slave process confirmed");
        }
    }
}
