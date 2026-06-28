//! Integration tests validating gromos-rs against GROMOS reference data
//!
//! These tests parse the reference topology/coordinates/parameters,
//! compute forces using the gromos-rs force kernels, and compare
//! against the expected forces/energies from GROMOS (double precision).
//!
//! Reference data is in GROMOS_references/ and must NOT be modified.
//!
//! Test levels (simplest first):
//!   Level 0: pair_lj, nacl_pair, pair_lj_mixed (2 atoms, vacuum, no PBC)
//!   Level 1: water_single, aladip_vacuum, benzene_vacuum
//!   Level 2: water_3_box, nacl_water_box
//!   Level 3: water_216_box, water_216_nvt, water_216_npt
//!   Level 4: aladip_solvated

use gromos::forces::nonbonded::{
    lj_crf_innerloop, lj_crf_interaction, CRFParameters, ForceStorage, LJParamMatrix,
};
use gromos::io::coordinate::read_coordinate_file;
use gromos::io::topology::{build_topology, read_topology_file};
use gromos::math::Vec3;
use gromos::topology::LJParameters;

use std::path::{Path, PathBuf};

// Initialize env_logger once for all tests (safe to call multiple times)
fn init_logging() {
    let _ = env_logger::try_init();
}

// ─── Helpers ────────────────────────────────────────────────────────────────

/// Path to GROMOS_references (in gromos-md/tests/)
fn ref_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("gromos-md/tests/GROMOS_references")
}

/// Parse the FREEFORCERED block from forces.trf at a given step
fn parse_forces_trf(path: &Path, step: usize) -> Vec<Vec3> {
    let content = std::fs::read_to_string(path).expect("Cannot read forces.trf");
    let mut forces = Vec::new();
    let mut current_step = None;
    let mut in_freeforce = false;

    for line in content.lines() {
        let trimmed = line.trim();

        if trimmed == "TIMESTEP" {
            // Next non-comment/non-empty line has step number
            current_step = None;
            continue;
        }

        if current_step.is_none()
            && trimmed != "END"
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
        {
            // Parse step number (first field on the line)
            if let Some(step_str) = trimmed.split_whitespace().next() {
                if let Ok(s) = step_str.parse::<usize>() {
                    current_step = Some(s);
                }
            }
            continue;
        }

        if trimmed == "FREEFORCERED" {
            if current_step == Some(step) {
                in_freeforce = true;
                forces.clear();
            }
            continue;
        }

        if trimmed == "END" {
            if in_freeforce {
                return forces;
            }
            continue;
        }

        if in_freeforce {
            let parts: Vec<f64> = trimmed
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            if parts.len() == 3 {
                forces.push(Vec3::new(parts[0], parts[1], parts[2]));
            }
        }
    }

    forces
}

/// Parse ENERGY03 block from energies.tre at a given step
/// Returns (e_total, e_kinetic, e_potential, e_vdw, e_crf)
fn parse_energies_tre(path: &Path, step: usize) -> EnergyValues {
    let content = std::fs::read_to_string(path).expect("Cannot read energies.tre");
    let mut current_step = None;
    let mut in_energy = false;
    let mut energy_values: Vec<f64> = Vec::new();

    for line in content.lines() {
        let trimmed = line.trim();

        if trimmed == "TIMESTEP" {
            current_step = None;
            continue;
        }

        if current_step.is_none()
            && trimmed != "END"
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && !trimmed.starts_with("ENEVERSION")
            && !trimmed.starts_with("2023")
        {
            if let Some(step_str) = trimmed.split_whitespace().next() {
                if let Ok(s) = step_str.parse::<usize>() {
                    current_step = Some(s);
                }
            }
            continue;
        }

        if trimmed == "ENERGY03" {
            if current_step == Some(step) {
                in_energy = true;
                energy_values.clear();
            }
            continue;
        }

        if trimmed == "END" {
            if in_energy {
                break;
            }
            continue;
        }

        if in_energy && !trimmed.starts_with('#') && !trimmed.is_empty() {
            if let Ok(val) = trimmed.parse::<f64>() {
                energy_values.push(val);
            }
        }
    }

    // ENERGY03 format (GROMOS standard):
    // [0] = E_total, [1] = E_kinetic, [2] = E_potential
    // ... more fields follow. VdW and CRF positions depend on the system.
    EnergyValues {
        e_total: energy_values.get(0).copied().unwrap_or(0.0),
        e_kinetic: energy_values.get(1).copied().unwrap_or(0.0),
        e_potential: energy_values.get(2).copied().unwrap_or(0.0),
    }
}

#[derive(Debug)]
struct EnergyValues {
    e_total: f64,
    e_kinetic: f64,
    e_potential: f64,
}

/// Convert topology LJ parameters to the nonbonded format
fn convert_lj_params(topo_lj: &[Vec<LJParameters>]) -> LJParamMatrix {
    let nested: Vec<Vec<gromos::forces::nonbonded::LJParameters>> = topo_lj
        .iter()
        .map(|row| {
            row.iter()
                .map(|p| gromos::forces::nonbonded::LJParameters {
                    c6: p.c6,
                    c12: p.c12,
                    cs6: p.cs6,
                    cs12: p.cs12,
                })
                .collect()
        })
        .collect();
    LJParamMatrix::from_nested(&nested)
}

/// Assert two Vec3 are approximately equal (absolute tolerance)
fn assert_vec3_approx(actual: Vec3, expected: Vec3, tol: f64, msg: &str) {
    let dx = (actual.x - expected.x).abs();
    let dy = (actual.y - expected.y).abs();
    let dz = (actual.z - expected.z).abs();
    assert!(
        dx < tol && dy < tol && dz < tol,
        "{}: expected ({:.9e}, {:.9e}, {:.9e}), got ({:.9e}, {:.9e}, {:.9e}), max_diff={:.2e}",
        msg,
        expected.x,
        expected.y,
        expected.z,
        actual.x,
        actual.y,
        actual.z,
        dx.max(dy).max(dz)
    );
}

/// Assert f64 values are approximately equal (relative tolerance)
fn assert_f64_rel(actual: f64, expected: f64, tol: f64, msg: &str) {
    let diff = (actual - expected).abs();
    let scale = expected.abs().max(1e-15);
    let rel = diff / scale;
    assert!(
        rel < tol,
        "{}: expected {:.10e}, got {:.10e}, rel_diff={:.2e}",
        msg,
        expected,
        actual,
        rel
    );
}

// ─── Level 0: pair_lj ───────────────────────────────────────────────────────

#[test]
fn test_pair_lj_forces_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj");

    // Read topology
    let parsed_topo = read_topology_file(base.join("pair_lj.topo")).unwrap();
    let topo = build_topology(parsed_topo);

    // Read coordinates
    let conf = read_coordinate_file(base.join("pair_lj.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    // Set up force calculation
    let n_atoms = positions.len();
    assert_eq!(n_atoms, 2);

    // Build pairlist: all pairs (no exclusions for this system)
    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];

    // IAC and charges
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();

    // LJ parameters
    let lj_params = convert_lj_params(&topo.lj_parameters);

    // CRF parameters (no electrostatics needed, charges are 0)
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    // Compute forces
    let bc = gromos::core::math::Vacuum;
    let mut storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        &mut storage,
    );

    // Load expected forces
    let expected_forces = parse_forces_trf(&base.join("expected/forces.trf"), 0);
    assert_eq!(expected_forces.len(), n_atoms);

    // Compare
    let tol = 1e-6; // kJ/(mol*nm)
    for i in 0..n_atoms {
        assert_vec3_approx(
            storage.forces[i],
            expected_forces[i],
            tol,
            &format!("pair_lj force atom {}", i),
        );
    }
}

#[test]
fn test_pair_lj_energy_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj");

    let parsed_topo = read_topology_file(base.join("pair_lj.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("pair_lj.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = gromos::core::math::Vacuum;
    let mut storage = ForceStorage::new(positions.len());
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        &mut storage,
    );

    // Expected VdW energy from input.toml: E_Vdw = -4.9395e-01
    // More precise from energies.tre: -4.939549393e-01
    let expected_e_lj = -4.939549393e-01;
    let tol = 1e-8;
    assert_f64_rel(storage.e_lj, expected_e_lj, tol, "pair_lj E_Vdw");
}

// ─── Level 0: nacl_pair ─────────────────────────────────────────────────────

#[test]
fn test_nacl_pair_forces_step0() {
    init_logging();
    let base = ref_dir().join("nacl_pair");

    let parsed_topo = read_topology_file(base.join("nacl_pair.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("nacl_pair.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let n_atoms = positions.len();
    assert_eq!(n_atoms, 2);

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);

    // For nacl_pair, electrostatics are active. Use CRF with large cutoff (vacuum, no PBC).
    // In vacuum with no PBC, cutoff > distance means plain Coulomb.
    // CRF formula: E = q*(1/r + crf_2cut3i*r^2 - crf_cut3i)
    // For no reaction field (epsrf=1, rcrf=rcutl): crf=0, so crf_2cut3i=0, crf_cut3i=1/rc
    // But for vacuum (no RF), just use large cutoff with crf terms = 0
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = gromos::core::math::Vacuum;
    let mut storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        &mut storage,
    );

    // Load expected forces
    let expected_forces = parse_forces_trf(&base.join("expected/forces.trf"), 0);
    assert_eq!(expected_forces.len(), n_atoms);

    let tol = 1e-6;
    for i in 0..n_atoms {
        assert_vec3_approx(
            storage.forces[i],
            expected_forces[i],
            tol,
            &format!("nacl_pair force atom {}", i),
        );
    }
}

#[test]
fn test_nacl_pair_energy_step0() {
    init_logging();
    let base = ref_dir().join("nacl_pair");

    let parsed_topo = read_topology_file(base.join("nacl_pair.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("nacl_pair.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = gromos::core::math::Vacuum;
    let mut storage = ForceStorage::new(positions.len());
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        &mut storage,
    );

    // Expected from input.toml: E_Vdw = 1.4529e+00 (repulsive LJ at short range)
    // E_Non_bonded = -4.6167e+02 (dominated by Coulomb attraction)
    let expected_e_vdw = 1.4529e+00;
    let tol = 1e-3; // Relaxed since we only have 4 sigfigs from input.toml
    assert_f64_rel(storage.e_lj, expected_e_vdw, tol, "nacl_pair E_Vdw");
}

// ─── Level 0: pair_lj_mixed ─────────────────────────────────────────────────

#[test]
fn test_pair_lj_mixed_forces_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj_mixed");

    let parsed_topo = read_topology_file(base.join("pair_lj_mixed.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("pair_lj_mixed.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let n_atoms = positions.len();
    assert_eq!(n_atoms, 2);

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = gromos::core::math::Vacuum;
    let mut storage = ForceStorage::new(n_atoms);
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        &mut storage,
    );

    // Load expected forces
    let expected_forces = parse_forces_trf(&base.join("expected/forces.trf"), 0);
    assert_eq!(expected_forces.len(), n_atoms);

    let tol = 1e-6;
    for i in 0..n_atoms {
        assert_vec3_approx(
            storage.forces[i],
            expected_forces[i],
            tol,
            &format!("pair_lj_mixed force atom {}", i),
        );
    }
}

#[test]
fn test_pair_lj_mixed_energy_step0() {
    init_logging();
    let base = ref_dir().join("pair_lj_mixed");

    let parsed_topo = read_topology_file(base.join("pair_lj_mixed.topo")).unwrap();
    let topo = build_topology(parsed_topo);
    let conf = read_coordinate_file(base.join("pair_lj_mixed.conf"), 1, 1).unwrap();
    let positions = &conf.current().pos;

    let pairlist: Vec<(u32, u32)> = vec![(0, 1)];
    let iac: Vec<u32> = topo.iac.iter().map(|&i| i as u32).collect();
    let charges: Vec<f64> = topo.charge.clone();
    let lj_params = convert_lj_params(&topo.lj_parameters);
    let crf = CRFParameters {
        crf_cut: 1.4,
        crf_2cut3i: 0.0,
        crf_cut3i: 0.0,
        cutoff_sq: 1.4_f64.powi(2),
    };

    let bc = gromos::core::math::Vacuum;
    let mut storage = ForceStorage::new(positions.len());
    lj_crf_innerloop(
        positions,
        &charges,
        &iac,
        &pairlist,
        &lj_params,
        &crf,
        &bc,
        &mut storage,
    );

    // Expected from input.toml: E_Vdw = -1.1219e+00
    let expected_e_vdw = -1.1219e+00;
    let tol = 1e-3;
    assert_f64_rel(storage.e_lj, expected_e_vdw, tol, "pair_lj_mixed E_Vdw");
}
