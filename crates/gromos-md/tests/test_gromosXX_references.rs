//! End-to-end integration tests against gromosXX double-precision references.
//!
//! Runs the `md` binary on each reference system (10 steps) and compares
//! energy trajectories (all steps) against expected gromosXX output.
//!
//! Output formats differ:
//!   - gromosXX expected: ENERGY03 blocks (one per step, one value per line)
//!   - gromos-rs actual:  ENERTRJ block  (one line per step, multi-column)
//!
//! Reference data: `.local/gromos_references/{system}/expected/`
//! Systems matching gromosXX are active; known mismatches are `#[ignore]`.
//!
//! Run passing:  cargo test -p gromos-md --test test_gromosXX_references
//! Run all:      cargo test -p gromos-md --test test_gromosXX_references -- --include-ignored

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─── Tolerances ─────────────────────────────────────────────────────────────

const ENERGY_REL_TOL: f64 = 1e-8;
const ENERGY_ABS_TOL: f64 = 1e-10; // for near-zero energies
const FORCE_ABS_TOL: f64 = 1e-6;   // kJ/(mol*nm)

// ─── Paths ──────────────────────────────────────────────────────────────────

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/gromosXX_references")
}

fn md_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_md"))
}

// ─── Minimal TOML field extraction ──────────────────────────────────────────

fn toml_str(content: &str, key: &str) -> String {
    let prefix = format!("{key} = ");
    for line in content.lines() {
        let t = line.trim();
        if t.starts_with(&prefix) {
            if let (Some(a), Some(b)) = (t.find('"'), t.rfind('"')) {
                if a < b {
                    return t[a + 1..b].to_string();
                }
            }
        }
    }
    panic!("{key} not found in input.toml");
}

fn toml_str_opt(content: &str, key: &str) -> Option<String> {
    let prefix = format!("{key} = ");
    for line in content.lines() {
        let t = line.trim();
        if t.starts_with(&prefix) {
            if let (Some(a), Some(b)) = (t.find('"'), t.rfind('"')) {
                if a < b {
                    return Some(t[a + 1..b].to_string());
                }
            }
        }
    }
    None
}

// ─── Energy frame ───────────────────────────────────────────────────────────

#[derive(Debug)]
struct EnergyFrame {
    e_total: f64,
    e_kinetic: f64,
    e_potential: f64,
}

// ─── Parser: gromosXX ENERGY03 format (expected) ────────────────────────────
//
// Structure: repeated TIMESTEP→END→ENERGY03→END blocks.
// ENERGY03 starts with "# totals", then one f64 per line:
//   [0]=e_total, [1]=e_kinetic, [2]=e_potential, ...

fn parse_energy03(path: &Path) -> Vec<EnergyFrame> {
    let content =
        fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut frames = Vec::new();
    let mut in_ene = false;
    let mut vals: Vec<f64> = Vec::new();

    for line in content.lines() {
        let t = line.trim();
        match t {
            "ENERGY03" => {
                in_ene = true;
                vals.clear();
            }
            "END" => {
                if in_ene {
                    if vals.len() >= 3 {
                        frames.push(EnergyFrame {
                            e_total: vals[0],
                            e_kinetic: vals[1],
                            e_potential: vals[2],
                        });
                    }
                    in_ene = false;
                }
            }
            _ if in_ene => {
                if !t.starts_with('#') && !t.is_empty() {
                    // Single-value lines only (totals); multi-value lines fail parse
                    if let Ok(v) = t.parse::<f64>() {
                        vals.push(v);
                    }
                }
            }
            _ => {}
        }
    }
    frames
}

// ─── Parser: gromos-rs ENERTRJ format (actual) ─────────────────────────────
//
// Structure: single ENERTRJ block, one line per step with columns:
//   [0]=time, [1]=E_kin, [2]=E_pot, [3]=E_tot, [4]=T, ...

fn parse_enertrj(path: &Path) -> Vec<EnergyFrame> {
    let content =
        fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut frames = Vec::new();
    let mut in_block = false;

    for line in content.lines() {
        let t = line.trim();
        if t == "ENERTRJ" {
            in_block = true;
            continue;
        }
        if t == "END" && in_block {
            break;
        }
        if in_block && !t.starts_with('#') && !t.is_empty() {
            let vals: Vec<f64> = t.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            // [0]=time [1]=Ekin [2]=Epot [3]=Etot [4]=T [5]=V [6]=P ...
            if vals.len() >= 4 {
                frames.push(EnergyFrame {
                    e_kinetic: vals[1],
                    e_potential: vals[2],
                    e_total: vals[3],
                });
            }
        }
    }
    frames
}

// ─── Force frame ────────────────────────────────────────────────────────────

/// One frame of per-atom forces: Vec of (fx, fy, fz) for each atom.
type ForceFrame = Vec<[f64; 3]>;

// ─── Parser: FREEFORCERED blocks from .trf files ────────────────────────────
//
// Both gromosXX and gromos-rs write FREEFORCERED blocks with 3 floats per line.
// Lines starting with '#' are comments (atom count markers).

fn parse_trf(path: &Path) -> Vec<ForceFrame> {
    let content =
        fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut frames = Vec::new();
    let mut in_block = false;
    let mut current_frame: ForceFrame = Vec::new();

    for line in content.lines() {
        let t = line.trim();
        if t == "FREEFORCERED" {
            in_block = true;
            current_frame.clear();
            continue;
        }
        if t == "END" && in_block {
            frames.push(current_frame.clone());
            in_block = false;
            continue;
        }
        if in_block && !t.starts_with('#') && !t.is_empty() {
            let vals: Vec<f64> = t.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            if vals.len() == 3 {
                current_frame.push([vals[0], vals[1], vals[2]]);
            }
        }
    }
    frames
}

fn assert_force_close(actual: f64, expected: f64, label: &str) {
    let diff = (actual - expected).abs();
    assert!(
        diff <= FORCE_ABS_TOL,
        "{label}: expected {expected:.10e}, got {actual:.10e}, diff={diff:.2e}, tol={FORCE_ABS_TOL:.2e}"
    );
}

// ─── Comparison ─────────────────────────────────────────────────────────────

fn assert_energy_close(actual: f64, expected: f64, label: &str) {
    let diff = (actual - expected).abs();
    let tol = (expected.abs() * ENERGY_REL_TOL).max(ENERGY_ABS_TOL);
    assert!(
        diff <= tol,
        "{label}: expected {expected:.10e}, got {actual:.10e}, diff={diff:.2e}, tol={tol:.2e}"
    );
}

// ─── Test driver ────────────────────────────────────────────────────────────

fn run_reference(system: &str) {
    let sys_dir = ref_root().join(system);
    if !sys_dir.exists() {
        eprintln!(
            "SKIP {system}: reference data not found at {}",
            sys_dir.display()
        );
        return;
    }

    // Read input file paths from input.toml
    let toml = fs::read_to_string(sys_dir.join("input.toml")).expect("input.toml missing");
    let topo = sys_dir.join(toml_str(&toml, "topology"));
    let conf = sys_dir.join(toml_str(&toml, "configuration"));
    let params = sys_dir.join(toml_str(&toml, "parameters"));

    // Output directory
    let out = std::env::temp_dir().join(format!(
        "gromos_reftest_{}_{}",
        system,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(&out).unwrap();

    let tre = out.join("energies.tre");
    let trf = out.join("forces.trf");

    // Run md binary
    let mut cmd = Command::new(md_bin());
    cmd.arg("@topo").arg(&topo)
        .arg("@conf").arg(&conf)
        .arg("@input").arg(&params)
        .arg("@fin").arg(out.join("final.conf"))
        .arg("@tre").arg(&tre)
        .arg("@trf").arg(&trf)
        .arg("@trc").arg(out.join("trajectory.trc"));

    // Optional position restraints
    if let Some(por) = toml_str_opt(&toml, "posresspec") {
        cmd.arg("@posresspec").arg(sys_dir.join(por));
    }
    if let Some(rpr) = toml_str_opt(&toml, "refpos") {
        cmd.arg("@refpos").arg(sys_dir.join(rpr));
    }

    let t0 = std::time::Instant::now();
    let result = cmd.output().expect("failed to execute md");
    let elapsed = t0.elapsed();

    eprintln!("  {system}: {:.3}ms", elapsed.as_secs_f64() * 1000.0);

    assert!(
        result.status.success(),
        "{system}: md exited with {}\nstdout: {}\nstderr: {}",
        result.status,
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );

    // Parse energies: expected (ENERGY03) vs actual (ENERTRJ)
    let expected = parse_energy03(&sys_dir.join("expected/energies.tre"));
    let actual = parse_enertrj(&tre);

    // gromos-rs writes step 0..NSTLIM (inclusive), gromosXX writes 0..NSTLIM-1;
    // compare the frames that exist in both outputs.
    assert!(
        actual.len() >= expected.len(),
        "{system}: too few frames (expected {}, got {})",
        expected.len(),
        actual.len()
    );

    for (i, (exp, act)) in expected.iter().zip(&actual).enumerate() {
        assert_energy_close(act.e_total, exp.e_total, &format!("{system}[{i}] E_total"));
        assert_energy_close(
            act.e_kinetic,
            exp.e_kinetic,
            &format!("{system}[{i}] E_kinetic"),
        );
        assert_energy_close(
            act.e_potential,
            exp.e_potential,
            &format!("{system}[{i}] E_potential"),
        );
    }

    // Compare forces if reference .trf exists
    let expected_trf = sys_dir.join("expected/forces.trf");
    if expected_trf.exists() && trf.exists() {
        let exp_forces = parse_trf(&expected_trf);
        let act_forces = parse_trf(&trf);

        let n_compare = exp_forces.len().min(act_forces.len());
        for (frame_idx, (ef, af)) in exp_forces.iter().zip(&act_forces).enumerate() {
            assert_eq!(
                ef.len(),
                af.len(),
                "{system}[frame {frame_idx}]: atom count mismatch (expected {}, got {})",
                ef.len(),
                af.len()
            );
            for (atom_idx, (e, a)) in ef.iter().zip(af.iter()).enumerate() {
                assert_force_close(
                    a[0], e[0],
                    &format!("{system}[{frame_idx}] atom {atom_idx} fx"),
                );
                assert_force_close(
                    a[1], e[1],
                    &format!("{system}[{frame_idx}] atom {atom_idx} fy"),
                );
                assert_force_close(
                    a[2], e[2],
                    &format!("{system}[{frame_idx}] atom {atom_idx} fz"),
                );
            }
        }
        if n_compare == 0 && !exp_forces.is_empty() {
            panic!("{system}: expected {} force frames but got none", exp_forces.len());
        }
    }

    let _ = fs::remove_dir_all(&out);
}

// ─── Test declarations ──────────────────────────────────────────────────────

macro_rules! ref_test {
    ($name:ident, $sys:literal) => {
        #[test]
        fn $name() {
            run_reference($sys);
        }
    };
    (ignore: $name:ident, $sys:literal) => {
        #[test]
        #[ignore = "known mismatch — see PLAN.md"]
        fn $name() {
            run_reference($sys);
        }
    };
}

// Level 0 — pair interactions, vacuum
ref_test!(pair_lj,       "pair_lj");
ref_test!(pair_lj_mixed, "pair_lj_mixed");
ref_test!(nacl_pair,     "nacl_pair");

// Level 1 — single molecule, bonded terms, PBC+RF
ref_test!(water_single,   "water_single");
ref_test!(water_single_genvel, "water_single_genvel");
ref_test!(benzene_vacuum, "benzene_vacuum");
ref_test!(nacl_pair_box,  "nacl_pair_box");

// Level 1 — dihedral + 1-4 LJ interaction
ref_test!(butane_vacuum,  "butane_vacuum");

// Level 2 — PBC, pairlist, SHAKE, solvent, twin-range
ref_test!(water_3_box,            "water_3_box");
ref_test!(nacl_1water_box,        "nacl_1water_box");
ref_test!(nacl_3water_box,        "nacl_3water_box");
ref_test!(water_3_box_twinrange,  "water_3_box_twinrange");
ref_test!(water_10_box,           "water_10_box");
ref_test!(nacl_3water_cutoff,     "nacl_3water_cutoff");
ref_test!(nacl_water_box,         "nacl_water_box");
ref_test!(nacl_water_box_shifted, "nacl_water_box_shifted");

// Level 3 — bulk water, pairlist scaling, long-range
ref_test!(water_216_box,       "water_216_box");
ref_test!(water_216_box_com,   "water_216_box_com");

// Level 3 — thermostats / barostats
ref_test!(water_216_nvt,       "water_216_nvt");
ref_test!(water_216_npt,       "water_216_npt");

ref_test!(aladip_vacuum,       "aladip_vacuum");
ref_test!(aladip_solvated,     "aladip_solvated");

// Energy minimization
ref_test!(aladip_vacuum_em,            "aladip_vacuum_em");
ref_test!(aladip_vacuum_em_shake,      "aladip_vacuum_em_shake");
ref_test!(aladip_solvated_em_noshake,  "aladip_solvated_em_noshake");
ref_test!(aladip_solvated_em_shake,    "aladip_solvated_em_shake");
ref_test!(aladip_solvated_em_posres,   "aladip_solvated_em_posres");
ref_test!(aladip_solvated_em,          "aladip_solvated_em");
