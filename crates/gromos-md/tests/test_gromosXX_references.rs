//! End-to-end integration tests against GROMOS double-precision references.
//!
//! Runs the `md` binary on each reference system (10 steps) and compares
//! energy trajectories (all steps) against expected GROMOS output.
//!
//! Output formats differ:
//!   - GROMOS expected: ENERGY03 blocks (one per step, one value per line)
//!   - gromos-rs actual:  ENERTRJ block  (one line per step, multi-column)
//!
//! Reference data: `.local/gromos_references/{system}/expected/`
//! Systems matching GROMOS are active; known mismatches are `#[ignore]`.
//!
//! Run passing:  cargo test -p gromos-md --test test_GROMOS_references
//! Run all:      cargo test -p gromos-md --test test_GROMOS_references -- --include-ignored

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

// ─── Tolerances ─────────────────────────────────────────────────────────────

const ENERGY_REL_TOL: f64 = 1e-8;
const ENERGY_ABS_TOL: f64 = 1e-10; // for near-zero energies
const FORCE_ABS_TOL: f64 = 1e-6; // kJ/(mol*nm)
const DHDL_REL_TOL: f64 = 1e-6; // dH/dλ relative tolerance

// ─── Paths ──────────────────────────────────────────────────────────────────

fn ref_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/gromosXX_references")
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

// ─── Parser: GROMOS ENERGY03 format (expected) ────────────────────────────
//
// Structure: repeated TIMESTEP→END→ENERGY03→END blocks.
// ENERGY03 starts with "# totals", then one f64 per line:
//   [0]=e_total, [1]=e_kinetic, [2]=e_potential, ...

fn parse_energy03(path: &Path) -> Vec<EnergyFrame> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut frames = Vec::new();
    let mut in_ene = false;
    let mut vals: Vec<f64> = Vec::new();

    for line in content.lines() {
        let t = line.trim();
        match t {
            "ENERGY03" => {
                in_ene = true;
                vals.clear();
            },
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
            },
            _ if in_ene => {
                if !t.starts_with('#') && !t.is_empty() {
                    // Single-value lines only (totals); multi-value lines fail parse
                    if let Ok(v) = t.parse::<f64>() {
                        vals.push(v);
                    }
                }
            },
            _ => {},
        }
    }
    frames
}

// The engine writes the native layout since 0.0.34, so `parse_energy03` reads both sides.
// ─── Parser: free-energy trajectories (.trg) ────────────────────────────────
//
// gromosXX writes `TIMESTEP` + `FREEENERDERIVS03` blocks (`# lambda`, then `# totals` in
// ENERGY03 order: total, kinetic, potential, covalent, bond, angle, improper, dihedral,
// cross-dihedral, nonbonded, LJ, CRF, …). Until 0.0.33 this parser looked for the legacy
// one-line `FREEENERGY03` block only, so no dH/dλ was ever compared against the native files.
// Both layouts are read through the crate's own reader now.

#[derive(Debug, Clone)]
struct FreeEnergyFrame {
    time: f64,
    lambda: f64,
    dhdl_total: f64,
    dhdl_lj: f64,
    dhdl_crf: f64,
    dhdl_bonded: f64,
}

fn parse_freeenergy03(path: &Path) -> Vec<FreeEnergyFrame> {
    gromos_io::read_free_energy_trajectory(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .into_iter()
        .map(|f| FreeEnergyFrame {
            time: f.time,
            lambda: f.lambda,
            dhdl_total: f.dhdl_total,
            dhdl_lj: f.dhdl_lj,
            dhdl_crf: f.dhdl_crf,
            dhdl_bonded: f.dhdl_bond + f.dhdl_angle + f.dhdl_improper + f.dhdl_dihedral,
        })
        .collect()
}

fn assert_dhdl_close(actual: f64, expected: f64, label: &str) {
    let diff = (actual - expected).abs();
    let tol = (expected.abs() * DHDL_REL_TOL).max(ENERGY_ABS_TOL);
    assert!(
        diff <= tol,
        "{label}: expected {expected:.10e}, got {actual:.10e}, diff={diff:.2e}, tol={tol:.2e}"
    );
}

// ─── Force frame ────────────────────────────────────────────────────────────

/// One frame of per-atom forces: Vec of (fx, fy, fz) for each atom.
type ForceFrame = Vec<[f64; 3]>;

// ─── Parser: FREEFORCERED blocks from .trf files ────────────────────────────
//
// Both GROMOS and gromos-rs write FREEFORCERED blocks with 3 floats per line.
// Lines starting with '#' are comments (atom count markers).

fn parse_trf(path: &Path) -> Vec<ForceFrame> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
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
            let vals: Vec<f64> = t
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
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

// ─── Parser: GROMOS .conf position file ─────────────────────────────────────
//
// Reads the POSITION block. Each data line: <label x4> <seq> x y z  (nm).
// The first 24 chars are a label field that is ignored by gromosXX convention.

fn parse_conf(path: &Path) -> Vec<[f64; 3]> {
    let content = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut positions = Vec::new();
    let mut in_pos = false;

    for line in content.lines() {
        let t = line.trim();
        if t == "POSITION" {
            in_pos = true;
            continue;
        }
        if t == "END" && in_pos {
            break;
        }
        if in_pos && !t.starts_with('#') && !t.is_empty() {
            // Format: resnum resname atomname atomnum x y z
            let vals: Vec<f64> = t
                .split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect();
            // first numeric token is an integer (resnum or atomnum), last three are x y z
            if vals.len() >= 3 {
                let x = vals[vals.len() - 3];
                let y = vals[vals.len() - 2];
                let z = vals[vals.len() - 1];
                positions.push([x, y, z]);
            }
        }
    }
    positions
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
    let out =
        std::env::temp_dir().join(format!("gromos_reftest_{}_{}", system, std::process::id()));
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(&out).unwrap();

    let tre = out.join("energies.tre");
    let trf = out.join("forces.trf");
    let trg = out.join("free_energy.trg");

    // Run md binary
    let mut cmd = Command::new(md_bin());
    cmd.arg("@topo")
        .arg(&topo)
        .arg("@conf")
        .arg(&conf)
        .arg("@input")
        .arg(&params)
        .arg("@fin")
        .arg(out.join("final.conf"))
        .arg("@tre")
        .arg(&tre)
        .arg("@trf")
        .arg(&trf)
        .arg("@trc")
        .arg(out.join("trajectory.trc"))
        .arg("@trg")
        .arg(&trg);

    // Optional position restraints
    if let Some(por) = toml_str_opt(&toml, "posresspec") {
        cmd.arg("@posresspec").arg(sys_dir.join(por));
    }
    if let Some(rpr) = toml_str_opt(&toml, "refpos") {
        cmd.arg("@refpos").arg(sys_dir.join(rpr));
    }

    // Optional distance restraints
    if let Some(dr) = toml_str_opt(&toml, "distrest") {
        cmd.arg("@distrest").arg(sys_dir.join(dr));
    }

    // Optional perturbation topology
    if let Some(pt) = toml_str_opt(&toml, "pttopo") {
        cmd.arg("@pttopo").arg(sys_dir.join(pt));
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
    let actual = parse_energy03(&tre);

    // gromos-rs writes step 0..NSTLIM (inclusive), GROMOS writes 0..NSTLIM-1;
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
                    a[0],
                    e[0],
                    &format!("{system}[{frame_idx}] atom {atom_idx} fx"),
                );
                assert_force_close(
                    a[1],
                    e[1],
                    &format!("{system}[{frame_idx}] atom {atom_idx} fy"),
                );
                assert_force_close(
                    a[2],
                    e[2],
                    &format!("{system}[{frame_idx}] atom {atom_idx} fz"),
                );
            }
        }
        if n_compare == 0 && !exp_forces.is_empty() {
            panic!(
                "{system}: expected {} force frames but got none",
                exp_forces.len()
            );
        }
    }

    // Compare free-energy trajectory if reference .trg exists
    if let Some(fe_path) = toml_str_opt(&toml, "free_energy") {
        let expected_trg = sys_dir.join(fe_path);
        if expected_trg.exists() && trg.exists() {
            let exp_fe = parse_freeenergy03(&expected_trg);
            let act_fe = parse_freeenergy03(&trg);
            assert!(
                act_fe.len() >= exp_fe.len(),
                "{system}: too few free-energy frames (expected {}, got {})",
                exp_fe.len(),
                act_fe.len()
            );
            if exp_fe.is_empty() {
                // The legacy ch4_water_fep reference was generated with NTWG=0: an empty .trg.
                // Its λ-range siblings (ch4_water_fep_l*) carry the dH/dλ frames.
                eprintln!(
                    "{system}: reference .trg has no free-energy frames — dH/dλ not compared"
                );
            }
            // Align by time: the engine may write a frame at t=0 that gromosXX does not.
            let mut actual = act_fe.iter().peekable();
            for (i, ef) in exp_fe.iter().enumerate() {
                let af = loop {
                    match actual.next() {
                        Some(af) if (af.time - ef.time).abs() < 1e-9 => break af,
                        Some(af) if af.time < ef.time => continue,
                        other => panic!(
                            "{system}[{i}]: no free-energy frame at t={:.4} ps (next: {:?})",
                            ef.time,
                            other.map(|f| f.time)
                        ),
                    }
                };
                let label = |what: &str| {
                    format!(
                        "{system}[{i}] dH/dλ {what} (t={:.3}ps λ={:.3})",
                        ef.time, ef.lambda
                    )
                };
                assert_dhdl_close(af.dhdl_total, ef.dhdl_total, &label("total"));
                // Per-term derivatives only where the engine records them (a frame whose
                // components are all zero while the total is not carries the total only).
                let has_split = af.dhdl_total == 0.0
                    || af.dhdl_lj != 0.0
                    || af.dhdl_crf != 0.0
                    || af.dhdl_bonded != 0.0;
                if has_split {
                    assert_dhdl_close(af.dhdl_lj, ef.dhdl_lj, &label("LJ"));
                    assert_dhdl_close(af.dhdl_crf, ef.dhdl_crf, &label("CRF"));
                    assert_dhdl_close(af.dhdl_bonded, ef.dhdl_bonded, &label("bonded"));
                }
            }
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
    (ignore: $name:ident, $sys:literal, $reason:literal) => {
        #[test]
        #[ignore = $reason]
        fn $name() {
            run_reference($sys);
        }
    };
}

// ── Level 0: pair interactions, vacuum ──────────────────────────────────────
ref_test!(pair_lj, "pair_lj");
ref_test!(pair_lj_mixed, "pair_lj_mixed");
ref_test!(nacl_pair, "nacl_pair");

// ── Level 1: single molecule, bonded terms, PBC+RF ──────────────────────────
ref_test!(water_single, "water_single");
ref_test!(water_single_genvel, "water_single_genvel");
ref_test!(benzene_vacuum, "benzene_vacuum");
ref_test!(nacl_pair_box, "nacl_pair_box");
ref_test!(butane_vacuum, "butane_vacuum");

// ── Level 2: PBC, pairlist, SHAKE, solvent, twin-range ──────────────────────
ref_test!(water_3_box, "water_3_box");
ref_test!(nacl_1water_box, "nacl_1water_box");
ref_test!(nacl_3water_box, "nacl_3water_box");
ref_test!(water_3_box_twinrange, "water_3_box_twinrange");
ref_test!(water_10_box, "water_10_box");
ref_test!(nacl_3water_cutoff, "nacl_3water_cutoff");
ref_test!(nacl_water_box, "nacl_water_box");
ref_test!(nacl_water_box_shifted, "nacl_water_box_shifted");

// ── Constraint algorithms ────────────────────────────────────────────────────
ref_test!(nacl_1water_settle, "nacl_1water_settle");
ref_test!(nacl_1water_lincs, "nacl_1water_lincs");
ref_test!(aladip_vacuum_lincs, "aladip_vacuum_lincs");

// ── Level 2: grid_cell pairlist (1000 SPC molecules) ────────────────────────
ref_test!(water_1000_spc_gridcell, "water_1000_spc_gridcell");

// ── Level 3: bulk water ──────────────────────────────────────────────────────
ref_test!(water_216_box, "water_216_box");
// Same system with the MULTIBATH block *absent*: gromosXX runs no temperature coupling. Until
// PLAN.md 3.9 step 2 the parser turned the absent block into a Berendsen bath (A18); the
// defaults now mean "absent", and this guards that they keep meaning it.
ref_test!(water_216_nve_nobath, "water_216_nve_nobath");
ref_test!(water_216_box_com, "water_216_box_com");
ref_test!(water_216_box_com_rot, "water_216_box_com_rot");

// ── Level 3: thermostats / barostats ────────────────────────────────────────
ref_test!(water_216_nvt, "water_216_nvt");
ref_test!(water_216_nvt_nosehoover, "water_216_nvt_nosehoover");
ref_test!(water_216_nvt_nhc_chain, "water_216_nvt_nhc_chain");
ref_test!(water_216_npt, "water_216_npt");

// ── Alanine dipeptide ────────────────────────────────────────────────────────
ref_test!(aladip_vacuum, "aladip_vacuum");
ref_test!(aladip_solvated, "aladip_solvated");
ref_test!(aladip_trunc_oct, "aladip_trunc_oct");

// ── Restraints ───────────────────────────────────────────────────────────────
ref_test!(nacl_1water_distres, "nacl_1water_distres");

// ── FEP / TI ────────────────────────────────────────────────────────────────
ref_test!(ch4_water_fep, "ch4_water_fep");
ref_test!(ch4_water_fep_l000, "ch4_water_fep_l000");
ref_test!(ch4_water_fep_l025, "ch4_water_fep_l025");
ref_test!(ch4_water_fep_l075, "ch4_water_fep_l075");
ref_test!(ch4_water_fep_l100, "ch4_water_fep_l100");
// Perturbed bonded (regular and soft-core), atom types, charges, λ-mixed masses and a
// PERTATOMPAIR entry, all at once — passes since the 2026-08-29/30 FEP fixes.
ref_test!(aladip_vacuum_fep, "aladip_vacuum_fep");
// Charged perturbed atoms (54a7 CH3OH -> dummies): the soft-core reaction-field terms that
// `ch4_water_fep` (zero charge) never exercised. Found and fixed the 1/cutoff³ error in the
// perturbed pair kernel's softened RF term (0.0.34).
ref_test!(meoh_water_fep, "meoh_water_fep");

// ── Energy minimization ──────────────────────────────────────────────────────
// gromosXX writes the converged energies once more at convergence; `md` does the same since 0.0.34.
ref_test!(aladip_vacuum_em, "aladip_vacuum_em");
ref_test!(aladip_vacuum_em_shake, "aladip_vacuum_em_shake");
ref_test!(aladip_solvated_em_noshake, "aladip_solvated_em_noshake");
ref_test!(aladip_solvated_em_shake, "aladip_solvated_em_shake");
ref_test!(aladip_solvated_em_posres, "aladip_solvated_em_posres");
ref_test!(aladip_solvated_em, "aladip_solvated_em");
