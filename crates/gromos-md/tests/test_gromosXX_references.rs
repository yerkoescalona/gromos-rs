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
/// The final configuration after the whole run: the frame-by-frame agreement is 1e-8 relative in
/// the energies, and the positions accumulate that over every step, so this is looser than the
/// per-step position tolerance.
const FINAL_POSITION_ABS_TOL: f64 = 1e-6; // nm
/// Trajectory frames are written with nine decimals, so agreement is exact up to that rounding.
const TRAJECTORY_POSITION_ABS_TOL: f64 = 1e-8; // nm
const FINAL_VELOCITY_ABS_TOL: f64 = 1e-4; // nm/ps

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
            _ if in_ene && !t.starts_with('#') && !t.is_empty() => {
                // Single-value lines only (totals); multi-value lines fail parse
                if let Ok(v) = t.parse::<f64>() {
                    vals.push(v);
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
    let traj = gromos_io::read_free_energy_trajectory(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    for w in &traj.warnings {
        eprintln!("{}: {w}", path.display());
    }
    traj.frames
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
        .arg("@trv")
        .arg(out.join("velocities.trv"))
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

    // Optional GaMD specification
    if let Some(g) = toml_str_opt(&toml, "gamd") {
        cmd.arg("@gamd").arg(sys_dir.join(g));
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

    // Both engines run exactly NSTLIM steps and write one frame per step, so the counts must
    // match. (They did not until 0.0.42: our loop ran NSTLIM+1 steps, and this assertion used to
    // allow the extra frame — which is precisely why the step count went unnoticed.)
    assert_eq!(
        actual.len(),
        expected.len(),
        "{system}: frame count differs from gromosXX (expected {}, got {})",
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

    // ── Trajectory ───────────────────────────────────────────────────────────
    // A frame's structure must belong with the energies of the same frame: gromosXX writes both at
    // the same point of the step, so its frame 0 is the input configuration. Ours wrote the state
    // *after* the step until 0.0.43, which left the two outputs one step out of step.
    let expected_trc = sys_dir.join("expected/trajectory.trc");
    if expected_trc.exists() {
        let (ours, _) = parse_configuration(&out.join("trajectory.trc"));
        let (theirs, _) = parse_configuration(&expected_trc);

        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: trajectory position count differs from gromosXX"
        );
        for (i, (x, y)) in ours.iter().zip(&theirs).enumerate() {
            let d = fold_into_box(&expected_trc, [x[0] - y[0], x[1] - y[1], x[2] - y[2]]);
            for (c, dc) in d.iter().enumerate() {
                assert!(
                    dc.abs() <= TRAJECTORY_POSITION_ABS_TOL,
                    "{system}: trajectory position {} component {c}: {:.10e} vs gromosXX {:.10e}",
                    i + 1,
                    x[c],
                    y[c]
                );
            }
        }
    }

    // ── Velocity trajectory (@trv / NTWV) ────────────────────────────────────
    let expected_trv = sys_dir.join("expected/velocities.trv");
    if expected_trv.exists() {
        let (_, ours) = parse_configuration(&out.join("velocities.trv"));
        let (_, theirs) = parse_configuration(&expected_trv);
        // `@trv` was accepted and ignored until 0.0.44; the count is what proves it is written.
        // The *values* are not compared yet: velocities are half-step quantities in leap-frog and
        // gromosXX's frame 0 is the input velocities, while ours is one half-step further on — the
        // residual from frame 1 onwards is ~1e-4 nm/ps, the size of one thermostat scaling. Which
        // half-step each frame carries is an open item (PLAN.md 1.5b).
        assert_eq!(
            ours.len(),
            theirs.len(),
            "{system}: velocity trajectory count differs from gromosXX (is @trv written at all?)"
        );
        assert!(
            !ours.is_empty(),
            "{system}: @trv produced no velocity frames"
        );
    }

    // ── Final configuration ──────────────────────────────────────────────────
    // The end state is what a continuation run starts from, and it is the only place a wrong
    // step count shows up — the per-step frames look right either way.
    let expected_fin = sys_dir.join("expected/final.conf");
    if expected_fin.exists() {
        let ours = parse_configuration(&out.join("final.conf"));
        let theirs = parse_configuration(&expected_fin);
        // Which periodic image a configuration file records is not physics: gromosXX puts charge
        // groups back into the box at the start of each step, so an atom that the last position
        // update pushed just outside is written on one side by one engine and on the other by the
        // other. Differences are compared modulo a box vector; anything else is a real deviation.

        for (what, a, b) in [
            ("position", &ours.0, &theirs.0),
            ("velocity", &ours.1, &theirs.1),
        ] {
            if b.is_empty() {
                continue;
            }
            assert_eq!(
                a.len(),
                b.len(),
                "{system}: final {what} count differs from gromosXX"
            );
            for (i, (x, y)) in a.iter().zip(b).enumerate() {
                for (c, (xc, yc)) in x.iter().zip(y).enumerate() {
                    let tol = if what == "position" {
                        FINAL_POSITION_ABS_TOL
                    } else {
                        FINAL_VELOCITY_ABS_TOL
                    };
                    let mut diff = xc - yc;
                    if what == "position" {
                        diff =
                            fold_into_box(&expected_fin, [x[0] - y[0], x[1] - y[1], x[2] - y[2]])
                                [c];
                    }
                    assert!(
                        diff.abs() <= tol,
                        "{system}: final {what} of atom {} component {c}: {xc:.10e} vs gromosXX {yc:.10e} (tol {tol:.1e})",
                        i + 1
                    );
                }
            }
        }
    }

    let _ = fs::remove_dir_all(&out);
}

/// Fold a coordinate difference into the periodic images of a reference `.cnf`'s box, so that a
/// comparison does not depend on which image an engine happened to write. Configurations are
/// written in the frame of the input file, so NTB = −1 is folded with the truncated-octahedron
/// rule *in the cube frame* (the cubic minimum image, then the body-centred (±L/2)³ shift), not
/// with the triclinic lattice the engine works in.
fn fold_into_box(path: &Path, d: [f64; 3]) -> [f64; 3] {
    let Some((ntb, edges)) = parse_box(path) else {
        return d;
    };
    if edges[0] <= 0.0 || ntb == 0 {
        return d;
    }
    let mut r = [0.0; 3];
    for c in 0..3 {
        let l = edges[c];
        r[c] = if l > 0.0 {
            d[c] - l * (d[c] / l).round()
        } else {
            d[c]
        };
    }
    if ntb == -1 {
        let l = edges[0];
        if r[0].abs() + r[1].abs() + r[2].abs() > 0.75 * l {
            for (c, rc) in r.iter_mut().enumerate() {
                let _ = c;
                *rc -= 0.5 * l * rc.signum();
            }
        }
    }
    r
}

/// NTB and the box edge lengths of a `.cnf`'s GENBOX block, if it has one.
fn parse_box(path: &Path) -> Option<(i32, [f64; 3])> {
    let text = fs::read_to_string(path).ok()?;
    let mut lines = text
        .lines()
        .skip_while(|l| l.trim() != "GENBOX")
        .skip(1)
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'));
    let ntb: i32 = lines.next()?.split_whitespace().next()?.parse().ok()?;
    let v: Vec<f64> = lines
        .next()?
        .split_whitespace()
        .filter_map(|x| x.parse().ok())
        .collect();
    (v.len() >= 3).then_some((ntb, [v[0], v[1], v[2]]))
}

/// POSITION and VELOCITY blocks of a `.cnf`, as `[x, y, z]` per atom.
fn parse_configuration(path: &Path) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let (mut pos, mut vel) = (Vec::new(), Vec::new());
    let mut block: Option<&str> = None;
    for line in text.lines() {
        let t = line.trim();
        if t == "POSITION" || t == "POSITIONRED" {
            block = Some("pos");
            continue;
        }
        if t == "VELOCITY" || t == "VELOCITYRED" {
            block = Some("vel");
            continue;
        }
        if t == "END" {
            block = None;
            continue;
        }
        let Some(kind) = block else { continue };
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let v: Vec<f64> = t
            .split_whitespace()
            .filter_map(|x| x.parse::<f64>().ok())
            .collect();
        if v.len() >= 3 {
            let xyz = [v[v.len() - 3], v[v.len() - 2], v[v.len() - 1]];
            if kind == "pos" {
                pos.push(xyz);
            } else {
                vel.push(xyz);
            }
        }
    }
    (pos, vel)
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

// ── Level 5: the input styles real GROMOS files come in ─────────────────────
// Added 2026-08-30 after the LiveCoMS tutorial suite exposed a class of defect the levels above
// cannot see: every system there uses one temperature bath, one value per line, and molecules that
// are never wrapped across the box. See README-livecoms.md in the reference directory.
ref_test!(aladip_multibath, "aladip_multibath");
ref_test!(aladip_multibath_collapsed, "aladip_multibath_collapsed");
ref_test!(aladip_wrapped, "aladip_wrapped");
ref_test!(aladip_harmonic_covalent, "aladip_harmonic_covalent");
// The classical part of this system already matches gromosXX exactly (its potential energy at
// frame 0 agrees to 1e-9); what differs is the GaMD boost itself. gromosXX accumulates the
// dihedral and total forces of every acceleration group during the force calculation and scales
// them by k·(V−E), recording ΔV = k(V−E)²/2 in the special energy — on this system ΔV ≈ 1.23
// kJ/mol, which our run neither applies the same way nor reports. The reference stays here as the
// specification for that work (PLAN.md 1.9).
ref_test!(
    ignore: aladip_gamd,
    "aladip_gamd",
    "GaMD boost unvalidated: gromosXX's per-acceleration-group force scaling and its ΔV in the special energy are not reproduced (PLAN.md 1.9)"
);
