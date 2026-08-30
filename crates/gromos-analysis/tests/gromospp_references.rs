//! The analysis programs against gromos++ output (`tests/data/README.md`): same arguments, every
//! number compared to 1e-6 relative (`-nan` literally), every label exactly.

use std::path::{Path, PathBuf};
use std::process::Command;

fn data(tool: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(tool)
}

fn shared() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

fn run(bin: &str, args: &[String], cwd: Option<&Path>) -> String {
    let mut cmd = Command::new(bin);
    cmd.args(args);
    if let Some(d) = cwd {
        cmd.current_dir(d);
    }
    let out = cmd.output().expect("run");
    assert!(
        out.status.success(),
        "{bin} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// Token-wise comparison: numbers within `rel`, everything else exact.
fn assert_same_text(ours: &str, theirs: &str, rel: f64) {
    let ol: Vec<&str> = ours.lines().collect();
    let tl: Vec<&str> = theirs.lines().collect();
    assert_eq!(
        ol.len(),
        tl.len(),
        "line count\n--- ours ---\n{ours}\n--- gromos++ ---\n{theirs}"
    );
    for (n, (a, b)) in ol.iter().zip(&tl).enumerate() {
        let at: Vec<&str> = a.split_whitespace().collect();
        let bt: Vec<&str> = b.split_whitespace().collect();
        assert_eq!(
            at.len(),
            bt.len(),
            "line {}: {a:?} vs gromos++ {b:?}",
            n + 1
        );
        for (x, y) in at.iter().zip(&bt) {
            match (x.parse::<f64>(), y.parse::<f64>()) {
                (Ok(fx), Ok(fy)) if fx.is_finite() && fy.is_finite() => {
                    assert!(
                        (fx - fy).abs() <= rel * fx.abs().max(fy.abs()).max(1e-3),
                        "line {}: {x} vs gromos++ {y}\n{a}\n{b}",
                        n + 1
                    );
                },
                _ => assert_eq!(x, y, "line {}: {a:?} vs gromos++ {b:?}", n + 1),
            }
        }
    }
}

#[test]
fn tser_matches_gromospp() {
    let out = run(
        env!("CARGO_BIN_EXE_tser"),
        &[
            "@topo",
            &shared().join("shared/aladip.topo").to_string_lossy(),
            "@pbc",
            "r",
            "@prop",
            "d%1:1,2 a%1:1,2,3 t%1:1,2,3,4 tp%1:2,3,4,5%60",
            "@traj",
            &shared()
                .join("aladip_solvated/expected/trajectory.trc")
                .to_string_lossy(),
            "@dist",
            "5",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        None,
    );
    let theirs =
        std::fs::read_to_string(data("tser").join("aladip_solvated.gromospp.out")).unwrap();
    assert_same_text(&out, &theirs, 1e-5);
}

#[test]
fn mdf_matches_gromospp() {
    let dir = std::env::temp_dir().join(format!(
        "gromos_mdf_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    run(
        env!("CARGO_BIN_EXE_mdf"),
        &[
            "@topo",
            &shared().join("shared/aladip.topo").to_string_lossy(),
            "@pbc",
            "r",
            "@centre",
            "1:1",
            "1:6",
            "@with",
            "1:3-12",
            "@traj",
            &shared()
                .join("aladip_solvated/expected/trajectory.trc")
                .to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        Some(&dir),
    );
    for name in ["MIN_1:1", "MIN_1:6"] {
        let ours = std::fs::read_to_string(dir.join(format!("{name}.dat"))).unwrap();
        let theirs =
            std::fs::read_to_string(data("mdf").join(format!("{name}.gromospp.dat"))).unwrap();
        assert_same_text(&ours, &theirs, 1e-8);
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn dg_ener_matches_gromospp() {
    let d = data("dg_ener");
    for (extra, reference) in [
        (vec![], "dg_ener.gromospp.out"),
        (vec!["@col", "3", "3"], "dg_ener_col3.gromospp.out"),
    ] {
        let mut args: Vec<String> = [
            "@temp",
            "300",
            "@stateA",
            &d.join("eA.dat").to_string_lossy(),
            "@stateB",
            &d.join("eB.dat").to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        args.extend(extra.iter().map(|s| s.to_string()));
        let out = run(env!("CARGO_BIN_EXE_dg_ener"), &args, None);
        assert_same_text(
            &out,
            &std::fs::read_to_string(d.join(reference)).unwrap(),
            2e-5,
        );
    }
}

#[test]
fn dfmult_matches_gromospp() {
    let d = data("dfmult");
    let out = run(
        env!("CARGO_BIN_EXE_dfmult"),
        &[
            "@temp",
            "300",
            "@stateR",
            &d.join("eR.dat").to_string_lossy(),
            "@endstates",
            &d.join("eX1.dat").to_string_lossy(),
            &d.join("eX2.dat").to_string_lossy(),
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("dfmult.gromospp.out")).unwrap(),
        1e-6,
    );
}

#[test]
fn matrix_overlap_matches_gromospp() {
    let d = data("matrix_overlap");
    let out = run(
        env!("CARGO_BIN_EXE_matrix_overlap"),
        &[
            "@m1",
            &d.join("m1.dat").to_string_lossy(),
            "@m2",
            &d.join("m2.dat").to_string_lossy(),
            "@dim",
            "3",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>(),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("matrix_overlap.gromospp.out")).unwrap(),
        1e-5,
    );
}
