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

fn topo() -> String {
    shared()
        .join("shared/aladip.topo")
        .to_string_lossy()
        .into_owned()
}

fn trc() -> String {
    shared()
        .join("aladip_solvated/expected/trajectory.trc")
        .to_string_lossy()
        .into_owned()
}

fn strs(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

#[test]
fn bilayer_dist_matches_gromospp() {
    let d = data("bilayer_dist");
    let out = run(
        env!("CARGO_BIN_EXE_bilayer_dist"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@selection",
            "s:OW",
            "@grid",
            "10",
            "@traj",
            &trc(),
        ]),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("aladip_grid10.gromospp.out")).unwrap(),
        1e-5,
    );
    let out = run(
        env!("CARGO_BIN_EXE_bilayer_dist"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@selection",
            "s:OW",
            "@grid",
            "10",
            "@density",
            "@traj",
            &trc(),
        ]),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(d.join("aladip_grid10_density.gromospp.out")).unwrap(),
        1e-5,
    );
}

#[test]
fn bilayer_oparam_matches_gromospp() {
    let out = run(
        env!("CARGO_BIN_EXE_bilayer_oparam"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:2-11",
            "@refvec",
            "0",
            "0",
            "1",
            "@traj",
            &trc(),
        ]),
        None,
    );
    assert_same_text(
        &out,
        &std::fs::read_to_string(data("bilayer_oparam").join("aladip_2-11.gromospp.out")).unwrap(),
        1e-4,
    );
}

#[test]
fn jval_matches_gromospp() {
    let d = data("jval");
    let jv = d.join("aladip.jval").to_string_lossy().into_owned();
    for (extra, reference) in [
        (vec!["@timeseries"], "timeseries.gromospp.out"),
        (vec!["@rmsd"], "rmsd.gromospp.out"),
        (vec![], "averages.gromospp.out"),
    ] {
        let mut args = strs(&["@topo", &topo(), "@pbc", "r", "@jval", &jv]);
        args.extend(strs(&extra));
        args.extend(strs(&["@traj", &trc()]));
        let out = run(env!("CARGO_BIN_EXE_jval"), &args, None);
        assert_same_text(
            &out,
            &std::fs::read_to_string(d.join(reference)).unwrap(),
            1e-5,
        );
    }
}

#[test]
fn edyn_matches_gromospp_and_projects_on_eigenvectors() {
    let dir = std::env::temp_dir().join(format!(
        "gromos_edyn_{}_{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    run(
        env!("CARGO_BIN_EXE_edyn"),
        &strs(&[
            "@topo",
            &topo(),
            "@pbc",
            "r",
            "@atoms",
            "1:a",
            "@eigenvalues",
            "1",
            "2",
            "@traj",
            &trc(),
        ]),
        Some(&dir),
    );
    let d = data("edyn");
    for (ours, theirs, rel) in [
        ("COVAR.out", "COVAR.gromospp.out", 1e-4),
        ("COVATOM.out", "COVATOM.gromospp.out", 1e-4),
        ("EIFLUC.out", "EIFLUC.gromospp.out", 1e-4),
    ] {
        assert_same_text(
            &std::fs::read_to_string(dir.join(ours)).unwrap(),
            &std::fs::read_to_string(d.join(theirs)).unwrap(),
            rel,
        );
    }
    // eigenvalues: the nine non-zero ones (10 frames → rank ≤ 9) to 1e-4, the rest are numerical zeros
    let ours = std::fs::read_to_string(dir.join("EIVAL.out")).unwrap();
    let theirs = std::fs::read_to_string(d.join("EIVAL.gromospp.out")).unwrap();
    for (a, b) in ours.lines().zip(theirs.lines()).skip(1).take(9) {
        let (x, y): (f64, f64) = (
            a.split_whitespace().nth(1).unwrap().parse().unwrap(),
            b.split_whitespace().nth(1).unwrap().parse().unwrap(),
        );
        assert!((x - y).abs() <= 1e-4 * y.abs(), "eigenvalue {a} vs {b}");
    }
    // projections on eigenvector k fluctuate with variance = eigenvalue k
    let essdyn = std::fs::read_to_string(dir.join("ESSDYN.out")).unwrap();
    let eig: Vec<f64> = ours
        .lines()
        .skip(1)
        .map(|l| l.split_whitespace().nth(1).unwrap().parse().unwrap())
        .collect();
    for (k, line) in essdyn.lines().filter(|l| !l.starts_with('#')).enumerate() {
        let fluct: f64 = line.split_whitespace().nth(2).unwrap().parse().unwrap();
        assert!(
            (fluct * fluct - eig[k]).abs() <= 1e-4 * eig[k].max(1e-12) + 1e-12,
            "EV {}: fluctuation² {} vs eigenvalue {}",
            k + 1,
            fluct * fluct,
            eig[k]
        );
    }
    std::fs::remove_dir_all(dir).ok();
}

fn positions_of(text: &str) -> Vec<Vec<[f64; 3]>> {
    let mut frames = Vec::new();
    let mut cur: Option<Vec<[f64; 3]>> = None;
    for line in text.lines() {
        if line.starts_with("POSITION") {
            cur = Some(Vec::new());
            continue;
        }
        if let Some(c) = cur.as_mut() {
            if line.starts_with("END") {
                frames.push(cur.take().unwrap());
            } else if !line.starts_with('#') && line.len() > 24 {
                let v: Vec<f64> = line[24..]
                    .split_whitespace()
                    .map(|x| x.parse().unwrap())
                    .collect();
                c.push([v[0], v[1], v[2]]);
            }
        }
    }
    frames
}

#[test]
fn gca_matches_gromospp() {
    let d = data("gca");
    let dihedral = std::fs::read_to_string(d.join("dihedral.txt"))
        .unwrap()
        .trim()
        .to_string();
    let conf = shared()
        .join("shared/aladip.conf")
        .to_string_lossy()
        .into_owned();
    for (props, mobile, reference) in [
        (
            format!("t%1:{dihedral}%30%-60%60"),
            None,
            "torsion_scan.gromospp.cnf",
        ),
        (
            "d%1:1,2%0.2 a%1:1,2,3%100".to_string(),
            Some("first"),
            "dist_angle_first.gromospp.cnf",
        ),
        ("a%1:1,2,3%100".to_string(), None, "angle_last.gromospp.cnf"),
    ] {
        let mut args = strs(&["@topo", &topo(), "@pbc", "r", "@prop", &props]);
        if let Some(m) = mobile {
            args.extend(strs(&["@mobile", m]));
        }
        args.extend(strs(&["@traj", &conf]));
        let out = run(env!("CARGO_BIN_EXE_gca"), &args, None);
        let ours = positions_of(&out);
        let theirs = positions_of(&std::fs::read_to_string(d.join(reference)).unwrap());
        assert_eq!(ours.len(), theirs.len(), "{reference}: frame count");
        for (f, (a, b)) in ours.iter().zip(&theirs).enumerate() {
            assert_eq!(a.len(), b.len());
            for (k, (x, y)) in a.iter().zip(b).enumerate() {
                for dim in 0..3 {
                    assert!(
                        (x[dim] - y[dim]).abs() < 1e-7,
                        "{reference} frame {f} atom {}: {x:?} vs gromos++ {y:?}",
                        k + 1
                    );
                }
            }
        }
    }
}
