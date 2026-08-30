//! `copy_box`, `inbox` and `ion` against gromos++ output for the alanine dipeptide in 20 SPC
//! waters (`tests/data/{copy_box,inbox,ion}/README.md`).

use std::path::{Path, PathBuf};
use std::process::Command;

fn data(tool: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(tool)
}

fn shared() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references/shared")
}

/// (residue name, atom name, position) per atom, plus the box lengths from GENBOX.
fn parse(text: &str) -> (Vec<(String, String, [f64; 3])>, [f64; 3]) {
    let mut atoms = Vec::new();
    let mut box_dims = [0.0; 3];
    let mut block = "";
    let mut n = 0;
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if block.is_empty() {
            block = t;
            n = 0;
            continue;
        }
        if t == "END" {
            block = "";
            continue;
        }
        n += 1;
        match block {
            "POSITION" => {
                let v: Vec<f64> = line[24..]
                    .split_whitespace()
                    .map(|x| x.parse().unwrap())
                    .collect();
                atoms.push((
                    line[5..11].trim().to_string(),
                    line[11..17].trim().to_string(),
                    [v[0], v[1], v[2]],
                ));
            },
            "GENBOX" if n == 2 => {
                let v: Vec<f64> = t.split_whitespace().map(|x| x.parse().unwrap()).collect();
                box_dims = [v[0], v[1], v[2]];
            },
            _ => {},
        }
    }
    (atoms, box_dims)
}

fn run(bin: &str, args: &[&str]) -> String {
    let out = Command::new(bin).args(args).output().expect("run tool");
    assert!(
        out.status.success(),
        "{bin} failed:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn assert_positions(ours: &[(String, String, [f64; 3])], theirs: &[(String, String, [f64; 3])]) {
    assert_eq!(ours.len(), theirs.len(), "atom count");
    for (k, (a, b)) in ours.iter().zip(theirs).enumerate() {
        assert_eq!((&a.0, &a.1), (&b.0, &b.1), "labels of atom {}", k + 1);
        for d in 0..3 {
            assert!(
                (a.2[d] - b.2[d]).abs() < 1e-8,
                "atom {}: {:?} vs gromos++ {:?}",
                k + 1,
                a.2,
                b.2
            );
        }
    }
}

#[test]
fn copy_box_along_x_matches_gromospp() {
    let topo = shared().join("aladip.topo");
    let pos = shared().join("aladip.conf");
    let out = run(
        env!("CARGO_BIN_EXE_copy_box"),
        &[
            "@topo",
            &topo.to_string_lossy(),
            "@pos",
            &pos.to_string_lossy(),
            "@dir",
            "x",
        ],
    );
    let (ours, obox) = parse(&out);
    let (theirs, tbox) =
        parse(&std::fs::read_to_string(data("copy_box").join("aladip_x.gromospp.cnf")).unwrap());
    assert_eq!(theirs.len(), 144);
    assert_positions(&ours, &theirs);
    // gromos++ run without @pbc (to avoid gathering) writes a zero GENBOX; ours keeps the
    // grown box, which is what the duplicated system needs.
    if tbox != [0.0; 3] {
        for d in 0..3 {
            assert!((obox[d] - tbox[d]).abs() < 1e-8, "box {obox:?} vs {tbox:?}");
        }
    }
    assert!(
        (obox[0] - 2.0 * 3.767055681).abs() < 1e-8 && (obox[1] - 3.767055681).abs() < 1e-8,
        "box {obox:?}"
    );
}

#[test]
fn inbox_matches_gromospp() {
    let topo = shared().join("aladip.topo");
    let pos = shared().join("aladip.conf");
    let out = run(
        env!("CARGO_BIN_EXE_inbox"),
        &[
            "@topo",
            &topo.to_string_lossy(),
            "@pos",
            &pos.to_string_lossy(),
            "@pbc",
            "r",
        ],
    );
    let (ours, _) = parse(&out);
    let (theirs, _) =
        parse(&std::fs::read_to_string(data("inbox").join("aladip.gromospp.cnf")).unwrap());
    assert_eq!(theirs.len(), 72);
    assert_positions(&ours, &theirs);
}

#[test]
fn ion_picks_the_same_water_as_gromospp() {
    let topo = shared().join("aladip.topo");
    let pos = shared().join("aladip.conf");
    let out = run(
        env!("CARGO_BIN_EXE_ion"),
        &[
            "--topo",
            &topo.to_string_lossy(),
            "--pbc",
            "r",
            "--pos",
            &pos.to_string_lossy(),
            "--positive",
            "1",
            "NA+",
            "--potential",
            "1.4",
        ],
    );
    let (ours, _) = parse(&out);
    let (theirs, _) =
        parse(&std::fs::read_to_string(data("ion").join("aladip_na.gromospp.cnf")).unwrap());
    assert_eq!(ours.len(), theirs.len(), "atom count");
    let na = |v: &[(String, String, [f64; 3])]| {
        v.iter()
            .find(|a| a.1 == "NA+")
            .map(|a| a.2)
            .expect("an NA+ atom")
    };
    let (a, b) = (na(&ours), na(&theirs));
    for d in 0..3 {
        assert!((a[d] - b[d]).abs() < 1e-8, "NA+ at {a:?}, gromos++ {b:?}");
    }
}
