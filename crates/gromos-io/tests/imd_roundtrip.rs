//! `.imd` read → write → read must be lossless (PLAN.md 3.9 G3, step 2).
//!
//! For every gromosXX reference input: parse, render with `write_imd`, parse again, and
//! compare every modelled field exactly. Blocks the parser does not model must survive
//! verbatim as passthrough.

use std::path::{Path, PathBuf};

use gromos_io::imd::{parse_imd_str, read_imd_file, ImdParameters};
use gromos_io::imd_write::{write_imd, MODELLED_BLOCKS};

fn refs() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../gromos-md/tests/gromosXX_references")
}

fn reference_inputs() -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = std::fs::read_dir(refs())
        .expect("reference dir")
        .flatten()
        .filter(|e| e.path().is_dir())
        .flat_map(|dir| {
            std::fs::read_dir(dir.path())
                .into_iter()
                .flatten()
                .flatten()
                .map(|e| e.path())
                .filter(|p| p.extension().is_some_and(|x| x == "in"))
        })
        .collect();
    files.sort();
    files
}

/// Modelled fields only: the raw copies of modelled blocks are regenerated on write and
/// legitimately differ in whitespace, so the comparison keeps only passthrough blocks.
fn comparable(p: &ImdParameters) -> ImdParameters {
    let mut q = p.clone();
    q.raw_blocks
        .retain(|name, _| !MODELLED_BLOCKS.contains(&name.as_str()));
    q
}

#[test]
fn every_reference_input_round_trips_exactly() {
    let inputs = reference_inputs();
    assert!(
        inputs.len() >= 40,
        "found only {} reference inputs",
        inputs.len()
    );
    for path in &inputs {
        let original = read_imd_file(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let text = write_imd(&original, None);
        let reparsed = parse_imd_str(&text).unwrap_or_else(|e| {
            panic!(
                "{}: re-parse of written text failed: {e}\n{text}",
                path.display()
            )
        });
        assert_eq!(
            comparable(&reparsed),
            comparable(&original),
            "{}: read → write → read changed a field\n--- written ---\n{text}",
            path.display()
        );
        // Writing the re-parsed parameters must be a fixed point.
        assert_eq!(
            write_imd(&reparsed, None),
            text,
            "{}: write is not idempotent",
            path.display()
        );
    }
}

#[test]
fn absent_optional_blocks_mean_gromosxx_defaults() {
    // No MULTIBATH → no bath (gromosXX: "Adding a bath, no temperature coupling");
    // no COMTRANSROT → no centre-of-mass motion removal (parameter.h: skip_step 0).
    let minimal = "\
TITLE
minimal
END
SYSTEM
#      NPM      NSM
         1        0
END
STEP
#   NSTLIM         T        DT
       10       0.0     0.002
END
BOUNDCOND
#      NTB    NDFMIN
         1         0
END
";
    let p = parse_imd_str(minimal).expect("parse");
    assert!(
        p.temp_bath.is_empty(),
        "absent MULTIBATH must not create a bath"
    );
    assert_eq!(p.nscm, 0, "absent COMTRANSROT must not remove COM motion");
    assert!(!p.couple_pressure);
    assert_eq!(p.ntg, 0);
    assert_eq!(p.ntem, 0);
    // The defaults table itself agrees with the absent-block meaning.
    let d = ImdParameters::default();
    assert!(d.temp_bath.is_empty());
    assert_eq!(d.nscm, 0);
}

#[test]
fn garbage_numbers_are_errors_not_zeros() {
    let bad = "\
STEP
#   NSTLIM         T        DT
       ten       0.0     0.002
END
";
    let err = parse_imd_str(bad).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("STEP"), "error must name the block: {msg}");
    assert!(
        msg.contains("ten"),
        "error must show the offending token: {msg}"
    );
}

#[test]
fn unmodelled_blocks_pass_through_verbatim() {
    let text = "\
SYSTEM
#      NPM      NSM
         1        0
END
GAMD
# a block this parser does not model
   1  2  3
END
";
    let p = parse_imd_str(text).expect("parse");
    let written = write_imd(&p, None);
    assert!(
        written.contains("GAMD\n# a block this parser does not model\n1  2  3\nEND\n"),
        "passthrough block missing or altered:\n{written}"
    );
    let q = parse_imd_str(&written).expect("re-parse");
    assert_eq!(q.raw_blocks.get("GAMD"), p.raw_blocks.get("GAMD"));
}

#[test]
fn multibath_without_dof_sets_is_synthesised_from_n_atoms() {
    let mut p = ImdParameters::nvt(0.002, 10, 300.0, "none").expect("factory");
    assert!(p.temp_bath[0].dof_sets.is_empty());
    let written = write_imd(&p, Some(648));
    assert!(written.contains("648"), "{written}");
    let q = parse_imd_str(&written).expect("re-parse");
    assert_eq!(q.temp_bath[0].dof_sets, vec![[648, 1, 1]]);
    // …and one energy group spanning the system in FORCE.
    p.nre.clear();
    let written = write_imd(&p, Some(648));
    let q = parse_imd_str(&written).expect("re-parse");
    assert_eq!((q.negr, q.nre.clone()), (1, vec![648]));
}
