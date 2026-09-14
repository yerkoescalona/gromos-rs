//! ext_ti_ana — Thermodynamic Integration analysis (gromos-rs-compatible).
//!
//! Reads one `.trg` file per λ window, computes ⟨dH/dλ⟩ ± ee() for each,
//! integrates with the trapezoidal rule, and reports ΔG.
//!
//! Usage (GROMOS style):
//!   ext_ti_ana @trg run_lam0.trg run_lam1.trg ... @lambda 0.0 0.25 0.5 0.75 1.0
//!
//!   @trg <file> [<file> ...]   One .trg file per λ window (in λ order).
//!   @lambda <val> [<val> ...]  Explicit λ values; must match number of @trg files.
//!                              If omitted, λ is read from the first frame of each file.
//!   @skip <n>                  Discard first n frames of each window (equilibration).
//!   @verb                      Verbose: print per-frame dH/dλ.

use crate::programs::io_err;
use gromos_core::stat::Stat;
use gromos_io::args::{Arguments, Program, Status};
use gromos_io::read_free_energy_trajectory;
use std::io::Write;

pub const USAGE: &str = "# ext_ti_ana
\t@trg           <.trg files, one per lambda window, in ascending lambda order>
\t[@lambda       <explicit lambda values (default: read from the file headers)>]
\t[@skip         <skip the first n frames of each window (default: 0)>]
\t[@verb         <print the per-window frame listing>]";

pub const KNOWN: &[&str] = &["trg", "lambda", "skip", "verb"];

pub const PROGRAM: Program = Program {
    name: "ext_ti_ana",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 3 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let trg_files: Vec<String> = args.values("trg").to_vec();
    let explicit_lambdas: Vec<f64> = args
        .values("lambda")
        .iter()
        .map(|v| {
            v.parse()
                .map_err(|_| format!("@lambda: cannot parse '{v}'"))
        })
        .collect::<Result<_, _>>()?;
    let skip: usize = args.get("skip", 0)?;
    let verbose = args.has("verb");

    if trg_files.is_empty() {
        return Err(format!("no @trg files specified\n{USAGE}"));
    }

    // Read each .trg file
    struct Window {
        lambda: f64,
        stat: Stat,
    }

    let mut windows: Vec<Window> = Vec::new();

    for (wi, path) in trg_files.iter().enumerate() {
        let frames = match read_free_energy_trajectory(path) {
            Ok(traj) => {
                for w in traj.warnings {
                    eprintln!("{w}");
                }
                traj.frames
            },
            Err(e) => return Err(format!("Error reading {path}: {e}")),
        };

        if frames.is_empty() {
            eprintln!("Warning: {path} contains no free-energy frames (FREEENERDERIVS03 / FREEENERGY03) — skipping");
            continue;
        }

        let lambda = if let Some(&lam) = explicit_lambdas.get(wi) {
            lam
        } else {
            frames[0].lambda
        };

        let mut stat = Stat::new();
        let n_skip = skip.min(frames.len());
        for frame in &frames[n_skip..] {
            stat.add(frame.dhdl_total);
        }

        if verbose {
            writeln!(
                out,
                "# λ={lambda:.6}  file={path}  frames={} (skipped={n_skip})",
                frames.len()
            )
            .map_err(io_err)?;
            for (fi, frame) in frames[n_skip..].iter().enumerate() {
                writeln!(
                    out,
                    "  frame {:5}  t={:.4}ps  dH/dλ={:.6e}",
                    fi, frame.time, frame.dhdl_total
                )
                .map_err(io_err)?;
            }
        }

        windows.push(Window { lambda, stat });
    }

    if windows.is_empty() {
        return Err("no valid λ windows loaded".into());
    }

    // Sort by λ
    windows.sort_by(|a, b| a.lambda.partial_cmp(&b.lambda).unwrap());

    // Print per-window summary
    writeln!(out, "#").map_err(io_err)?;
    writeln!(out, "# ext_ti_ana — Thermodynamic Integration").map_err(io_err)?;
    writeln!(
        out,
        "# {:>10}  {:>14}  {:>14}  {:>10}",
        "lambda", "<dH/dλ>", "ee(dH/dλ)", "N"
    )
    .map_err(io_err)?;
    writeln!(out, "#").map_err(io_err)?;

    let mut lambdas: Vec<f64> = Vec::new();
    let mut aves: Vec<f64> = Vec::new();

    for w in &windows {
        let n = w.stat.n();
        if n == 0 {
            eprintln!("Warning: window λ={:.6} has 0 samples", w.lambda);
            continue;
        }
        let ave = w.stat.ave();
        let ee = if n >= 4 { w.stat.ee() } else { w.stat.rmsd() };
        writeln!(
            out,
            "  {:10.6}  {:14.6e}  {:14.6e}  {:10}",
            w.lambda, ave, ee, n
        )
        .map_err(io_err)?;
        lambdas.push(w.lambda);
        aves.push(ave);
    }

    if lambdas.len() < 2 {
        writeln!(out, "#").map_err(io_err)?;
        writeln!(out, "# Need ≥2 λ windows for integration").map_err(io_err)?;
        return Ok(());
    }

    // Trapezoidal integration: ΔG = ∫⟨dH/dλ⟩dλ
    let mut delta_g = 0.0f64;
    for k in 0..lambdas.len() - 1 {
        let dl = lambdas[k + 1] - lambdas[k];
        delta_g += 0.5 * (aves[k] + aves[k + 1]) * dl;
    }

    writeln!(out, "#").map_err(io_err)?;
    writeln!(
        out,
        "# TI integration (trapezoidal):  ΔG = {delta_g:.6} kJ/mol"
    )
    .map_err(io_err)?;
    writeln!(
        out,
        "#   λ range: {:.4} – {:.4}  ({} windows)",
        lambdas[0],
        lambdas[lambdas.len() - 1],
        lambdas.len()
    )
    .map_err(io_err)?;
    Ok(())
}
