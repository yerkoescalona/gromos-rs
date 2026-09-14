//! bar — Bennett Acceptance Ratio (BAR) free-energy estimator.
//!
//! Computes ΔG between two states i and j using the BAR method
//! (Bennett 1976, Shirts et al. 2003).
//!
//! Input: two text files, each with one energy difference per line (kJ/mol):
//!   - forward: ΔE = E_j(x) − E_i(x)   sampled from state i
//!   - backward: ΔE = E_i(x) − E_j(x)  sampled from state j
//!
//! The iteration uses the numerically stable log-sum-exp form from
//! Shirts et al. Phys. Rev. Lett. 91 (2003) 140601, matching the GROMOS
//! implementation in bar.cc.
//!
//! Usage:
//!   bar @forward <file> @backward <file> [@temp <K>]
//!       [@maxiter <n>] [@conv <eps>] [@bootstrap <n>]

use crate::programs::io_err;
use gromos_core::stat::Stat;
use gromos_core::units::kB;
use gromos_io::args::{Arguments, Program, Status};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader};

pub const USAGE: &str = "# bar
\t@forward       <energy differences E_j(x)-E_i(x) from the state-i samples>
\t@backward      <energy differences E_i(x)-E_j(x) from the state-j samples>
\t[@temp         <temperature in K (default: 300)>]
\t[@maxiter      <maximum BAR iterations (default: 500)>]
\t[@conv         <convergence threshold in kBT (default: 1e-5)>]
\t[@bootstrap    <bootstrap error-estimate iterations (default: 0)>]";

pub const KNOWN: &[&str] = &[
    "forward",
    "backward",
    "temp",
    "maxiter",
    "conv",
    "bootstrap",
];

pub const PROGRAM: Program = Program {
    name: "bar",
    known: KNOWN,
    usage: USAGE,
    status: Status::Ported,
    run,
};

/// Numerically stable BAR iteration (log-sum-exp form, Shirts 2003 / GROMOS bar.cc).
///
/// Returns the updated ΔG estimate in units of kBT.
fn bar_step(
    dw_fwd: &[f64], // (E_j - E_i)/kBT for samples from i
    dw_bwd: &[f64], // (E_i - E_j)/kBT for samples from j
    dg: f64,        // current ΔG estimate in kBT units
) -> f64 {
    let n_i = dw_fwd.len() as f64;
    let n_j = dw_bwd.len() as f64;
    let m = (n_i / n_j).ln(); // log ratio of sample counts

    // max of forward energy differences (shift for numerical stability)
    let max_ij = dw_fwd.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let max_ji = dw_bwd.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    // Pre-compute exponentials (constant across iterations)
    let exp_ij: Vec<f64> = dw_fwd.iter().map(|&e| (e - max_ij).exp()).collect();
    let exp_ji: Vec<f64> = dw_bwd.iter().map(|&e| (e - max_ji).exp()).collect();

    // Forward: ln⟨f(x_i)⟩_i  with x_i = dw_fwd - dg + m
    let mx_ij = max_ij - dg + m;
    let exp_mx_ij = (-mx_ij).exp();
    let lnsum_ij = {
        let s: f64 = exp_ij.iter().map(|&e| (exp_mx_ij + e).ln()).sum::<f64>();
        -mx_ij - s / n_i
    };

    // Backward: ln⟨f(x_j)⟩_j  with x_j = dw_bwd - dg - m  (sign flip from xj identity)
    let mx_ji = max_ji + dg - m;
    let exp_mx_ji = (-mx_ji).exp();
    let lnsum_ji = {
        let sum_log: f64 = dw_bwd
            .iter()
            .zip(&exp_ji)
            .map(|(&e_orig, &e)| -(mx_ji + (exp_mx_ji + e).ln()) - e_orig)
            .sum::<f64>();
        sum_log / n_j
    };

    // New estimate: DG = (lnsum_ij - lnsum_ji) + m
    lnsum_ij - lnsum_ji + m
}

fn read_values(path: &str) -> Result<Vec<f64>, String> {
    let file = File::open(path).map_err(|e| format!("Error opening '{path}': {e}"))?;
    Ok(BufReader::new(file)
        .lines()
        .map_while(Result::ok)
        .filter(|l| !l.trim_start().starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| l.trim().parse::<f64>().ok())
        .collect())
}

pub fn run(argv: &[String], out: &mut dyn Write) -> Result<(), String> {
    if argv.len() < 2 {
        return Err(USAGE.to_string());
    }
    let args = Arguments::from_vec(argv.to_vec(), KNOWN, USAGE)?;
    let fwd_path = args.value("forward")?;
    let bwd_path = args.value("backward")?;
    let temp: f64 = args.get("temp", 300.0)?;
    let max_iter: usize = args.get("maxiter", 500)?;
    let conv_eps: f64 = args.get("conv", 1e-5)?;
    let bootstrap: usize = args.get("bootstrap", 0)?;

    let kbt = kB * temp;
    let beta = 1.0 / kbt;

    // Read raw energy differences (kJ/mol) and convert to kBT units
    let dw_fwd_kj = read_values(fwd_path)?;
    let dw_bwd_kj = read_values(bwd_path)?;

    if dw_fwd_kj.is_empty() || dw_bwd_kj.is_empty() {
        return Err("empty input files".into());
    }

    let dw_fwd: Vec<f64> = dw_fwd_kj.iter().map(|&e| e * beta).collect();
    let dw_bwd: Vec<f64> = dw_bwd_kj.iter().map(|&e| e * beta).collect();

    eprintln!(
        "# BAR: T={temp} K  N_fwd={}  N_bwd={}",
        dw_fwd.len(),
        dw_bwd.len()
    );

    // Initial guess: half the mean forward energy difference
    let mut dg = dw_fwd.iter().sum::<f64>() / dw_fwd.len() as f64 / 2.0;
    let mut converged = false;

    for k in 0..max_iter {
        let dg_new = bar_step(&dw_fwd, &dw_bwd, dg);
        let reldiff = (dg_new - dg).abs();
        dg = dg_new;
        if reldiff < conv_eps {
            eprintln!("# Converged after {} iterations (Δ={reldiff:.2e})", k + 1);
            converged = true;
            break;
        }
    }
    if !converged {
        eprintln!("# Warning: did not converge within {max_iter} iterations");
    }

    let dg_kj = dg * kbt;

    // Bootstrap error estimate
    let error_kj = if bootstrap > 0 {
        let mut boot_stat = Stat::new();
        let n_i = dw_fwd.len();
        let n_j = dw_bwd.len();
        let mut rng_state = 12345u64;

        for _ in 0..bootstrap {
            // Simple LCG random resampling
            let sample_fwd: Vec<f64> = (0..n_i)
                .map(|_| {
                    rng_state = rng_state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    dw_fwd[(rng_state >> 33) as usize % n_i]
                })
                .collect();
            let sample_bwd: Vec<f64> = (0..n_j)
                .map(|_| {
                    rng_state = rng_state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    dw_bwd[(rng_state >> 33) as usize % n_j]
                })
                .collect();

            let mut dg_b = dg;
            for _ in 0..max_iter {
                let dg_new = bar_step(&sample_fwd, &sample_bwd, dg_b);
                if (dg_new - dg_b).abs() < conv_eps {
                    dg_b = dg_new;
                    break;
                }
                dg_b = dg_new;
            }
            boot_stat.add(dg_b * kbt);
        }
        boot_stat.rmsd()
    } else {
        0.0
    };

    writeln!(out, "# BAR free-energy estimate").map_err(io_err)?;
    writeln!(out, "# T = {temp:.1} K   kBT = {kbt:.4} kJ/mol").map_err(io_err)?;
    writeln!(out, "# N_fwd = {}   N_bwd = {}", dw_fwd.len(), dw_bwd.len()).map_err(io_err)?;
    if bootstrap > 0 {
        writeln!(out, "#").map_err(io_err)?;
        writeln!(
            out,
            "# ΔG = {dg_kj:.6} ± {error_kj:.6} kJ/mol  (bootstrap N={bootstrap})"
        )
        .map_err(io_err)?;
    }
    writeln!(out, "#").map_err(io_err)?;
    writeln!(out, "ΔG = {dg_kj:.6} kJ/mol").map_err(io_err)?;
    if bootstrap > 0 {
        writeln!(out, "err = {error_kj:.6} kJ/mol").map_err(io_err)?;
    }
    Ok(())
}
