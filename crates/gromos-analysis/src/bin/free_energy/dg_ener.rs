//! dg_ener — free-energy difference by exponential averaging of ΔE (gromos++ `dg_ener`).
//!
//! Usage: dg_ener @temp <T> @stateA <file> @stateB <file> [@col <colA colB>]
//!
//! Both files are time series (time in column 1); the energy is the last column unless `@col`
//! names it (2-based). Prints the running −kT ln<exp(−βΔE)> and the final ΔG_BA = G_B − G_A.

use gromos_analysis::args::{fail, Arguments};
use gromos_analysis::distribution::cpp_g;
use gromos_io::table::read_columns;

const USAGE: &str = "# dg_ener
\t@temp <temperature for perturbation>
\t@stateA <energy file for state A>
\t@stateB <energy file for state B>
\t@col <numbers of the columns to use from file A and B [default: last]>";

fn main() {
    if let Err(e) = run() {
        fail(&e);
    }
}

fn run() -> Result<(), String> {
    let args = Arguments::parse(&["temp", "stateA", "stateB", "col"], USAGE)?;
    let temp: f64 = args.require("temp")?;
    let a = read_columns(
        args.value("stateA")
            .map_err(|_| "energy file for state A missing".to_string())?,
    )
    .map_err(|e| format!("state A: {e}"))?;
    let b = read_columns(
        args.value("stateB")
            .map_err(|_| "energy file for state B missing".to_string())?,
    )
    .map_err(|e| format!("state B: {e}"))?;
    let cols = args.values("col");
    let (col_a, col_b) = if cols.is_empty() {
        (0usize, 0usize)
    } else {
        let ca: usize = cols.first().and_then(|s| s.parse().ok()).unwrap_or(0);
        let cb: usize = cols.get(1).and_then(|s| s.parse().ok()).unwrap_or(ca);
        if ca < 2 || cb < 2 {
            return Err(
                "@col numbers have to be >1, the first column is reserved for the time".into(),
            );
        }
        (ca, cb)
    };
    if a.len() != b.len() {
        return Err("Error while reading file for state A or state B: check if number of lines are identical and number of columns are correct.".into());
    }
    let pick = |row: &Vec<f64>, col: usize, what: &str| -> Result<f64, String> {
        if row.len() < 2 {
            return Err(format!("State{what} file does not have enough columns."));
        }
        if col == 0 {
            Ok(*row.last().unwrap())
        } else {
            row.get(col - 1)
                .copied()
                .ok_or_else(|| format!("State{what} file does not have enough columns."))
        }
    };
    let mut delta_v = Vec::with_capacity(a.len());
    let mut t = Vec::with_capacity(a.len());
    let mut time_warning = false;
    for (ra, rb) in a.iter().zip(&b) {
        if ra[0] != rb[0] {
            time_warning = true;
        }
        delta_v.push(pick(rb, col_b, "B")? - pick(ra, col_a, "A")?);
        t.push(ra[0]);
    }
    if delta_v.is_empty() {
        return Err("no data".into());
    }
    let beta = 1.0 / (gromos_analysis::GROMOSPP_BOLTZMANN * temp);
    let mut out = String::new();
    out.push_str(&format!(
        "# Time{:>12}{:>12}{:>12}\n",
        "DE_tot", "probability", "DG_BA"
    ));
    let mut sum = -delta_v[0] * beta;
    let mut dg = 0.0;
    for i in 1..delta_v.len() {
        let p = -delta_v[i] * beta;
        sum = sum.max(p) + (1.0 + (sum.min(p) - sum.max(p)).exp()).ln();
        dg = -(sum - ((i + 1) as f64).ln()) / beta;
        out.push_str(&format!(
            "{:>6}{:>12}{:>12}{:>12}\n",
            cpp_g(t[i], 5),
            cpp_g(delta_v[i], 5),
            cpp_g(p.exp(), 5),
            cpp_g(dg, 5)
        ));
    }
    out.push_str(&format!(
        "# final result dG_BA = G_B - G_A: {}\n",
        cpp_g(dg, 5)
    ));
    if time_warning {
        out.push_str("#\n# WARNING: time was not equal in state A and state B!\n");
        eprintln!("\nWARNING: time was not equal in state A and state B!");
    }
    print!("{out}");
    Ok(())
}
