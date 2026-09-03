//! Free energy trajectory (.trg) — gromosXX's `TIMESTEP` + `FREEENERDERIVS03` layout.
//!
//! The writer emits what the native binary emits (`out_configuration.cc::_print_free_energyred`
//! via `_print_energyred_helper`): a `TIMESTEP` block, then `FREEENERDERIVS03` with `# lambda`,
//! the `# totals` of the derivative in ENERGY03 order (total, kinetic, potential, covalent,
//! bond, angle, improper, dihedral, cross-dihedral, nonbonded, LJ, CRF, then the lattice-sum,
//! special and restraint slots), and the per-bath / per-energy-group sections with the counts
//! of the run — so gromos++ `ene_ana`/`ext_ti_ana` read the file, and the reference suite
//! compares dH/dλ per term. Per-group derivatives are written as zeros (the engine tracks
//! totals per term, not per energy group); numbers are printed with Rust's `e` exponent
//! (`1.5e-2`, not `1.500000000e-02`), which every GROMOS reader parses.
//!
//! The reader accepts this layout and the one-line `FREEENERGY03` block older files carry.
//!
//! dH/dλ (kJ/mol) is the thermodynamic-integration integrand at the current λ. Integrate
//! ⟨dH/dλ⟩_λ over λ with `ext_ti_ana` to get ΔG.

use crate::IoError;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;

/// One frame of free-energy derivative data.
#[derive(Debug, Clone, Default)]
pub struct FreeEnergyFrame {
    /// MD step number (the `TIMESTEP` block).
    pub step: usize,
    pub time: f64,
    pub lambda: f64,
    // Per-component dH/dλ (kJ/mol); zero if not computed separately.
    pub dhdl_bond: f64,
    pub dhdl_angle: f64,
    pub dhdl_improper: f64,
    pub dhdl_dihedral: f64,
    pub dhdl_lj: f64,
    pub dhdl_crf: f64,
    pub dhdl_special: f64,
    /// Total dH/dλ — the quantity consumed by ext_ti_ana / TI integration.
    pub dhdl_total: f64,
}

impl FreeEnergyFrame {
    pub fn new(time: f64, lambda: f64, dhdl_total: f64) -> Self {
        Self {
            time,
            lambda,
            dhdl_total,
            ..Default::default()
        }
    }
}

/// Read every frame of a `.trg` file, in either layout:
///
/// * gromosXX's `FREEENERDERIVS03` (`# lambda`, then `# totals` in ENERGY03 order — total,
///   kinetic, potential, covalent, bond, angle, improper, dihedral, cross-dihedral, nonbonded,
///   LJ, CRF, …; the time comes from the preceding `TIMESTEP` block), which is what the
///   native binary writes and what this crate's writer now emits;
/// * the one-line `FREEENERGY03` block this crate wrote before 0.0.33 (time, λ, per-term
///   dH/dλ, total), so older files keep reading.
///
/// Tolerates extra comment lines and blank lines. Returns frames in order.
pub fn read_free_energy_trajectory<P: AsRef<Path>>(
    path: P,
) -> Result<Vec<FreeEnergyFrame>, IoError> {
    let reader = crate::gz::open_text(path.as_ref())
        .map_err(|e| IoError::FileNotFound(format!("{}: {e}", path.as_ref().display())))?;

    #[derive(PartialEq)]
    enum Block {
        None,
        Timestep,
        Legacy,
        Derivs,
    }
    let mut frames = Vec::new();
    let mut block = Block::None;
    let mut section = String::new();
    let mut step = 0usize;
    let mut time = 0.0;
    let mut lambda = 0.0;
    let mut totals: Vec<f64> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(IoError::Io)?;
        let t = line.trim();
        match t {
            "TIMESTEP" => {
                block = Block::Timestep;
                continue;
            },
            "FREEENERGY03" => {
                block = Block::Legacy;
                continue;
            },
            "FREEENERDERIVS03" => {
                block = Block::Derivs;
                section.clear();
                totals.clear();
                continue;
            },
            "END" => {
                if block == Block::Derivs {
                    let at = |i: usize| totals.get(i).copied().unwrap_or(0.0);
                    frames.push(FreeEnergyFrame {
                        step,
                        time,
                        lambda,
                        dhdl_total: at(0),
                        dhdl_bond: at(4),
                        dhdl_angle: at(5),
                        dhdl_improper: at(6),
                        dhdl_dihedral: at(7),
                        dhdl_lj: at(10),
                        dhdl_crf: at(11),
                        dhdl_special: 0.0,
                    });
                }
                block = Block::None;
                continue;
            },
            _ => {},
        }
        if t.is_empty() {
            continue;
        }
        match block {
            Block::Timestep => {
                let vals: Vec<f64> = t
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 2 {
                    step = vals[0] as usize;
                    time = vals[1];
                }
            },
            Block::Legacy => {
                if t.starts_with('#') {
                    continue;
                }
                let vals: Vec<f64> = t
                    .split_whitespace()
                    .filter_map(|s| s.parse().ok())
                    .collect();
                if vals.len() >= 9 {
                    frames.push(FreeEnergyFrame {
                        step: 0,
                        time: vals[0],
                        lambda: vals[1],
                        dhdl_bond: vals[2],
                        dhdl_angle: vals[3],
                        dhdl_improper: vals[4],
                        dhdl_dihedral: vals[5],
                        dhdl_lj: vals[6],
                        dhdl_crf: vals[7],
                        dhdl_special: 0.0,
                        dhdl_total: vals[8],
                    });
                } else if vals.len() == 3 {
                    frames.push(FreeEnergyFrame::new(vals[0], vals[1], vals[2]));
                }
            },
            Block::Derivs => {
                if let Some(name) = t.strip_prefix('#') {
                    section = name.trim().to_string();
                    continue;
                }
                match section.as_str() {
                    "lambda" => {
                        if let Some(v) = t.split_whitespace().next().and_then(|s| s.parse().ok()) {
                            lambda = v;
                        }
                    },
                    "totals" => {
                        totals.extend(t.split_whitespace().filter_map(|s| s.parse::<f64>().ok()))
                    },
                    _ => {},
                }
            },
            Block::None => {},
        }
    }
    Ok(frames)
}

/// Writer for gromosXX-layout free-energy trajectory files (.trg).
pub struct FreeEnergyWriter {
    writer: BufWriter<File>,
    n_baths: usize,
    n_groups: usize,
}

/// Number of entries in a `# totals` section (gromosXX `Energy` totals, 2023-04-15 layout).
const N_TOTALS: usize = 52;

impl FreeEnergyWriter {
    /// Create the file with its `TITLE` and `ENEVERSION` blocks. The per-bath / per-group
    /// sections default to one bath and one energy group; set the run's counts with
    /// [`with_layout`](Self::with_layout).
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "TITLE")?;
        writeln!(writer, "\t{}", title)?;
        writeln!(writer, "\n\tfree energy trajectory")?;
        writeln!(writer, "END")?;
        writeln!(writer, "ENEVERSION")?;
        writeln!(writer, "\t2023-04-15")?;
        writeln!(writer, "END")?;
        Ok(Self {
            writer,
            n_baths: 1,
            n_groups: 1,
        })
    }

    /// The number of temperature baths and energy groups of the run (sizes of the per-bath
    /// and per-group sections, as in the native file).
    pub fn with_layout(mut self, n_baths: usize, n_groups: usize) -> Self {
        self.n_baths = n_baths.max(1);
        self.n_groups = n_groups.max(1);
        self
    }

    /// Append one frame: a `TIMESTEP` block and a `FREEENERDERIVS03` block.
    pub fn write_frame(&mut self, frame: &FreeEnergyFrame) -> Result<(), IoError> {
        let w = &mut self.writer;
        writeln!(w, "TIMESTEP")?;
        writeln!(w, "{:>15}{:>15.9}", frame.step, frame.time)?;
        writeln!(w, "END")?;
        writeln!(w, "FREEENERDERIVS03")?;
        writeln!(w, "# lambda")?;
        writeln!(w, "{:>18.9e}", frame.lambda)?;
        writeln!(w, "# totals")?;
        let bonded = frame.dhdl_bond + frame.dhdl_angle + frame.dhdl_improper + frame.dhdl_dihedral;
        let nonbonded = frame.dhdl_lj + frame.dhdl_crf;
        let mut totals = [0.0f64; N_TOTALS];
        totals[0] = frame.dhdl_total;
        totals[2] = frame.dhdl_total; // potential: no kinetic derivative is tracked
        totals[3] = bonded;
        totals[4] = frame.dhdl_bond;
        totals[5] = frame.dhdl_angle;
        totals[6] = frame.dhdl_improper;
        totals[7] = frame.dhdl_dihedral;
        totals[9] = nonbonded;
        totals[10] = frame.dhdl_lj;
        totals[11] = frame.dhdl_crf;
        totals[18] = frame.dhdl_special;
        for v in totals {
            writeln!(w, "{:>18.9e}", v)?;
        }
        let zeros = |n: usize| {
            std::iter::repeat_n("   0.000000000e0", n)
                .collect::<Vec<_>>()
                .join("")
        };
        writeln!(w, "# baths")?;
        writeln!(w, "{}", self.n_baths)?;
        for _ in 0..self.n_baths {
            writeln!(w, "{}", zeros(3))?;
        }
        writeln!(w, "# bonded")?;
        writeln!(w, "{}", self.n_groups)?;
        for _ in 0..self.n_groups {
            writeln!(w, "{}", zeros(5))?;
        }
        writeln!(w, "# nonbonded")?;
        for _ in 0..(self.n_groups * (self.n_groups + 1) / 2) {
            writeln!(w, "{}", zeros(6))?;
        }
        writeln!(w, "# special")?;
        for _ in 0..self.n_groups {
            writeln!(w, "{}", zeros(13))?;
        }
        writeln!(w, "# eds\n# numstates\n0")?;
        writeln!(
            w,
            "           # total         nonbonded           special            offset        total_orig       total_phys"
        )?;
        writeln!(w, "# gamd\n# numaccelgroups\n0")?;
        writeln!(
            w,
            "      # E_dihedral       E_potential        K_dihedral       K_potential         dV_group"
        )?;
        writeln!(w, "# precalclam\n# nr_lambdas\n0")?;
        writeln!(
            w,
            "          # A_e_lj            B_e_lj           A_e_crf           B_e_crf        AB_kinetic           AB_bond          AB_angle       AB_improper         AB_disres         AB_angres         AB_dihres        AB_disfld"
        )?;
        writeln!(w, "# ABdih\n      # A_dihedral       B_dihedral")?;
        writeln!(w, "{}", zeros(2))?;
        writeln!(w, "END")?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush().map_err(IoError::Io)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_survive_write_then_read() {
        let p = std::env::temp_dir().join(format!("gromos_trg_test_{}.trg", std::process::id()));
        {
            let mut w = FreeEnergyWriter::new(&p, "round trip")
                .unwrap()
                .with_layout(2, 3);
            let mut f = FreeEnergyFrame::new(0.002, 0.25, -12.5);
            f.dhdl_lj = -10.0;
            f.dhdl_crf = -2.5;
            w.write_frame(&f).unwrap();
            w.write_frame(&FreeEnergyFrame::new(0.004, 0.25, -11.0))
                .unwrap();
            w.flush().unwrap();
        }
        let frames = read_free_energy_trajectory(&p).unwrap();
        assert_eq!(frames.len(), 2);
        assert!((frames[0].time - 0.002).abs() < 1e-12);
        assert!((frames[0].lambda - 0.25).abs() < 1e-12);
        assert!((frames[0].dhdl_total + 12.5).abs() < 1e-9);
        assert!((frames[0].dhdl_lj + 10.0).abs() < 1e-9);
        assert!((frames[1].dhdl_total + 11.0).abs() < 1e-9);
        std::fs::remove_file(p).ok();
    }
}
