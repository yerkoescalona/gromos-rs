//! Free energy trajectory writer (.trg) — gromosXX FREEENERGY03 block format.
//!
//! One FREEENERGY03 block per frame, containing time, lambda, and dH/dλ components.
//!
//! # Format (gromosXX out_energy.cc::_print_free_energy_trajectory)
//! ```text
//! TITLE
//!   GROMOS-RS free energy trajectory
//! END
//! FREEENERGY03
//!   time[ps]  lambda  dH/dL_bond  dH/dL_angle  dH/dL_improper  dH/dL_dihedral
//!             dH/dL_lj  dH/dL_crf  dH/dL_special  dH/dL_total
//! END
//! ```
//!
//! dH/dλ (kJ/mol) is the thermodynamic-integration integrand at the current λ.
//! Integrate ⟨dH/dλ⟩_λ over λ with ext_ti_ana to get ΔG.

use crate::IoError;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// One frame of free-energy derivative data.
#[derive(Debug, Clone, Default)]
pub struct FreeEnergyFrame {
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
        Self { time, lambda, dhdl_total, ..Default::default() }
    }
}

/// Read all FREEENERGY03 frames from a `.trg` file.
///
/// Tolerates extra comment lines and blank lines.  Returns frames in order.
pub fn read_free_energy_trajectory<P: AsRef<Path>>(path: P) -> Result<Vec<FreeEnergyFrame>, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|e| IoError::FileNotFound(format!("{}: {e}", path.as_ref().display())))?;
    let reader = BufReader::new(file);

    let mut frames = Vec::new();
    let mut in_block = false;

    for line in reader.lines() {
        let line = line.map_err(IoError::Io)?;
        let t = line.trim();
        if t == "FREEENERGY03" {
            in_block = true;
            continue;
        }
        if t == "END" && in_block {
            in_block = false;
            continue;
        }
        if in_block && !t.starts_with('#') && !t.is_empty() {
            let vals: Vec<f64> = t.split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            // columns: time lambda dhdl_bond dhdl_angle dhdl_improper dhdl_dihedral dhdl_lj dhdl_crf dhdl_total
            if vals.len() >= 9 {
                frames.push(FreeEnergyFrame {
                    time:          vals[0],
                    lambda:        vals[1],
                    dhdl_bond:     vals[2],
                    dhdl_angle:    vals[3],
                    dhdl_improper: vals[4],
                    dhdl_dihedral: vals[5],
                    dhdl_lj:       vals[6],
                    dhdl_crf:      vals[7],
                    dhdl_special:  0.0,
                    dhdl_total:    vals[8],
                });
            }
        }
    }
    Ok(frames)
}

/// Writer for gromosXX-compatible free-energy trajectory files (.trg).
///
/// Each call to `write_frame` appends one FREEENERGY03 block.
pub struct FreeEnergyWriter {
    writer: BufWriter<File>,
}

impl FreeEnergyWriter {
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "TITLE")?;
        writeln!(writer, "  {}", title)?;
        writeln!(writer, "END")?;
        Ok(Self { writer })
    }

    /// Append one FREEENERGY03 block for the given frame.
    pub fn write_frame(&mut self, frame: &FreeEnergyFrame) -> Result<(), IoError> {
        writeln!(self.writer, "FREEENERGY03")?;
        writeln!(
            self.writer,
            " {:20.9e} {:20.9e} {:20.9e} {:20.9e} {:20.9e} {:20.9e} {:20.9e} {:20.9e} {:20.9e}",
            frame.time,
            frame.lambda,
            frame.dhdl_bond,
            frame.dhdl_angle,
            frame.dhdl_improper,
            frame.dhdl_dihedral,
            frame.dhdl_lj,
            frame.dhdl_crf,
            frame.dhdl_total,
        )?;
        writeln!(self.writer, "END")?;
        Ok(())
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush().map_err(IoError::Io)
    }
}
