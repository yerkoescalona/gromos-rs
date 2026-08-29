//! GROMOS energy trajectory (.tre) — gromosXX's `TIMESTEP` + `ENERGY03` + `VOLUMEPRESSURE03` layout.
//!
//! The writer emits what the native binary emits (`out_configuration.cc::_print_energyred` via
//! `_print_energyred_helper`, and `_print_volumepressurered`), so gromos++ `ene_ana` reads the file
//! with the current `ene_ana.md++.lib` (`totene`, `totkin`, `totpot`, `totlj`, `totcrf`, …):
//!
//! * `ENERGY03`: the 52 `# totals` in gromosXX order (total, kinetic, potential, covalent, bond,
//!   angle, improper, dihedral, cross-dihedral, nonbonded, LJ, CRF, lattice-sum slots, …, special,
//!   SASA, constraints, distance restraints, …), then `# baths` (kinetic energy per temperature
//!   bath — the engine tracks the total per bath, so the COM / internal-rotation split is written
//!   as zeros), `# bonded` per energy group (the engine tracks bonded totals, not per group: group
//!   1 carries them), `# nonbonded` per group pair (LJ and CRF from the engine's own per-group
//!   accumulators; lattice-sum slots zero), `# special` per group (zeros — restraint energies are
//!   in the totals), and the empty EDS / GaMD / precalculated-λ / AB-dihedral sections.
//! * `VOLUMEPRESSURE03`: total mass, temperature per bath (total; COM / internal split and the
//!   scaling factor are not tracked: 0, 0, 1), volume, the box vectors, the pressure (the two
//!   scalars gromosXX prints after it — virial and kinetic contributions — are not tracked: 0),
//!   the pressure tensor (not tracked: zeros), the virial tensor and the molecular kinetic-energy
//!   tensor.
//!
//! Numbers use Rust's `e` exponent (`1.5e-2`, not `1.500000000e-02`); every GROMOS reader parses it.
//! [`read_energy_trajectory_native`] reads this layout (and gromosXX's own files); [`EnergyReader`]
//! still reads the one-block `ENERTRJ` layout this crate wrote before 0.0.34.

use crate::IoError;
use gromos_core::configuration::Energy;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Energy component identifiers (GROMOS convention)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyBlock {
    Time = 0,
    Kinetic = 1,
    Potential = 2,
    Total = 3,
    Temperature = 4,
    Volume = 5,
    Pressure = 6,

    // Bonded
    Bond = 7,
    Angle = 8,
    ImproperDihedral = 9,
    ProperDihedral = 10,

    // Nonbonded
    LennardJones = 11,
    ElectrostaticReal = 12,
    ElectrostaticReciprocal = 13,
    ElectrostaticSelf = 14,

    // Constraints
    Shake = 15,
    DistanceRestraint = 16,

    // Kinetic components
    KineticTrans = 17,
    KineticRot = 18,
    KineticInt = 19,
}

/// Energy data for a single timestep
#[derive(Debug, Clone)]
pub struct EnergyFrame {
    pub time: f64,
    pub kinetic: f64,
    pub potential: f64,
    pub total: f64,
    pub temperature: f64,
    pub volume: f64,
    pub pressure: f64,

    // Bonded energies
    pub bond: f64,
    pub angle: f64,
    pub improper: f64,
    pub dihedral: f64,

    // Nonbonded energies
    pub lj: f64,
    pub coul_real: f64,
    pub coul_recip: f64,
    pub coul_self: f64,

    // Constraints
    pub shake: f64,
    pub restraint: f64,

    // Additional components (optional)
    pub extra: Vec<f64>,
    /// MD step number (`TIMESTEP` block).
    pub step: usize,
    /// Special-interaction total (position restraints etc.; `ENERGY03` slot 18).
    pub special: f64,
    /// SASA solvation energy (slot 19).
    pub sasa: f64,
    /// Kinetic energy per temperature bath (`# baths`).
    pub bath_kinetic: Vec<f64>,
    /// LJ energy per energy-group pair, `group_lj[i][j]` (`# nonbonded`).
    pub group_lj: Vec<Vec<f64>>,
    /// CRF energy per energy-group pair, `group_crf[i][j]`.
    pub group_crf: Vec<Vec<f64>>,
    /// Total mass of the system (`VOLUMEPRESSURE03`).
    pub mass: f64,
    /// Box vectors K, L, M as rows.
    pub box_vectors: [[f64; 3]; 3],
    /// Virial tensor, `[a][b]`.
    pub virial_tensor: [[f64; 3]; 3],
    /// Molecular kinetic-energy tensor, `[a][b]`.
    pub kinetic_tensor: [[f64; 3]; 3],
}

impl Default for EnergyFrame {
    fn default() -> Self {
        Self {
            time: 0.0,
            kinetic: 0.0,
            potential: 0.0,
            total: 0.0,
            temperature: 0.0,
            volume: 0.0,
            pressure: 0.0,
            bond: 0.0,
            angle: 0.0,
            improper: 0.0,
            dihedral: 0.0,
            lj: 0.0,
            coul_real: 0.0,
            coul_recip: 0.0,
            coul_self: 0.0,
            shake: 0.0,
            restraint: 0.0,
            extra: Vec::new(),
            step: 0,
            special: 0.0,
            sasa: 0.0,
            bath_kinetic: Vec::new(),
            group_lj: Vec::new(),
            group_crf: Vec::new(),
            mass: 0.0,
            box_vectors: [[0.0; 3]; 3],
            virial_tensor: [[0.0; 3]; 3],
            kinetic_tensor: [[0.0; 3]; 3],
        }
    }
}

impl EnergyFrame {
    /// Create energy frame from components
    pub fn new(time: f64, kinetic: f64, potential: f64, temperature: f64) -> Self {
        Self {
            time,
            kinetic,
            potential,
            total: kinetic + potential,
            temperature,
            ..Default::default()
        }
    }

    /// Update total energy (sum of kinetic and potential)
    pub fn update_total(&mut self) {
        self.total = self.kinetic + self.potential;
    }

    /// Update potential energy (sum of all components)
    pub fn update_potential(&mut self) {
        self.potential = self.bond
            + self.angle
            + self.improper
            + self.dihedral
            + self.lj
            + self.coul_real
            + self.coul_recip
            + self.coul_self
            + self.shake
            + self.restraint;
    }

    /// Build a `.tre`-shaped frame from the runtime [`Energy`] accumulator.
    ///
    /// The single conversion point from `gromos-core`'s `Energy` to the GROMOS
    /// `.tre` block layout — reused by the `md` binary's file writer and by the
    /// `pyo3-gromos` in-memory `Simulation.run()` array, so the two never drift
    /// apart on which component landed in which column (they used to: `angle`,
    /// `dihedral`, and `improper` were separately hardcoded to `0.0` in each).
    pub fn from_energy(
        energy: &Energy,
        time: f64,
        temperature: f64,
        volume: f64,
        pressure: f64,
    ) -> Self {
        Self {
            time,
            kinetic: energy.kinetic_total,
            potential: energy.potential_total,
            total: energy.total(),
            temperature,
            volume,
            pressure,
            bond: energy.bond_total,
            angle: energy.angle_total,
            improper: energy.improper_total,
            dihedral: energy.dihedral_total,
            lj: energy.lj_total,
            coul_real: energy.crf_total,
            coul_recip: energy.ls_kspace_total,
            coul_self: 0.0,
            shake: energy.constraint_total,
            restraint: energy.distanceres_total,
            special: energy.special_total,
            sasa: energy.sasa_total,
            bath_kinetic: energy.kinetic_energy.clone(),
            group_lj: energy.lj_energy.clone(),
            group_crf: energy.crf_energy.clone(),
            ..Default::default()
        }
    }
}

/// Writer for gromosXX-layout energy trajectories (`.tre`).
pub struct EnergyWriter {
    writer: BufWriter<File>,
    frame_count: usize,
    n_baths: usize,
    n_groups: usize,
}

/// Number of entries in an `ENERGY03` `# totals` section (gromosXX `Energy`, 2023-04-15 layout).
const N_TOTALS: usize = 52;

impl EnergyWriter {
    /// Create the file with its `TITLE` and `ENEVERSION` blocks. Per-bath / per-group sections
    /// default to one bath and one energy group; set the run's counts with
    /// [`with_layout`](Self::with_layout).
    pub fn new<P: AsRef<Path>>(path: P, title: &str) -> Result<Self, IoError> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "TITLE")?;
        writeln!(writer, "\t{}", title)?;
        writeln!(writer, "\n\tenergy trajectory")?;
        writeln!(writer, "END")?;
        writeln!(writer, "ENEVERSION")?;
        writeln!(writer, "\t2023-04-15")?;
        writeln!(writer, "END")?;
        Ok(Self {
            writer,
            frame_count: 0,
            n_baths: 1,
            n_groups: 1,
        })
    }

    /// The number of temperature baths and energy groups of the run.
    pub fn with_layout(mut self, n_baths: usize, n_groups: usize) -> Self {
        self.n_baths = n_baths.max(1);
        self.n_groups = n_groups.max(1);
        self
    }

    /// Append one frame: `TIMESTEP`, `ENERGY03` and `VOLUMEPRESSURE03` blocks.
    pub fn write_frame(&mut self, frame: &EnergyFrame) -> Result<(), IoError> {
        let (nb, ng) = (self.n_baths, self.n_groups);
        let w = &mut self.writer;
        let num = |v: f64| format!("{:>18.9e}", v);
        let row = |vals: &[f64]| vals.iter().map(|&v| num(v)).collect::<Vec<_>>().join("");
        writeln!(w, "TIMESTEP")?;
        writeln!(w, "{:>15}{:>15.9}", frame.step, frame.time)?;
        writeln!(w, "END")?;

        writeln!(w, "ENERGY03")?;
        writeln!(w, "# totals")?;
        let mut totals = [0.0f64; N_TOTALS];
        let bonded = frame.bond + frame.angle + frame.improper + frame.dihedral;
        totals[0] = frame.total;
        totals[1] = frame.kinetic;
        totals[2] = frame.potential;
        totals[3] = bonded;
        totals[4] = frame.bond;
        totals[5] = frame.angle;
        totals[6] = frame.improper;
        totals[7] = frame.dihedral;
        totals[9] = frame.lj + frame.coul_real + frame.coul_recip;
        totals[10] = frame.lj;
        totals[11] = frame.coul_real;
        totals[13] = frame.coul_recip;
        totals[18] = frame.special;
        totals[19] = frame.sasa;
        totals[21] = frame.shake;
        totals[22] = frame.restraint;
        for v in totals {
            writeln!(w, "{}", num(v))?;
        }
        writeln!(w, "# baths")?;
        writeln!(w, "{nb}")?;
        for b in 0..nb {
            let kin = frame.bath_kinetic.get(b).copied().unwrap_or(0.0);
            writeln!(w, "{}", row(&[kin, 0.0, 0.0]))?;
        }
        writeln!(w, "# bonded")?;
        writeln!(w, "{ng}")?;
        for g in 0..ng {
            let v = if g == 0 {
                [frame.bond, frame.angle, frame.improper, frame.dihedral, 0.0]
            } else {
                [0.0; 5]
            };
            writeln!(w, "{}", row(&v))?;
        }
        writeln!(w, "# nonbonded")?;
        for i in 0..ng {
            for j in i..ng {
                let lj = frame
                    .group_lj
                    .get(i)
                    .and_then(|r| r.get(j))
                    .copied()
                    .unwrap_or(0.0);
                let crf = frame
                    .group_crf
                    .get(i)
                    .and_then(|r| r.get(j))
                    .copied()
                    .unwrap_or(0.0);
                writeln!(w, "{}", row(&[lj, crf, 0.0, 0.0, 0.0, 0.0]))?;
            }
        }
        writeln!(w, "# special")?;
        for _ in 0..ng {
            writeln!(w, "{}", row(&[0.0; 13]))?;
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
        writeln!(w, "{}", row(&[0.0, 0.0]))?;
        writeln!(w, "END")?;

        writeln!(w, "VOLUMEPRESSURE03")?;
        writeln!(w, "# mass")?;
        writeln!(w, "{}", num(frame.mass))?;
        writeln!(w, "# temperature")?;
        writeln!(w, "{nb}")?;
        for _ in 0..nb {
            writeln!(w, "{}", row(&[frame.temperature, 0.0, 0.0, 1.0]))?;
        }
        writeln!(w, "# volume")?;
        writeln!(w, "{}", num(frame.volume))?;
        for r in &frame.box_vectors {
            writeln!(w, "{}", row(r))?;
        }
        writeln!(w, "# pressure")?;
        writeln!(w, "{}", num(frame.pressure))?;
        writeln!(w, "{}", num(0.0))?;
        writeln!(w, "{}", num(0.0))?;
        for _ in 0..3 {
            writeln!(w, "{}", row(&[0.0; 3]))?;
        }
        for r in &frame.virial_tensor {
            writeln!(w, "{}", row(r))?;
        }
        for r in &frame.kinetic_tensor {
            writeln!(w, "{}", row(r))?;
        }
        writeln!(w, "END")?;
        self.frame_count += 1;
        Ok(())
    }

    /// Kept for callers of the pre-0.0.34 API; the native layout has no trailing block.
    pub fn finalize(&mut self) -> Result<(), IoError> {
        self.flush()
    }

    pub fn flush(&mut self) -> Result<(), IoError> {
        self.writer.flush().map_err(IoError::Io)
    }

    pub fn frame_count(&self) -> usize {
        self.frame_count
    }
}

/// Read a native-layout energy trajectory (this crate's writer since 0.0.34, or gromosXX's own
/// files): the `# totals` of every `ENERGY03` block, the time from the preceding `TIMESTEP`,
/// volume / pressure / temperature from `VOLUMEPRESSURE03`. Comment lines are skipped.
pub fn read_energy_trajectory_native<P: AsRef<Path>>(path: P) -> Result<Vec<EnergyFrame>, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|e| IoError::FileNotFound(format!("{}: {e}", path.as_ref().display())))?;
    let reader = BufReader::new(file);
    #[derive(PartialEq)]
    enum Block {
        None,
        Timestep,
        Energy,
        VolPres,
    }
    let mut frames: Vec<EnergyFrame> = Vec::new();
    let mut block = Block::None;
    let mut section = String::new();
    let mut step = 0usize;
    let mut time = 0.0;
    let mut totals: Vec<f64> = Vec::new();
    let mut vp: Vec<(String, Vec<f64>)> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(IoError::Io)?;
        let t = line.trim();
        match t {
            "TIMESTEP" => {
                block = Block::Timestep;
                continue;
            },
            "ENERGY03" => {
                block = Block::Energy;
                section.clear();
                totals.clear();
                continue;
            },
            "VOLUMEPRESSURE03" => {
                block = Block::VolPres;
                section.clear();
                vp.clear();
                continue;
            },
            "END" => {
                match block {
                    Block::Energy => {
                        let at = |i: usize| totals.get(i).copied().unwrap_or(0.0);
                        frames.push(EnergyFrame {
                            step,
                            time,
                            total: at(0),
                            kinetic: at(1),
                            potential: at(2),
                            bond: at(4),
                            angle: at(5),
                            improper: at(6),
                            dihedral: at(7),
                            lj: at(10),
                            coul_real: at(11),
                            coul_recip: at(13),
                            special: at(18),
                            sasa: at(19),
                            shake: at(21),
                            restraint: at(22),
                            ..Default::default()
                        });
                    },
                    Block::VolPres => {
                        if let Some(f) = frames.last_mut() {
                            for (name, vals) in &vp {
                                match (name.as_str(), vals.first()) {
                                    ("volume", Some(v)) => f.volume = *v,
                                    ("pressure", Some(v)) => f.pressure = *v,
                                    ("mass", Some(v)) => f.mass = *v,
                                    ("temperature", _)
                                        // count, then one line of 4 per bath: the first bath's total
                                        if vals.len() >= 2 => {
                                            f.temperature = vals[1];
                                        },
                                    _ => {},
                                }
                            }
                        }
                    },
                    _ => {},
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
                    .filter_map(|x| x.parse().ok())
                    .collect();
                if vals.len() >= 2 {
                    step = vals[0] as usize;
                    time = vals[1];
                }
            },
            Block::Energy => {
                if let Some(name) = t.strip_prefix('#') {
                    section = name.trim().to_string();
                    continue;
                }
                if section.is_empty() || section == "totals" {
                    totals.extend(t.split_whitespace().filter_map(|x| x.parse::<f64>().ok()));
                }
            },
            Block::VolPres => {
                if let Some(name) = t.strip_prefix('#') {
                    vp.push((name.trim().to_string(), Vec::new()));
                    continue;
                }
                if let Some((_, vals)) = vp.last_mut() {
                    vals.extend(t.split_whitespace().filter_map(|x| x.parse::<f64>().ok()));
                }
            },
            Block::None => {},
        }
    }
    Ok(frames)
}

/// Write a single energy frame to a file (convenience function)
pub fn write_energy_frame<P: AsRef<Path>>(
    path: P,
    frame: &EnergyFrame,
    title: &str,
) -> Result<(), IoError> {
    let mut writer = EnergyWriter::new(path, title)?;
    writer.write_frame(frame)?;
    writer.finalize()?;
    Ok(())
}

/// GROMOS energy file reader
pub struct EnergyReader {
    reader: BufReader<File>,
    title: String,
    frames_read: usize,
}

impl EnergyReader {
    /// Open an energy trajectory file for reading
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, IoError> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read TITLE block
        let title = Self::read_title_block(&mut reader)?;

        Ok(Self {
            reader,
            title,
            frames_read: 0,
        })
    }

    /// Get the title from the energy file
    pub fn title(&self) -> &str {
        &self.title
    }

    /// Get the number of frames read so far
    pub fn frames_read(&self) -> usize {
        self.frames_read
    }

    /// Read a single energy frame
    ///
    /// Returns Ok(Some(frame)) if a frame was read, Ok(None) if end of file reached
    pub fn read_frame(&mut self) -> Result<Option<EnergyFrame>, IoError> {
        let mut line = String::new();

        loop {
            line.clear();
            let bytes_read = self.reader.read_line(&mut line)?;

            if bytes_read == 0 {
                // End of file
                return Ok(None);
            }

            let trimmed = line.trim();

            // Skip empty lines and comments
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }

            // Check for END block
            if trimmed == "END" {
                return Ok(None);
            }

            // Check for ENERTRJ block
            if trimmed == "ENERTRJ" {
                continue;
            }

            // Parse energy values
            let values: Result<Vec<f64>, _> = trimmed
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect();

            match values {
                Ok(vals) if vals.len() >= 17 => {
                    // Successfully parsed energy frame
                    let frame = EnergyFrame {
                        time: vals[0],
                        kinetic: vals[1],
                        potential: vals[2],
                        total: vals[3],
                        temperature: vals[4],
                        volume: vals[5],
                        pressure: vals[6],
                        bond: vals[7],
                        angle: vals[8],
                        improper: vals[9],
                        dihedral: vals[10],
                        lj: vals[11],
                        coul_real: vals[12],
                        coul_recip: vals[13],
                        coul_self: vals[14],
                        shake: vals[15],
                        restraint: vals[16],
                        extra: if vals.len() > 17 {
                            vals[17..].to_vec()
                        } else {
                            Vec::new()
                        },
                        ..Default::default()
                    };

                    self.frames_read += 1;
                    return Ok(Some(frame));
                },
                Ok(vals) if vals.is_empty() => {
                    // Empty line after splitting, continue
                    continue;
                },
                Ok(vals) => {
                    return Err(IoError::ParseError(format!(
                        "Incomplete energy frame: expected at least 17 values, got {}",
                        vals.len()
                    )));
                },
                Err(e) => {
                    return Err(IoError::ParseError(format!(
                        "Failed to parse energy values: {}",
                        e
                    )));
                },
            }
        }
    }

    /// Read all energy frames from the file
    pub fn read_all_frames(&mut self) -> Result<Vec<EnergyFrame>, IoError> {
        let mut frames = Vec::new();
        while let Some(frame) = self.read_frame()? {
            frames.push(frame);
        }
        Ok(frames)
    }

    /// Read TITLE block
    fn read_title_block(reader: &mut BufReader<File>) -> Result<String, IoError> {
        let mut line = String::new();
        let mut title_lines = Vec::new();
        let mut in_title = false;

        loop {
            line.clear();
            let bytes_read = reader.read_line(&mut line)?;

            if bytes_read == 0 {
                return Err(IoError::ParseError(
                    "Unexpected EOF while reading TITLE".to_string(),
                ));
            }

            let trimmed = line.trim();

            if trimmed == "TITLE" {
                in_title = true;
                continue;
            }

            if trimmed == "END" {
                if in_title {
                    return Ok(title_lines.join(" ").trim().to_string());
                }
                continue;
            }

            if in_title {
                title_lines.push(trimmed.to_string());
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_energy_frame_creation() {
        let frame = EnergyFrame::new(0.0, 100.0, -200.0, 300.0);
        assert_eq!(frame.time, 0.0);
        assert_eq!(frame.kinetic, 100.0);
        assert_eq!(frame.potential, -200.0);
        assert_eq!(frame.total, -100.0);
        assert_eq!(frame.temperature, 300.0);
    }

    #[test]
    fn test_energy_writer_creation() {
        let temp_file = "/tmp/test_energy.tre";
        let writer = EnergyWriter::new(temp_file, "Test energy");
        assert!(writer.is_ok());
    }

    #[test]
    fn test_write_energy_frame() {
        let temp_file = "/tmp/test_energy2.tre";
        let mut writer = EnergyWriter::new(temp_file, "Test").unwrap();

        let mut frame = EnergyFrame::new(0.0, 150.0, -300.0, 298.15);
        frame.bond = -50.0;
        frame.angle = -30.0;
        frame.lj = -150.0;
        frame.coul_real = -70.0;
        frame.update_potential();

        let result = writer.write_frame(&frame);
        assert!(result.is_ok());
        assert_eq!(writer.frame_count(), 1);

        let result = writer.finalize();
        assert!(result.is_ok());
    }
}
