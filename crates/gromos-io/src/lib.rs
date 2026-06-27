//! GROMOS-RS I/O: File format readers and writers
//!
//! This crate handles all file I/O for GROMOS-RS:
//! - .topo/.top - Topology files
//! - .conf/.cnf - Coordinate files  
//! - .imd - Input parameter files (simulation settings)
//! - .trc/.trj - Trajectory files (coordinates over time)
//! - .tre - Energy files (ENE/ENA format)
//! - .trf - Force files
//! - .ptp - Perturbation topology files (FEP)

pub mod coordinate;
pub mod distanceres;
pub mod dlg;
pub mod energy;
pub mod free_energy;
pub mod ifp;
pub mod jobs;
pub mod mk_script;
pub mod mtb;
pub mod energy_binary;
pub mod force;
pub mod g96;
pub mod imd;
pub mod input;
pub mod output;
pub mod pdb;
pub mod posres;
pub mod ptp;
pub mod script_template;
pub mod topology;
pub mod trajectory;
pub mod trajectory_binary;

// Re-export commonly used types
pub use coordinate::{CoordinateData, G96Atom, LabeledCoordinateData, read_coordinates, read_g96_labeled};
pub use dlg::{DlgWriter, LambdaDerivativeFrame};
pub use free_energy::{FreeEnergyWriter, FreeEnergyFrame, read_free_energy_trajectory};
pub use energy::{EnergyBlock, EnergyFrame, EnergyReader, EnergyWriter};
pub use energy_binary::{BinaryEnergyReader, BinaryEnergyWriter};
pub use force::ForceWriter;
pub use g96::{write_por, write_rpr};
pub use imd::{ImdParameters, PressureParameters, TempBathParameters};
pub use input::{EdsBlock, GamdBlock, ReplicaBlock};
pub use output::{EdsStatsWriter, EdsVrWriter, GamdBoostWriter, GamdStatsWriter};
pub use ptp::{PtpWriter, PerturbedTopology, read_pttopo};
pub use pdb::write_pdb_positions;
pub use trajectory::TrajectoryWriter;
pub use trajectory_binary::{
    BinaryFrame, BinaryTrajectoryReader, BinaryTrajectoryWriter, DcdReader, DcdWriter,
};

use std::io;

/// Common error type for I/O operations
#[derive(Debug)]
pub enum IoError {
    FileNotFound(String),
    ParseError(String),
    FormatError(String),
    WriteError(String),
    Io(io::Error),
}

impl From<io::Error> for IoError {
    fn from(err: io::Error) -> Self {
        IoError::Io(err)
    }
}

impl std::fmt::Display for IoError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            IoError::FileNotFound(path) => write!(f, "File not found: {}", path),
            IoError::ParseError(msg) => write!(f, "Parse error: {}", msg),
            IoError::FormatError(msg) => write!(f, "Format error: {}", msg),
            IoError::WriteError(msg) => write!(f, "Write error: {}", msg),
            IoError::Io(err) => write!(f, "I/O error: {}", err),
        }
    }
}

impl std::error::Error for IoError {}

/// Pre-process argv to support GROMOS `@key` conventions.
///
/// Converts `@key` → `--key` and expands `@f argfile` contents inline.
/// Blank lines and `#` comments in argfiles are skipped.
/// This allows using clap with the standard GROMOS argument syntax.
pub fn gromos_args() -> Vec<String> {
    let raw: Vec<String> = std::env::args().collect();
    let mut expanded = vec![raw[0].clone()];

    let mut i = 1;
    while i < raw.len() {
        if raw[i] == "@f" {
            i += 1;
            if i >= raw.len() {
                eprintln!("Error: @f requires a filename");
                std::process::exit(1);
            }
            match std::fs::read_to_string(&raw[i]) {
                Ok(content) => {
                    for line in content.lines() {
                        let line = line.trim();
                        if line.is_empty() || line.starts_with('#') {
                            continue;
                        }
                        for token in line.split_whitespace() {
                            if let Some(key) = token.strip_prefix('@') {
                                expanded.push(format!("--{}", key));
                            } else {
                                expanded.push(token.to_string());
                            }
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Error: cannot read argfile '{}': {}", raw[i], e);
                    std::process::exit(1);
                }
            }
        } else if let Some(key) = raw[i].strip_prefix('@') {
            expanded.push(format!("--{}", key));
        } else {
            expanded.push(raw[i].clone());
        }
        i += 1;
    }
    expanded
}
