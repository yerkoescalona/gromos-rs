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
pub mod dlg;
pub mod energy;
pub mod energy_binary;
pub mod force;
pub mod g96;
pub mod imd;
pub mod input;
pub mod output;
pub mod pdb;
pub mod ptp;
pub mod topology;
pub mod trajectory;
pub mod trajectory_binary;

// Re-export commonly used types
pub use dlg::{DlgWriter, LambdaDerivativeFrame};
pub use energy::{EnergyBlock, EnergyFrame, EnergyReader, EnergyWriter};
pub use energy_binary::{BinaryEnergyReader, BinaryEnergyWriter};
pub use force::ForceWriter;
pub use imd::{ImdParameters, PressureParameters, TempBathParameters};
pub use input::{EdsBlock, GamdBlock, ReplicaBlock};
pub use output::{EdsStatsWriter, EdsVrWriter, GamdBoostWriter, GamdStatsWriter};
pub use ptp::PtpWriter;
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
