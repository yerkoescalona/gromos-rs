//! Parser for position restraint specification (.por) and reference position (.rpr) files.
//!
//! gromosXX convention:
//! - POSRESSPEC block (.por): specifies which atoms are restrained (atom indices only)
//!   First 17 characters per line are ignored, then atom number (1-based) is read.
//! - REFPOSITION block (.rpr or in .conf): reference positions for all atoms
//!   Same format as POSITION block.
//!
//! When NTPORB=0, reference positions come from the startup configuration.
//! When NTPORB=1, reference positions come from the @refpos file.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use gromos_core::math::Vec3;

use crate::IoError;

/// A single position restraint specification entry.
#[derive(Debug, Clone)]
pub struct PosResEntry {
    /// Atom index (0-based)
    pub atom: usize,
    /// Reference position (nm)
    pub reference_pos: Vec3,
}

/// Parse a POSRESSPEC block from a .por file.
///
/// Returns 0-based atom indices of restrained atoms.
/// gromosXX convention: first 17 characters per line are ignored,
/// then atom number (1-based) is read.
pub fn read_posresspec<P: AsRef<Path>>(path: P) -> Result<Vec<usize>, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut atoms = Vec::new();
    let mut in_block = false;

    for line in reader.lines() {
        let line = line.map_err(|e| IoError::ParseError(e.to_string()))?;
        let trimmed = line.trim();

        if trimmed == "POSRESSPEC" {
            in_block = true;
            continue;
        }
        if trimmed == "END" {
            if in_block {
                break;
            }
            continue;
        }
        if !in_block || trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // gromosXX: first 17 chars ignored, then read atom number
        // Our format also has: res_num res_name atom_name atom_index [x y z]
        // Parse atom index from the 4th token (both formats work)
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() >= 4 {
            let atom_index: usize = tokens[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom index: {}", tokens[3])))?;
            atoms.push(atom_index - 1); // Convert to 0-based
        }
    }

    Ok(atoms)
}

/// Parse a REFPOSITION block from a .rpr file.
///
/// Returns reference positions for ALL atoms in the system, indexed by atom number.
/// Same format as POSITION block in coordinate files.
pub fn read_refpos<P: AsRef<Path>>(path: P) -> Result<Vec<Vec3>, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut positions = Vec::new();
    let mut in_block = false;

    for line in reader.lines() {
        let line = line.map_err(|e| IoError::ParseError(e.to_string()))?;
        let trimmed = line.trim();

        if trimmed == "REFPOSITION" {
            in_block = true;
            continue;
        }
        if trimmed == "END" {
            if in_block {
                break;
            }
            continue;
        }
        if !in_block || trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        // Same format as POSITION: first 24 chars ignored, then x y z
        let tokens: Vec<&str> = trimmed.split_whitespace().collect();
        if tokens.len() >= 7 {
            // res_num res_name atom_name atom_index x y z
            let x: f64 = tokens[4]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid x: {}", tokens[4])))?;
            let y: f64 = tokens[5]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid y: {}", tokens[5])))?;
            let z: f64 = tokens[6]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid z: {}", tokens[6])))?;
            positions.push(Vec3::new(x, y, z));
        }
    }

    Ok(positions)
}

/// Build PosResEntry list by combining atom indices from POSRESSPEC with
/// reference positions from REFPOSITION (or startup configuration).
pub fn build_posres_entries(
    restrained_atoms: &[usize],
    reference_positions: &[Vec3],
) -> Vec<PosResEntry> {
    restrained_atoms
        .iter()
        .map(|&atom| PosResEntry {
            atom,
            reference_pos: reference_positions[atom],
        })
        .collect()
}
