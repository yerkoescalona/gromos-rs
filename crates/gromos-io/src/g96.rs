//! GROMOS96 (.g96) file format reader/writer
//!
//! GROMOS96 format is a block-based format with keywords like TITLE, POSITION, VELOCITY, BOX

use gromos_core::configuration::Configuration;
use gromos_core::math::Vec3;
use gromos_core::topology::Topology;
use crate::coordinate::G96Atom;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Write a GROMOS96 format file
pub struct G96Writer {
    pub writer: BufWriter<File>,
}

impl G96Writer {
    /// Create a new G96 writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self, String> {
        let file = File::create(path).map_err(|e| format!("Cannot create .g96 file: {}", e))?;
        Ok(G96Writer {
            writer: BufWriter::new(file),
        })
    }

    /// Write title block
    pub fn write_title(&mut self, title: &str) -> Result<(), String> {
        writeln!(self.writer, "TITLE").map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "{}", title).map_err(|e| format!("Write error: {}", e))?;
        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write position block
    ///
    /// GROMOS96 position format:
    /// ```text
    /// POSITION
    /// %5d %-5s %-5s%7d%15.9f%15.9f%15.9f
    /// (resnum, resname, atomname, atomnum, x, y, z)
    /// END
    /// ```
    pub fn write_positions(
        &mut self,
        positions: &[Vec3],
        topology: Option<&Topology>,
    ) -> Result<(), String> {
        writeln!(self.writer, "POSITION").map_err(|e| format!("Write error: {}", e))?;

        for (i, pos) in positions.iter().enumerate() {
            let (resnum, resname, atomname) = if let Some(topo) = topology {
                let n_solute = topo.solute.num_atoms();
                if i < n_solute {
                    let atom = &topo.solute.atoms[i];
                    (atom.residue_nr + 1, atom.residue_name.as_str(), atom.name.as_str())
                } else if !topo.solvents.is_empty() {
                    let sv = &topo.solvents[0];
                    let atoms_per_mol = sv.atoms_per_molecule();
                    let sv_idx = i - n_solute;
                    let mol_nr = sv_idx / atoms_per_mol;
                    let atom_in_mol = sv_idx % atoms_per_mol;
                    let sv_atom = &sv.atoms[atom_in_mol];
                    (n_solute / atoms_per_mol + mol_nr + 1, sv_atom.residue_name.as_str(), sv_atom.name.as_str())
                } else {
                    (1, "UNK", "AT")
                }
            } else {
                (1, "UNK", "AT")
            };
            let atomnum = i + 1;

            writeln!(
                self.writer,
                "{:>5} {:<5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
                resnum, resname, atomname, atomnum, pos.x, pos.y, pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write velocity block (optional)
    pub fn write_velocities(&mut self, velocities: &[Vec3]) -> Result<(), String> {
        writeln!(self.writer, "VELOCITY").map_err(|e| format!("Write error: {}", e))?;

        for vel in velocities {
            let vx = vel.x;
            let vy = vel.y;
            let vz = vel.z;

            writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", vx, vy, vz)
                .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write box dimensions (for periodic boundary conditions)
    pub fn write_box(&mut self, box_dims: Vec3) -> Result<(), String> {
        writeln!(self.writer, "BOX").map_err(|e| format!("Write error: {}", e))?;

        let x = box_dims.x;
        let y = box_dims.y;
        let z = box_dims.z;

        writeln!(self.writer, "{:15.9}{:15.9}{:15.9}", x, y, z)
            .map_err(|e| format!("Write error: {}", e))?;

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write position block with full atom labels (residue number, name, atom name).
    ///
    /// Atoms are renumbered sequentially starting from 1.
    pub fn write_labeled_positions(&mut self, atoms: &[G96Atom]) -> Result<(), String> {
        writeln!(self.writer, "POSITION").map_err(|e| format!("Write error: {}", e))?;

        for (i, atom) in atoms.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
                atom.res_num, atom.res_name, atom.atom_name, i + 1,
                atom.pos.x, atom.pos.y, atom.pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write POSRESSPEC block — specifies which atoms are positionally restrained.
    ///
    /// Same format as POSITION but with block name POSRESSPEC.
    /// Only the atom sequence numbers matter to the MD engine (first 17 chars
    /// and coordinates are ignored), but we write full labels for readability.
    pub fn write_posresspec(&mut self, atoms: &[G96Atom]) -> Result<(), String> {
        writeln!(self.writer, "POSRESSPEC").map_err(|e| format!("Write error: {}", e))?;

        for (i, atom) in atoms.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
                atom.res_num, atom.res_name, atom.atom_name, i + 1,
                atom.pos.x, atom.pos.y, atom.pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Write REFPOSITION block — reference positions for position restraints.
    ///
    /// Same format as POSITION but with block name REFPOSITION.
    /// Must contain all atoms; the MD engine uses coordinates of atoms
    /// listed in the POSRESSPEC file as restraint targets.
    pub fn write_refposition(&mut self, atoms: &[G96Atom]) -> Result<(), String> {
        writeln!(self.writer, "REFPOSITION").map_err(|e| format!("Write error: {}", e))?;

        for (i, atom) in atoms.iter().enumerate() {
            writeln!(
                self.writer,
                "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
                atom.res_num, atom.res_name, atom.atom_name, i + 1,
                atom.pos.x, atom.pos.y, atom.pos.z
            )
            .map_err(|e| format!("Write error: {}", e))?;
        }

        writeln!(self.writer, "END").map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Flush and close the writer
    pub fn close(mut self) -> Result<(), String> {
        self.writer
            .flush()
            .map_err(|e| format!("Flush error: {}", e))
    }
}

/// Write a configuration to a GROMOS96 file
pub fn write_g96<P: AsRef<Path>>(
    path: P,
    title: &str,
    positions: &[Vec3],
    velocities: Option<&[Vec3]>,
    box_dims: Option<Vec3>,
    topology: Option<&Topology>,
) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;

    writer.write_title(title)?;
    writer.write_positions(positions, topology)?;

    if let Some(vels) = velocities {
        writer.write_velocities(vels)?;
    }

    if let Some(box_size) = box_dims {
        writer.write_box(box_size)?;
    }

    writer.close()
}

/// Generate a .por (position restraint specification) file from coordinate data.
///
/// The .por file contains a POSRESSPEC block with only the solute atoms
/// (first `n_solute_atoms` atoms). The MD engine uses only the atom indices
/// from this block to determine which atoms are restrained.
pub fn write_por<P: AsRef<Path>>(
    path: P,
    atoms: &[G96Atom],
    n_solute_atoms: usize,
) -> Result<(), String> {
    let solute = &atoms[..n_solute_atoms];
    let mut writer = G96Writer::new(path)?;
    writer.write_title("Position restraint specification: solute atoms")?;
    writer.write_posresspec(solute)?;
    writer.close()
}

/// Generate a .rpr (reference positions) file from coordinate data.
///
/// The .rpr file contains a REFPOSITION block with all atoms.
/// The MD engine reads reference coordinates for restrained atoms from this file.
pub fn write_rpr<P: AsRef<Path>>(
    path: P,
    atoms: &[G96Atom],
) -> Result<(), String> {
    let mut writer = G96Writer::new(path)?;
    writer.write_title("Reference positions for position restraints")?;
    writer.write_refposition(atoms)?;
    writer.close()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_g96_write() {
        let positions = vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 5.0, 6.0)];

        let result = write_g96(
            "/tmp/test.g96",
            "Test structure",
            &positions,
            None,
            Some(Vec3::new(10.0, 10.0, 10.0)),
            None,
        );

        assert!(result.is_ok());
    }
}
