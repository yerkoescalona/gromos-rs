//! Parser for GROMOS topology files (.topo/.top)
//!
//! Format blocks:
//! - SOLUTEATOM: atom definitions (mass, charge, IAC)
//! - BOND/BONDSTRETCHTYPE: covalent bonds
//! - BONDANGLE/BONDANGLEBENDTYPE: bond angles
//! - DIHEDRAL/TORSDIHEDRALTYPE: torsional angles
//! - CGPARAMETERS: LJ interaction parameters
//! - TEMPERATUREGROUPS: temperature coupling groups

use crate::IoError;
use gromos_core::topology::{Angle, AngleParameters, Bond, BondParameters, Dihedral, DihedralParameters, ImproperDihedralParameters, LJParameters, Topology};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Parsed solvent atom template
#[derive(Debug, Clone)]
pub struct ParsedSolventAtom {
    pub iac: usize,
    pub name: String,
    pub mass: f64,
    pub charge: f64,
}

/// Parsed solvent constraint
#[derive(Debug, Clone)]
pub struct ParsedSolventConstraint {
    pub i: usize, // 0-indexed within solvent molecule
    pub j: usize,
    pub length: f64,
}

/// Parsed topology data from GROMOS .topo file
#[derive(Debug)]
pub struct ParsedTopology {
    pub n_atoms: usize,
    pub masses: Vec<f64>,
    pub charges: Vec<f64>,
    pub iac: Vec<usize>,
    pub chargegroup_codes: Vec<usize>, // CGC values: 1 = end of chargegroup
    pub exclusions: Vec<Vec<usize>>,
    pub bonds: Vec<(usize, usize, usize)>, // (i, j, type)
    pub bond_parameters: Vec<BondParameters>,
    pub angles: Vec<(usize, usize, usize, usize)>, // (i, j, k, type)
    pub angle_parameters: Vec<AngleParameters>,
    pub proper_dihedrals: Vec<(usize, usize, usize, usize, usize)>, // (i, j, k, l, type)
    pub dihedral_parameters: Vec<DihedralParameters>,
    pub improper_dihedrals: Vec<(usize, usize, usize, usize, usize)>, // (i, j, k, l, type)
    pub improper_dihedral_parameters: Vec<ImproperDihedralParameters>,
    pub lj_parameters: HashMap<(usize, usize), LJParameters>,
    pub temperature_groups: Vec<usize>, // Last atom index of each group
    pub solvent_atoms: Vec<ParsedSolventAtom>,
    pub solvent_constraints: Vec<ParsedSolventConstraint>,
}

/// Read GROMOS topology file
pub fn read_topology_file<P: AsRef<Path>>(path: P) -> Result<ParsedTopology, IoError> {
    let file = File::open(path.as_ref())
        .map_err(|_| IoError::FileNotFound(path.as_ref().display().to_string()))?;
    let reader = BufReader::new(file);

    let mut lines = reader.lines();

    let mut n_atoms = 0;
    let mut masses = Vec::new();
    let mut charges = Vec::new();
    let mut iac = Vec::new();
    let mut chargegroup_codes = Vec::new();
    let mut exclusions = Vec::new();
    let mut bonds = Vec::new();
    let mut bond_parameters = Vec::new();
    let mut angles = Vec::new();
    let mut angle_parameters = Vec::new();
    let mut lj_parameters = HashMap::new();
    let mut temperature_groups = Vec::new();
    let mut proper_dihedrals = Vec::new();
    let mut dihedral_parameters = Vec::new();
    let mut improper_dihedrals = Vec::new();
    let mut improper_dihedral_parameters = Vec::new();
    let mut solvent_atoms = Vec::new();
    let mut solvent_constraints = Vec::new();

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        // Skip comments and empty lines
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        match trimmed {
            "SOLUTEATOM" => {
                parse_soluteatom(
                    &mut lines,
                    &mut n_atoms,
                    &mut masses,
                    &mut charges,
                    &mut iac,
                    &mut chargegroup_codes,
                    &mut exclusions,
                )?;
            },
            "BONDSTRETCHTYPE" => {
                parse_bond_types(&mut lines, &mut bond_parameters)?;
            },
            "BOND" | "BONDH" => {
                parse_bonds(&mut lines, &mut bonds)?;
            },
            "BONDANGLEBENDTYPE" => {
                parse_angle_types(&mut lines, &mut angle_parameters)?;
            },
            "BONDANGLE" | "BONDANGLEH" => {
                parse_angles(&mut lines, &mut angles)?;
            },
            "TORSDIHEDRALTYPE" => {
                parse_dihedral_types(&mut lines, &mut dihedral_parameters)?;
            },
            "DIHEDRAL" | "DIHEDRALH" => {
                parse_dihedrals(&mut lines, &mut proper_dihedrals)?;
            },
            "IMPDIHEDRALTYPE" => {
                parse_improper_dihedral_types(&mut lines, &mut improper_dihedral_parameters)?;
            },
            "IMPDIHEDRAL" | "IMPDIHEDRALH" => {
                parse_dihedrals(&mut lines, &mut improper_dihedrals)?;
            },
            "CGPARAMETERS" | "LJPARAMETERS" => {
                parse_lj_parameters(&mut lines, &mut lj_parameters)?;
            },
            "TEMPERATUREGROUPS" => {
                parse_temperature_groups(&mut lines, &mut temperature_groups)?;
            },
            "SOLVENTATOM" => {
                parse_solventatom(&mut lines, &mut solvent_atoms)?;
            },
            "SOLVENTCONSTR" => {
                parse_solventconstr(&mut lines, &mut solvent_constraints)?;
            },
            _ => {
                // Skip unknown blocks
                continue;
            },
        }
    }

    Ok(ParsedTopology {
        n_atoms,
        masses,
        charges,
        iac,
        chargegroup_codes,
        exclusions,
        bonds,
        bond_parameters,
        angles,
        angle_parameters,
        proper_dihedrals,
        dihedral_parameters,
        improper_dihedrals,
        improper_dihedral_parameters,
        lj_parameters,
        temperature_groups,
        solvent_atoms,
        solvent_constraints,
    })
}

fn parse_soluteatom<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    n_atoms: &mut usize,
    masses: &mut Vec<f64>,
    charges: &mut Vec<f64>,
    iac: &mut Vec<usize>,
    chargegroup_codes: &mut Vec<usize>,
    exclusions: &mut Vec<Vec<usize>>,
) -> Result<(), IoError> {
    // First line after SOLUTEATOM is atom count
    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First non-comment line is atom count
        if *n_atoms == 0 {
            *n_atoms = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid atom count: {}", trimmed)))?;
            continue;
        }

        // Parse atom line (format: ATNM MRES PANM IAC MASS CG CGC INE ...)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 8 {
            let atom_iac: usize = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid IAC: {}", parts[3])))?;
            let mass: f64 = parts[4]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid mass: {}", parts[4])))?;
            let charge: f64 = parts[5]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid charge: {}", parts[5])))?;
            let cgc: usize = parts[6]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CGC: {}", parts[6])))?;
            let n_exclusions: usize = parts[7].parse().map_err(|_| {
                IoError::ParseError(format!("Invalid exclusion count: {}", parts[7]))
            })?;

            iac.push(atom_iac);
            masses.push(mass);
            charges.push(charge);
            chargegroup_codes.push(cgc);

            // Exclusion indices are on the SAME line, starting at parts[8]
            // Format: ATNM MRES PANM IAC MASS CG CGC INE [excl1 excl2 ...]
            let mut atom_exclusions = Vec::new();
            for idx in 8..8 + n_exclusions {
                if idx < parts.len() {
                    if let Ok(excl_idx) = parts[idx].parse::<usize>() {
                        // Convert from 1-based (GROMOS topology) to 0-based
                        atom_exclusions.push(excl_idx - 1);
                    }
                }
            }
            exclusions.push(atom_exclusions);

            // Skip the INE14 line (1-4 pairs) — always present after each atom
            // Format: INE14 [pair1 pair2 ...]
            // We consume it but don't use it (1-4 pairs are built from dihedrals)
            while let Some(Ok(next_line)) = lines.next() {
                let next_trimmed = next_line.trim();
                if next_trimmed.is_empty() || next_trimmed.starts_with('#') {
                    continue;
                }
                // This is the INE14 line, consumed and discarded
                break;
            }
        }
    }

    Ok(())
}

fn parse_bond_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    bond_parameters: &mut Vec<BondParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is type count
        if n_types == 0 {
            n_types = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid bond type count: {}", trimmed))
            })?;
            continue;
        }

        // Parse bond type (format: CB CHB B0)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let k_quartic: f64 = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CB: {}", parts[0])))?;
            let k_harmonic: f64 = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CHB: {}", parts[1])))?;
            let r0: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid B0: {}", parts[2])))?;

            bond_parameters.push(BondParameters {
                k_quartic,
                k_harmonic,
                r0,
            });
        }
    }

    Ok(())
}

fn parse_bonds<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    bonds: &mut Vec<(usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_bonds = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is bond count
        if n_bonds == 0 {
            n_bonds = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond count: {}", trimmed)))?;
            continue;
        }

        // Parse bond (format: IB JB ICB)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let i: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond atom i: {}", parts[0])))?;
            let j: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond atom j: {}", parts[1])))?;
            let bond_type: usize = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid bond type: {}", parts[2])))?;

            // Convert from 1-indexed to 0-indexed
            bonds.push((i - 1, j - 1, bond_type - 1));
        }
    }

    Ok(())
}

fn parse_angle_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    angle_parameters: &mut Vec<AngleParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is type count
        if n_types == 0 {
            n_types = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid angle type count: {}", trimmed))
            })?;
            continue;
        }

        // Parse angle type (format: CT CHT T0)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 3 {
            let k_cosine: f64 = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CT: {}", parts[0])))?;
            let k_harmonic: f64 = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CHT: {}", parts[1])))?;
            let theta0_deg: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid T0: {}", parts[2])))?;

            angle_parameters.push(AngleParameters {
                k_cosine,
                k_harmonic,
                theta0: theta0_deg.to_radians(),
            });
        }
    }

    Ok(())
}

fn parse_angles<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    angles: &mut Vec<(usize, usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_angles = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is angle count
        if n_angles == 0 {
            n_angles = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle count: {}", trimmed)))?;
            continue;
        }

        // Parse angle (format: IT JT KT ICT)
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 4 {
            let i: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom i: {}", parts[0])))?;
            let j: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom j: {}", parts[1])))?;
            let k: usize = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle atom k: {}", parts[2])))?;
            let angle_type: usize = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid angle type: {}", parts[3])))?;

            // Convert from 1-indexed to 0-indexed
            angles.push((i - 1, j - 1, k - 1, angle_type - 1));
        }
    }

    Ok(())
}

fn parse_dihedrals<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    dihedrals: &mut Vec<(usize, usize, usize, usize, usize)>,
) -> Result<(), IoError> {
    let mut n_dihedrals = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        if trimmed == "END" { break; }

        if n_dihedrals == 0 {
            n_dihedrals = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid dihedral count: {}", trimmed)))?;
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 5 {
            let i: usize = parts[0].parse().map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[0])))?;
            let j: usize = parts[1].parse().map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[1])))?;
            let k: usize = parts[2].parse().map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[2])))?;
            let l: usize = parts[3].parse().map_err(|_| IoError::ParseError(format!("Invalid atom: {}", parts[3])))?;
            let dtype: usize = parts[4].parse().map_err(|_| IoError::ParseError(format!("Invalid type: {}", parts[4])))?;
            dihedrals.push((i - 1, j - 1, k - 1, l - 1, dtype - 1));
        }
    }
    Ok(())
}

fn parse_dihedral_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    params: &mut Vec<DihedralParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        if trimmed == "END" { break; }

        if n_types == 0 {
            n_types = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid dihedral type count: {}", trimmed)))?;
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 3 {
            let cp: f64 = parts[0].parse().map_err(|_| IoError::ParseError(format!("Invalid CP: {}", parts[0])))?;
            let pd: f64 = parts[1].parse().map_err(|_| IoError::ParseError(format!("Invalid PD: {}", parts[1])))?;
            let np: i32 = parts[2].parse().map_err(|_| IoError::ParseError(format!("Invalid NP: {}", parts[2])))?;
            params.push(DihedralParameters {
                k: cp,
                pd: pd * std::f64::consts::PI / 180.0,
                cospd: (pd * std::f64::consts::PI / 180.0).cos(),
                m: np,
            });
        }
    }
    Ok(())
}

fn parse_improper_dihedral_types<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    params: &mut Vec<ImproperDihedralParameters>,
) -> Result<(), IoError> {
    let mut n_types = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        if trimmed == "END" { break; }

        if n_types == 0 {
            n_types = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid improper type count: {}", trimmed)))?;
            continue;
        }

        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 2 {
            let cq: f64 = parts[0].parse().map_err(|_| IoError::ParseError(format!("Invalid CQ: {}", parts[0])))?;
            let q0: f64 = parts[1].parse().map_err(|_| IoError::ParseError(format!("Invalid Q0: {}", parts[1])))?;
            params.push(ImproperDihedralParameters {
                q0: q0 * std::f64::consts::PI / 180.0, // Convert to radians
                k: cq,
            });
        }
    }
    Ok(())
}

fn parse_lj_parameters<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    lj_params: &mut HashMap<(usize, usize), LJParameters>,
) -> Result<(), IoError> {
    let mut n_pairs = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is pair count
        if n_pairs == 0 {
            n_pairs = trimmed
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid LJ pair count: {}", trimmed)))?;
            continue;
        }

        // Parse LJ parameters (format: IAC JAC C12 C6 [CS12 CS6])
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        if parts.len() >= 4 {
            let iac: usize = parts[0]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid IAC: {}", parts[0])))?;
            let jac: usize = parts[1]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid JAC: {}", parts[1])))?;
            let c12: f64 = parts[2]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid C12: {}", parts[2])))?;
            let c6: f64 = parts[3]
                .parse()
                .map_err(|_| IoError::ParseError(format!("Invalid C6: {}", parts[3])))?;

            lj_params.insert((iac, jac), LJParameters::new(c6, c12));

            // LJ matrix is symmetric
            if iac != jac {
                lj_params.insert((jac, iac), LJParameters::new(c6, c12));
            }
        }
    }

    Ok(())
}

fn parse_temperature_groups<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    temp_groups: &mut Vec<usize>,
) -> Result<(), IoError> {
    let mut n_groups = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }

        if trimmed == "END" {
            break;
        }

        // First line is group count
        if n_groups == 0 {
            n_groups = trimmed.parse().map_err(|_| {
                IoError::ParseError(format!("Invalid temp group count: {}", trimmed))
            })?;
            continue;
        }

        // Next line(s) contain last atom indices for each group
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        for part in parts {
            if let Ok(last_atom) = part.parse::<usize>() {
                temp_groups.push(last_atom - 1); // Convert to 0-indexed
            }
        }
    }

    Ok(())
}

fn parse_solventatom<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    solvent_atoms: &mut Vec<ParsedSolventAtom>,
) -> Result<(), IoError> {
    let mut nram: usize = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        if trimmed == "END" { break; }

        if nram == 0 {
            nram = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid NRAM: {}", trimmed)))?;
            if nram == 0 { break; }
            continue;
        }

        // Format: I ANMS IACS MASS CGS
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 5 {
            let atom_iac: usize = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid solvent IAC: {}", parts[2])))?;
            let mass: f64 = parts[3].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid solvent mass: {}", parts[3])))?;
            let charge: f64 = parts[4].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid solvent charge: {}", parts[4])))?;
            solvent_atoms.push(ParsedSolventAtom {
                iac: atom_iac,
                name: parts[1].to_string(),
                mass,
                charge,
            });
        }
    }
    Ok(())
}

fn parse_solventconstr<I: Iterator<Item = Result<String, std::io::Error>>>(
    lines: &mut I,
    constraints: &mut Vec<ParsedSolventConstraint>,
) -> Result<(), IoError> {
    let mut ncons: usize = 0;

    while let Some(Ok(line)) = lines.next() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') { continue; }
        if trimmed == "END" { break; }

        if ncons == 0 {
            ncons = trimmed.parse()
                .map_err(|_| IoError::ParseError(format!("Invalid NCONS: {}", trimmed)))?;
            if ncons == 0 { break; }
            continue;
        }

        // Format: ICONS JCONS CONS
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() >= 3 {
            let i: usize = parts[0].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid ICONS: {}", parts[0])))?;
            let j: usize = parts[1].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid JCONS: {}", parts[1])))?;
            let length: f64 = parts[2].parse()
                .map_err(|_| IoError::ParseError(format!("Invalid CONS: {}", parts[2])))?;
            constraints.push(ParsedSolventConstraint {
                i: i - 1, // 0-indexed within molecule
                j: j - 1,
                length,
            });
        }
    }
    Ok(())
}

/// Convert ParsedTopology to Topology.
///
/// Builds solute topology and stores the solvent template.
/// To expand solvent molecules, call `topo.solvate(nsm)` afterwards
/// (gromosXX convention: NSM comes from the IMD SYSTEM block).
pub fn build_topology(parsed: ParsedTopology) -> Topology {
    use gromos_core::topology::{Atom, ChargeGroup, SolventAtomTemplate, SolventConstraintTemplate};

    let mut topo = Topology::new();
    let n_solute = parsed.n_atoms;

    // --- Populate solute atoms ---
    topo.mass = parsed.masses.clone();
    topo.charge = parsed.charges.clone();
    topo.iac = parsed.iac.clone();

    for i in 0..n_solute {
        topo.solute.atoms.push(Atom {
            name: format!("ATOM{}", i + 1),
            residue_nr: 1,
            residue_name: "UNK".to_string(),
            iac: parsed.iac[i],
            mass: parsed.masses[i],
            charge: parsed.charges[i],
            is_perturbed: false,
            is_polarisable: false,
            is_coarse_grained: true,
        });
    }

    // --- Store solvent template (expanded later by topo.solvate(nsm)) ---
    for sa in &parsed.solvent_atoms {
        topo.solvent_atom_template.push(SolventAtomTemplate {
            iac: sa.iac,
            name: sa.name.clone(),
            mass: sa.mass,
            charge: sa.charge,
        });
    }
    for sc in &parsed.solvent_constraints {
        topo.solvent_constraint_template.push(SolventConstraintTemplate {
            i: sc.i,
            j: sc.j,
            length: sc.length,
        });
    }

    topo.compute_inverse_masses();

    // --- Store chargegroup codes and build solute chargegroups ---
    topo.chargegroup_codes = parsed.chargegroup_codes.clone();
    {
        let mut current_cg_atoms = Vec::new();
        for i in 0..n_solute {
            current_cg_atoms.push(i);
            let cgc = if i < parsed.chargegroup_codes.len() {
                parsed.chargegroup_codes[i]
            } else {
                1 // default: each atom is its own chargegroup
            };
            if cgc == 1 {
                topo.chargegroups.push(ChargeGroup { atoms: current_cg_atoms.clone() });
                current_cg_atoms.clear();
            }
        }
        if !current_cg_atoms.is_empty() {
            topo.chargegroups.push(ChargeGroup { atoms: current_cg_atoms });
        }

        // Build atom_to_chargegroup mapping (solute only at this point)
        topo.atom_to_chargegroup = vec![0; n_solute];
        for (cg_idx, cg) in topo.chargegroups.iter().enumerate() {
            for &atom in &cg.atoms {
                if atom < n_solute {
                    topo.atom_to_chargegroup[atom] = cg_idx;
                }
            }
        }
    }

    // --- Initialize exclusions for solute atoms ---
    topo.exclusions = vec![std::collections::HashSet::new(); n_solute];

    for (i, excl_list) in parsed.exclusions.iter().enumerate() {
        for &j in excl_list {
            topo.exclusions[i].insert(j);
            if j < n_solute {
                topo.exclusions[j].insert(i);
            }
        }
    }

    // Build bonds
    for (i, j, bond_type) in parsed.bonds {
        topo.solute.bonds.push(Bond { i, j, bond_type });
        topo.exclusions[i].insert(j);
        topo.exclusions[j].insert(i);
    }

    topo.bond_parameters = parsed.bond_parameters;

    // Build angles
    for (i, j, k, angle_type) in parsed.angles {
        topo.solute.angles.push(Angle {
            i, j, k, angle_type,
        });
        topo.exclusions[i].insert(k);
        topo.exclusions[k].insert(i);
    }

    topo.angle_parameters = parsed.angle_parameters;

    // Build proper dihedrals
    for (i, j, k, l, dihedral_type) in parsed.proper_dihedrals {
        topo.solute.proper_dihedrals.push(Dihedral {
            i, j, k, l, dihedral_type,
        });
    }
    topo.dihedral_parameters = parsed.dihedral_parameters;

    // Build improper dihedrals
    for (i, j, k, l, dihedral_type) in parsed.improper_dihedrals {
        topo.solute.improper_dihedrals.push(Dihedral {
            i, j, k, l, dihedral_type,
        });
    }
    topo.improper_dihedral_parameters = parsed.improper_dihedral_parameters;

    // Build LJ parameter matrix
    let max_iac = topo.iac.iter()
        .chain(topo.solvent_atom_template.iter().map(|sa| &sa.iac))
        .max().copied().unwrap_or(0);
    let n_types = max_iac + 1;
    topo.lj_parameters = vec![vec![LJParameters::default(); n_types]; n_types];

    for ((iac, jac), params) in parsed.lj_parameters {
        if iac < n_types && jac < n_types {
            topo.lj_parameters[iac][jac] = params;
            topo.lj_parameters[jac][iac] = params;
            log::debug!("LJ params[{}][{}]: c6={:.6e}, c12={:.6e}", iac, jac, params.c6, params.c12);
        }
    }

    log::debug!("Built topology: {} solute atoms, {} chargegroups, LJ {}x{}, solvent template: {} atoms",
        n_solute, topo.chargegroups.len(), n_types, n_types,
        topo.solvent_atom_template.len());

    topo
}

/// Write GROMOS topology file
pub fn write_topology_file<P: AsRef<Path>>(
    path: P,
    topo: &Topology,
    title: &str,
) -> Result<(), IoError> {
    use std::io::Write;

    let file = File::create(path.as_ref())
        .map_err(|e| IoError::WriteError(format!("Cannot create topology file: {}", e)))?;
    let mut writer = std::io::BufWriter::new(file);

    // TITLE block
    writeln!(writer, "TITLE").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "{}", title).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // PHYSICALCONSTANTS block (optional but good practice)
    writeln!(writer, "PHYSICALCONSTANTS").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(
        writer,
        "# FPEPSI: 1.0/(4.0*PI*EPS0) (EPS0 is the permittivity of vacuum)"
    )
    .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  138.9354").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# HBAR: Planck's constant HBAR = H/(2* PI)")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  0.0635078").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# SPDL: Speed of light (nm/ps)")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  299792.458").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# BOLTZ: Boltzmann's constant kB")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "  0.00831441").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // TOPVERSION block
    writeln!(writer, "TOPVERSION").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "2.0").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // ATOMTYPENAME block
    writeln!(writer, "ATOMTYPENAME").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# NTYP: number of atom types")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    let n_types = topo.num_atom_types();
    writeln!(writer, "{}", n_types).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# TYPE: atom type names").map_err(|e| IoError::WriteError(e.to_string()))?;
    for i in 1..=n_types {
        writeln!(writer, "TYPE{}", i).map_err(|e| IoError::WriteError(e.to_string()))?;
    }
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // RESNAME block
    writeln!(writer, "RESNAME").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# NRES: number of residues")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "1").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# residue name").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "MOL").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // SOLUTEATOM block
    let n_atoms = topo.solute.num_atoms();
    writeln!(writer, "SOLUTEATOM").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "#   NRP: number of solute atoms")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "{}", n_atoms).map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# ATNM: atom number").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# MRES: residue number").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# PANM: atom name").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# IAC: integer atom code").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# MASS: mass").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# CG: charge").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# CGC: charge group code").map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "# INE: number of exclusions")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    writeln!(writer, "#ATNM MRES PANM IAC  MASS      CG     CGC INE")
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    for (i, atom) in topo.solute.atoms.iter().enumerate() {
        let n_exclusions = topo.exclusions.get(i).map_or(0, |e| e.len());
        let cg_code = topo.atom_to_chargegroup.get(i).map_or(1, |&c| c + 1);

        writeln!(
            writer,
            "{:5} {:4} {:>4} {:3} {:8.4} {:8.5} {:4} {:3}",
            i + 1,
            atom.residue_nr,
            atom.name,
            atom.iac,
            atom.mass,
            atom.charge,
            cg_code,
            n_exclusions
        )
        .map_err(|e| IoError::WriteError(e.to_string()))?;

        // Write exclusions if any
        if let Some(exclusions) = topo.exclusions.get(i) {
            if !exclusions.is_empty() {
                let mut excl_vec: Vec<_> = exclusions.iter().collect();
                excl_vec.sort();
                for (j, excl) in excl_vec.iter().enumerate() {
                    write!(writer, "{:5}", *excl + 1)
                        .map_err(|e| IoError::WriteError(e.to_string()))?;
                    if (j + 1) % 10 == 0 {
                        writeln!(writer).map_err(|e| IoError::WriteError(e.to_string()))?;
                    }
                }
                writeln!(writer).map_err(|e| IoError::WriteError(e.to_string()))?;
            }
        }
    }
    writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;

    // BONDSTRETCHTYPE block
    if !topo.bond_parameters.is_empty() {
        writeln!(writer, "BONDSTRETCHTYPE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NBTY: number of bond types")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.bond_parameters.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  CB     CHB       B0")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        for params in &topo.bond_parameters {
            writeln!(
                writer,
                "{:10.5e} {:10.5e} {:10.7}",
                params.k_quartic, params.k_harmonic, params.r0
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BOND block
    if !topo.solute.bonds.is_empty() {
        writeln!(writer, "BOND").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NBH: number of bonds")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.solute.bonds.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  IB   JB  ICB").map_err(|e| IoError::WriteError(e.to_string()))?;
        for bond in &topo.solute.bonds {
            writeln!(
                writer,
                "{:5}{:5}{:5}",
                bond.i + 1,
                bond.j + 1,
                bond.bond_type + 1
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BONDANGLEBENDTYPE block
    if !topo.angle_parameters.is_empty() {
        writeln!(writer, "BONDANGLEBENDTYPE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NTTY: number of angle types")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.angle_parameters.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#   CT      CHT     T0[deg]")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        for params in &topo.angle_parameters {
            writeln!(
                writer,
                "{:10.5e} {:10.5e} {:10.4}",
                params.k_cosine,
                params.k_harmonic,
                params.theta0.to_degrees()
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    // BONDANGLE block
    if !topo.solute.angles.is_empty() {
        writeln!(writer, "BONDANGLE").map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "# NTHEH: number of angles")
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "{}", topo.solute.angles.len())
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        writeln!(writer, "#  IT   JT   KT  ICT").map_err(|e| IoError::WriteError(e.to_string()))?;
        for angle in &topo.solute.angles {
            writeln!(
                writer,
                "{:5}{:5}{:5}{:5}",
                angle.i + 1,
                angle.j + 1,
                angle.k + 1,
                angle.angle_type + 1
            )
            .map_err(|e| IoError::WriteError(e.to_string()))?;
        }
        writeln!(writer, "END").map_err(|e| IoError::WriteError(e.to_string()))?;
    }

    writer
        .flush()
        .map_err(|e| IoError::WriteError(e.to_string()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Minimal 4-atom GROMOS topology (cg16-style).
    /// SOLUTEATOM format: ATNM MRES PANM IAC MASS CG CGC INE IEXCL...
    const CG16_TOPO: &str = "\
TITLE
  4-atom test topology (cg16)
END
SOLUTEATOM
# ATNM MRES PANM IAC   MASS      CG   CGC INE IEXCL
4
    1    1   C     6  62.000000  0.000000  0  2  2  3
    2    1   H     1  12.000000  0.000000  0  1  1
    3    1   H     1  12.000000  0.000000  0  1  1
    4    1   H     1  12.000000  0.000000  0  0
END
BONDSTRETCHTYPE
# CB        CHB       B0
1
  1000.000  500.000   0.150
END
BOND
# IB  JB  ICB
3
  1   2   1
  2   3   1
  3   4   1
END
BONDANGLEBENDTYPE
# CT        CHT       T0
1
   50.000  100.000  109.500
END
BONDANGLE
# IT  JT  KT  ICT
2
  1   2   3   1
  2   3   4   1
END
";

    fn write_tmp(content: &str, suffix: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!("gromos_test_{suffix}.tmp"));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn test_parse_cg16_topology() {
        let path = write_tmp(CG16_TOPO, "cg16_topo");
        let parsed = read_topology_file(&path).expect("Failed to parse inline topology");

        assert_eq!(parsed.n_atoms, 4);
        assert_eq!(parsed.bonds.len(), 3);
        assert_eq!(parsed.angles.len(), 2);

        // First atom: IAC=6, mass=62.0, charge=0.0
        assert_eq!(parsed.iac[0], 6);
        assert!((parsed.masses[0] - 62.0).abs() < 1e-6);
        assert!((parsed.charges[0] - 0.0).abs() < 1e-6);

        // First bond: atoms 1-2 in file → 0-1 (0-indexed), type 1 → 0
        assert_eq!(parsed.bonds[0], (0, 1, 0));

        // Bond parameters
        assert!((parsed.bond_parameters[0].k_harmonic - 500.0).abs() < 1e-6);
        assert!((parsed.bond_parameters[0].r0 - 0.15).abs() < 1e-6);

        std::fs::remove_file(path).ok();
    }

    #[test]
    fn test_parse_missing_topology_returns_error() {
        let result = read_topology_file("/nonexistent/path/file.topo");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_topology_bonds_are_zero_indexed() {
        let path = write_tmp(CG16_TOPO, "cg16_topo_idx");
        let parsed = read_topology_file(&path).unwrap();

        // All bond indices must be 0-indexed (< n_atoms)
        for (i, j, _t) in &parsed.bonds {
            assert!(*i < parsed.n_atoms);
            assert!(*j < parsed.n_atoms);
        }

        std::fs::remove_file(path).ok();
    }
}
