//! pdb2g96 - Convert PDB files to GROMOS96 format (topology-guided)
//!
//! Usage: pdb2g96 @topo <top> @pdb <pdb> @lib <lib> [@gch]
//!    or: pdb2g96 @f <argfile>
//!
//! Reads a topology and a PDB file, then writes GROMOS96 coordinates
//! with atoms in topology order. Uses a library file for name mapping.
//! The @gch flag requests hydrogen coordinate generation (missing H atoms
//! get (0,0,0) and a warning; future: geometry-based generation).

use clap::Parser;
use gromos_io::gromos_args;
use gromos_io::pdb::PDBStructure;
use gromos_io::topology::read_topology_file;
use std::collections::HashMap;
use std::process;

#[derive(Parser)]
#[command(name = "pdb2g96", version, about)]
struct Args {
    /// Molecular topology file
    #[arg(long)]
    topo: String,

    /// Input PDB file
    #[arg(long)]
    pdb: String,

    /// Library for atom and residue name mapping
    #[arg(long)]
    lib: Option<String>,

    /// (Re)generate hydrogen coordinates
    #[arg(long, action = clap::ArgAction::SetTrue)]
    gch: bool,
}

/// Library mapping: PDB residue/atom names → topology names.
struct Library {
    residue_map: HashMap<String, String>, // pdb_res → topo_res
    atom_map: HashMap<String, HashMap<String, String>>, // topo_res → {pdb_atom → topo_atom}
}

impl Library {
    fn new() -> Self {
        Self {
            residue_map: HashMap::new(),
            atom_map: HashMap::new(),
        }
    }

    fn read(path: &str) -> Result<Self, String> {
        let content = std::fs::read_to_string(path)
            .map_err(|e| format!("Cannot read library file '{}': {}", path, e))?;
        let mut lib = Self::new();
        let mut in_block = "";

        for line in content.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if trimmed == "END" {
                in_block = "";
                continue;
            }
            if trimmed == "TITLE" {
                in_block = "TITLE";
                continue;
            }
            if trimmed == "RESIDUENAMELIB" || trimmed == "RESIDUES" {
                in_block = "RESIDUES";
                continue;
            }
            if trimmed == "ATOMNAMELIB" || trimmed == "ATOMS" {
                in_block = "ATOMS";
                continue;
            }
            match in_block {
                "RESIDUES" => {
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 2 {
                        lib.residue_map
                            .insert(parts[0].to_string(), parts[1].to_string());
                    }
                },
                "ATOMS" => {
                    // Format: topo_residue pdb_atom topo_atom
                    let parts: Vec<&str> = trimmed.split_whitespace().collect();
                    if parts.len() >= 3 {
                        lib.atom_map
                            .entry(parts[0].to_string())
                            .or_default()
                            .insert(parts[1].to_string(), parts[2].to_string());
                    }
                },
                _ => {},
            }
        }
        Ok(lib)
    }

    /// Map a PDB residue name to topology residue name.
    fn map_residue<'a>(&'a self, pdb_name: &'a str) -> &'a str {
        self.residue_map
            .get(pdb_name)
            .map(|s| s.as_str())
            .unwrap_or(pdb_name)
    }

    /// Map a PDB atom name to topology atom name for a given topology residue.
    fn map_atom<'a>(&'a self, topo_res: &str, pdb_atom: &'a str) -> &'a str {
        self.atom_map
            .get(topo_res)
            .and_then(|m| m.get(pdb_atom))
            .map(|s| s.as_str())
            .unwrap_or(pdb_atom)
    }
}

/// A PDB atom with its coordinates (already in nm from the PDB parser).
struct PdbAtomEntry {
    name: String,
    x: f64,
    y: f64,
    z: f64,
}

/// A PDB residue: name + list of atoms.
struct PdbResidue {
    name: String,
    atoms: Vec<PdbAtomEntry>,
}

/// Match topology atoms against a PDB residue's atoms, assigning coordinates.
/// Returns (matched, missing) counts.
fn match_atoms(
    topo_atoms: &[(usize, String)],
    topo_res_name: &str,
    pdb_res: &PdbResidue,
    used: &mut Vec<bool>,
    lib: &Library,
    coords: &mut [(f64, f64, f64)],
    masses: &[f64],
    gch: bool,
) -> (usize, usize) {
    let mut matched = 0;
    let mut missing = 0;
    for (global_idx, topo_atom_name) in topo_atoms {
        let mut found = false;
        for (j, pdb_atom) in pdb_res.atoms.iter().enumerate() {
            if used[j] {
                continue;
            }
            let mapped_name = lib.map_atom(topo_res_name, &pdb_atom.name);
            if mapped_name == topo_atom_name || pdb_atom.name == *topo_atom_name {
                coords[*global_idx] = (pdb_atom.x, pdb_atom.y, pdb_atom.z);
                used[j] = true;
                found = true;
                matched += 1;
                break;
            }
        }
        if !found {
            let is_hydrogen = masses[*global_idx] < 2.0;
            if !gch || !is_hydrogen {
                eprintln!(
                    "Warning: Could not find atom '{}' in residue '{}' — set to (0,0,0)",
                    topo_atom_name, topo_res_name
                );
            }
            missing += 1;
        }
    }
    (matched, missing)
}

fn main() {
    let args = Args::parse_from(gromos_args());

    // Read topology
    let parsed_topo = match read_topology_file(&args.topo) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error reading topology '{}': {:?}", args.topo, e);
            process::exit(1);
        },
    };

    eprintln!(
        "Topology: {} atoms, {} residues",
        parsed_topo.n_atoms,
        parsed_topo.residue_names.len()
    );

    // Read PDB
    let pdb = match PDBStructure::read_pdb(&args.pdb) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error reading PDB '{}': {}", args.pdb, e);
            process::exit(1);
        },
    };

    eprintln!("PDB: {} residues", pdb.residues.len());

    // Read library
    let lib = if let Some(ref lib_path) = args.lib {
        match Library::read(lib_path) {
            Ok(l) => l,
            Err(e) => {
                eprintln!("Error: {}", e);
                process::exit(1);
            },
        }
    } else {
        Library::new()
    };

    // Group PDB atoms by residue
    let mut pdb_residues: Vec<PdbResidue> = Vec::new();
    for res in &pdb.residues {
        let mapped_name = lib.map_residue(res.name.trim());
        pdb_residues.push(PdbResidue {
            name: mapped_name.to_string(),
            atoms: res
                .atoms
                .iter()
                .map(|a| PdbAtomEntry {
                    name: a.name.trim().to_string(),
                    x: a.coord.x,
                    y: a.coord.y,
                    z: a.coord.z,
                })
                .collect(),
        });
    }

    // Build topology atom list: for each atom, (residue_number, residue_name, atom_name)
    // Group topology atoms by residue for matching
    let n_atoms = parsed_topo.n_atoms;
    let mut topo_residues: Vec<(usize, String, Vec<(usize, String)>)> = Vec::new();
    // (res_number, res_name, [(global_atom_idx, atom_name)])
    {
        let mut current_res_nr = 0;
        for i in 0..n_atoms {
            let res_nr = parsed_topo.residue_numbers[i];
            let atom_name = &parsed_topo.atom_names[i];
            if res_nr != current_res_nr {
                let res_name = if res_nr > 0 && res_nr <= parsed_topo.residue_names.len() {
                    parsed_topo.residue_names[res_nr - 1].clone()
                } else {
                    "UNK".to_string()
                };
                topo_residues.push((res_nr, res_name, Vec::new()));
                current_res_nr = res_nr;
            }
            if let Some(last) = topo_residues.last_mut() {
                last.2.push((i, atom_name.clone()));
            }
        }
    }

    // Match: for each topology residue, find matching PDB residue and assign coordinates
    // GROMOS topologies can have end-group residues (NH3+, COO-) that are not
    // separate residues in the PDB. These atoms belong to the adjacent PDB residue.
    let mut coords: Vec<(f64, f64, f64)> = vec![(0.0, 0.0, 0.0); n_atoms];
    let mut pdb_idx = 0; // current PDB residue index
    let mut missing_count = 0;
    let mut matched_count = 0;
    // Track which PDB atoms have been used (by pdb_idx)
    let mut used_map: HashMap<usize, Vec<bool>> = HashMap::new();

    for (_, topo_res_name, topo_atoms) in &topo_residues {
        if pdb_idx >= pdb_residues.len() {
            // No more PDB residues — try the last one for C-terminal end-groups
            if pdb_idx > 0 {
                let prev_idx = pdb_idx - 1;
                let pdb_res = &pdb_residues[prev_idx];
                let used = used_map
                    .entry(prev_idx)
                    .or_insert_with(|| vec![false; pdb_res.atoms.len()]);
                eprintln!("Note: Topology residue '{}' — looking for atoms in previous PDB residue '{}' (C-terminal end-group)",
                    topo_res_name, pdb_res.name);
                let (m, mi) = match_atoms(
                    topo_atoms,
                    topo_res_name,
                    pdb_res,
                    used,
                    &lib,
                    &mut coords,
                    &parsed_topo.masses,
                    args.gch,
                );
                matched_count += m;
                missing_count += mi;
            } else {
                eprintln!(
                    "Warning: No more PDB residues for topology residue '{}'",
                    topo_res_name
                );
                missing_count += topo_atoms.len();
            }
            continue;
        }

        let pdb_res = &pdb_residues[pdb_idx];

        // Check if this topology residue matches the current PDB residue
        let names_match = pdb_res.name == *topo_res_name;

        // If names don't match, this might be an end-group residue (NH3+, COO-)
        // whose atoms are part of the current PDB residue. Don't advance pdb_idx.
        // If names match, we'll advance pdb_idx after matching.
        let advance_pdb = names_match;

        if !names_match {
            eprintln!(
                "Note: Topology residue '{}' — looking for atoms in PDB residue '{}' (end-group)",
                topo_res_name, pdb_res.name
            );
        }

        // Get or create the used-atoms tracker for this PDB residue
        let used = used_map
            .entry(pdb_idx)
            .or_insert_with(|| vec![false; pdb_res.atoms.len()]);

        for (global_idx, topo_atom_name) in topo_atoms {
            let mut found = false;

            for (j, pdb_atom) in pdb_res.atoms.iter().enumerate() {
                if used[j] {
                    continue;
                }

                let mapped_name = lib.map_atom(topo_res_name, &pdb_atom.name);
                if mapped_name == topo_atom_name || pdb_atom.name == *topo_atom_name {
                    coords[*global_idx] = (pdb_atom.x, pdb_atom.y, pdb_atom.z);
                    used[j] = true;
                    found = true;
                    matched_count += 1;
                    break;
                }
            }

            if !found {
                let is_hydrogen = parsed_topo.masses[*global_idx] < 2.0;
                if !args.gch || !is_hydrogen {
                    eprintln!(
                        "Warning: Could not find atom '{}' in residue '{}' — set to (0,0,0)",
                        topo_atom_name, topo_res_name
                    );
                }
                missing_count += 1;
            }
        }

        if advance_pdb {
            // Warn about unused PDB atoms in this residue
            let used = used_map.get(&pdb_idx).unwrap();
            for (j, pdb_atom) in pdb_res.atoms.iter().enumerate() {
                if !used[j] {
                    eprintln!(
                        "Warning: Ignored PDB atom '{}' in residue '{}'",
                        pdb_atom.name, pdb_res.name
                    );
                }
            }
            pdb_idx += 1;
        }
    }

    eprintln!(
        "Matched {} atoms, {} missing (set to 0,0,0)",
        matched_count, missing_count
    );

    // Write GROMOS96 POSITION block to stdout
    println!("TITLE");
    println!("    pdb2g96: Reordered atoms from {}", args.pdb);
    if missing_count > 0 {
        println!("    {} atoms could not be matched", missing_count);
    }
    println!("END");

    println!("POSITION");
    for i in 0..n_atoms {
        let res_nr = parsed_topo.residue_numbers[i];
        let res_name = if res_nr > 0 && res_nr <= parsed_topo.residue_names.len() {
            &parsed_topo.residue_names[res_nr - 1]
        } else {
            "UNK"
        };
        let atom_name = &parsed_topo.atom_names[i];
        let (x, y, z) = coords[i];
        println!(
            "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
            res_nr,
            res_name,
            atom_name,
            i + 1,
            x,
            y,
            z
        );
    }
    println!("END");
}
