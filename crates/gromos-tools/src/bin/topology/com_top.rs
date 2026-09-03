//! com_top - Combine multiple GROMOS topology files
//!
//! Usage: com_top @topo <file1> [N:file2] ... @param <n> @solv <n>
//!    or: com_top @f <argfile>
//!
//! Combines multiple topology files into a single topology.
//! Use prefix 'N:' before a file to repeat that topology N times.
//! Parameters and solvent are taken from the topology indicated by @param/@solv.

use clap::Parser;
use gromos_io::gromos_args;
use gromos_io::topology::{read_topology_file, write_parsed_topology, ParsedTopology};
use std::collections::HashMap;
use std::process;
/// Combine multiple GROMOS topology files.
///
/// Supports GROMOS @-prefix args: com_top @topo GB3.top 2:Na.top @param 1 @solv 1
#[derive(Parser)]
#[command(name = "com_top", version, about)]
struct Args {
    /// Topology files (prefix N: to repeat, e.g. 2:Na.top)
    #[arg(long, num_args = 1..)]
    topo: Vec<String>,

    /// Index (1-based) of topology to take force field parameters from
    #[arg(long, default_value = "1")]
    param: usize,

    /// Index (1-based) of topology to take solvent from
    #[arg(long, default_value = "1")]
    solv: usize,
}

/// Parse a topology spec: "file.top" or "N:file.top"
fn parse_topo_spec(spec: &str) -> (usize, String) {
    if let Some(pos) = spec.find(':') {
        let count: usize = spec[..pos].parse().unwrap_or_else(|_| {
            eprintln!("Error: invalid repeat count in '{}'", spec);
            process::exit(1);
        });
        if count == 0 {
            eprintln!("Error: repeat count 0 not allowed in '{}'", spec);
            process::exit(1);
        }
        (count, spec[pos + 1..].to_string())
    } else {
        (1, spec.to_string())
    }
}

fn main() {
    let args = Args::parse_from(gromos_args());

    if args.topo.is_empty() {
        eprintln!("Error: @topo requires at least one topology file");
        process::exit(1);
    }

    // Parse topology specs and count total
    let specs: Vec<(usize, String)> = args.topo.iter().map(|s| parse_topo_spec(s)).collect();
    let total: usize = specs.iter().map(|(n, _)| n).sum();

    if args.param < 1 || args.param > total {
        eprintln!("Error: @param {} out of range 1..{}", args.param, total);
        process::exit(1);
    }
    if args.solv < 1 || args.solv > total {
        eprintln!("Error: @solv {} out of range 1..{}", args.solv, total);
        process::exit(1);
    }

    // Read topology files (read each unique file once)
    let mut file_cache: HashMap<String, ParsedTopology> = HashMap::new();
    for (_, filename) in &specs {
        if !file_cache.contains_key(filename) {
            let parsed = match read_topology_file(filename) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!("Error reading '{}': {:?}", filename, e);
                    process::exit(1);
                },
            };
            eprintln!(
                "  Read '{}': {} atoms, {} residues",
                filename,
                parsed.n_atoms,
                parsed.residue_names.len()
            );
            file_cache.insert(filename.clone(), parsed);
        }
    }

    // Determine which file provides params and solvent
    let mut count = 0;
    let mut param_file = String::new();
    let mut solv_file = String::new();
    for (repeat, filename) in &specs {
        for _ in 0..*repeat {
            count += 1;
            if count == args.param {
                param_file = filename.clone();
            }
            if count == args.solv {
                solv_file = filename.clone();
            }
        }
    }

    // Build title
    let mut title = String::from("COM_TOP: Combined topology using:\n");
    let mut c = 0;
    for (repeat, filename) in &specs {
        let start = c + 1;
        c += repeat;
        if *repeat > 1 {
            title.push_str(&format!("  {} .. {} : {}\n", start, c, filename));
        } else {
            title.push_str(&format!("  {} : {}\n", start, filename));
        }
    }

    // Combine topologies
    let combined = combine_topologies(&specs, &file_cache, &param_file, &solv_file);

    // Write to stdout
    print!("{}", write_parsed_topology(&combined, &title));
}

fn combine_topologies(
    specs: &[(usize, String)],
    file_cache: &HashMap<String, ParsedTopology>,
    param_file: &str,
    solv_file: &str,
) -> ParsedTopology {
    let param_topo = &file_cache[param_file];
    let solv_topo = &file_cache[solv_file];

    let mut combined = ParsedTopology {
        physical_constants: param_topo.physical_constants,
        n_atoms: 0,
        atom_names: Vec::new(),
        residue_numbers: Vec::new(),
        residue_names: Vec::new(),
        atom_type_names: param_topo.atom_type_names.clone(),
        masses: Vec::new(),
        charges: Vec::new(),
        iac: Vec::new(),
        chargegroup_codes: Vec::new(),
        exclusions: Vec::new(),
        one_four_pairs: Vec::new(),
        bonds: Vec::new(),
        bond_parameters: param_topo.bond_parameters.clone(),
        angles: Vec::new(),
        angle_parameters: param_topo.angle_parameters.clone(),
        proper_dihedrals: Vec::new(),
        dihedral_parameters: param_topo.dihedral_parameters.clone(),
        improper_dihedrals: Vec::new(),
        improper_dihedral_parameters: param_topo.improper_dihedral_parameters.clone(),
        lj_parameters: param_topo.lj_parameters.clone(),
        temperature_groups: Vec::new(),
        pressure_groups: Vec::new(),
        solvent_atoms: solv_topo.solvent_atoms.clone(),
        solvent_constraints: solv_topo.solvent_constraints.clone(),
        solute_molecules: Vec::new(),
    };

    let mut atom_offset: usize = 0;
    let mut residue_offset: usize = 0;

    for (repeat, filename) in specs {
        let topo = &file_cache[filename];

        for _ in 0..*repeat {
            // Merge residue names
            combined
                .residue_names
                .extend(topo.residue_names.iter().cloned());

            // Merge atoms
            for i in 0..topo.n_atoms {
                combined.atom_names.push(topo.atom_names[i].clone());
                combined
                    .residue_numbers
                    .push(topo.residue_numbers[i] + residue_offset);
                combined.masses.push(topo.masses[i]);
                combined.charges.push(topo.charges[i]);
                combined.iac.push(topo.iac[i]);
                combined.chargegroup_codes.push(topo.chargegroup_codes[i]);

                // Renumber exclusions
                let excl = topo.exclusions.get(i).cloned().unwrap_or_default();
                combined
                    .exclusions
                    .push(excl.iter().map(|&e| e + atom_offset).collect());

                // Renumber 1-4 pairs
                let pairs = topo.one_four_pairs.get(i).cloned().unwrap_or_default();
                combined
                    .one_four_pairs
                    .push(pairs.iter().map(|&p| p + atom_offset).collect());
            }

            // Merge bonds with renumbered atoms
            for &(i, j, t) in &topo.bonds {
                combined.bonds.push((i + atom_offset, j + atom_offset, t));
            }

            // Merge angles
            for &(i, j, k, t) in &topo.angles {
                combined
                    .angles
                    .push((i + atom_offset, j + atom_offset, k + atom_offset, t));
            }

            // Merge proper dihedrals
            for &(i, j, k, l, t) in &topo.proper_dihedrals {
                combined.proper_dihedrals.push((
                    i + atom_offset,
                    j + atom_offset,
                    k + atom_offset,
                    l + atom_offset,
                    t,
                ));
            }

            // Merge improper dihedrals
            for &(i, j, k, l, t) in &topo.improper_dihedrals {
                combined.improper_dihedrals.push((
                    i + atom_offset,
                    j + atom_offset,
                    k + atom_offset,
                    l + atom_offset,
                    t,
                ));
            }

            // Temperature / pressure groups: renumber
            for &tg in &topo.temperature_groups {
                combined.temperature_groups.push(tg + atom_offset);
            }
            for &pg in &topo.pressure_groups {
                combined.pressure_groups.push(pg + atom_offset);
            }

            // Solute molecules: renumber
            for &mol_last in &topo.solute_molecules {
                combined.solute_molecules.push(mol_last + atom_offset);
            }

            residue_offset += topo.residue_names.len();
            atom_offset += topo.n_atoms;
        }
    }

    combined.n_atoms = combined.masses.len();

    eprintln!(
        "Combined topology: {} atoms, {} residues, {} molecules",
        combined.n_atoms,
        combined.residue_names.len(),
        combined.solute_molecules.len()
    );

    combined
}
