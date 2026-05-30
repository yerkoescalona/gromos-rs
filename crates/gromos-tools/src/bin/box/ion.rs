//! ion - Replace solvent molecules by ions
//!
//! Usage: ion @topo <top> @pbc <type> @pos <cnf> @positive <N> <name> [@negative <N> <name>]
//!            @potential <cutoff> [@mindist <d>]
//!    or: ion @f <argfile>
//!
//! Places counter-ions by replacing solvent molecules at positions with the
//! lowest (for positive ions) or highest (for negative ions) electrostatic
//! potential from solute charges.

use clap::Parser;
use gromos_io::{gromos_args, read_g96_labeled, G96Atom};
use gromos_io::topology::read_topology_file;
use std::collections::HashSet;
use std::process;

#[derive(Parser)]
#[command(name = "ion", version, about = "Replace solvent molecules by ions")]
struct Args {
    /// Molecular topology file
    #[arg(long)]
    topo: String,

    /// Periodic boundary conditions (r = rectangular)
    #[arg(long)]
    pbc: String,

    /// Input coordinate file (solvated system)
    #[arg(long)]
    pos: String,

    /// Place positive ions: <number> <ion_name>
    #[arg(long, num_args = 2)]
    positive: Option<Vec<String>>,

    /// Place negative ions: <number> <ion_name>
    #[arg(long, num_args = 2)]
    negative: Option<Vec<String>>,

    /// Cutoff for electrostatic potential calculation (nm)
    #[arg(long)]
    potential: Option<f64>,

    /// Minimum distance between ions (nm)
    #[arg(long, default_value = "0.0")]
    mindist: f64,
}

/// Minimum image distance squared under rectangular PBC.
fn dist_sq_pbc(a: (f64, f64, f64), b: (f64, f64, f64), box_dim: &[f64; 3]) -> f64 {
    let mut dx = a.0 - b.0;
    let mut dy = a.1 - b.1;
    let mut dz = a.2 - b.2;
    dx -= (dx / box_dim[0]).round() * box_dim[0];
    dy -= (dy / box_dim[1]).round() * box_dim[1];
    dz -= (dz / box_dim[2]).round() * box_dim[2];
    dx * dx + dy * dy + dz * dz
}

fn main() {
    let args = Args::parse_from(gromos_args());

    // Parse ion counts
    let (n_positive, pos_name) = if let Some(ref v) = args.positive {
        let n: usize = v[0].parse().unwrap_or_else(|_| {
            eprintln!("Error: invalid positive ion count '{}'", v[0]);
            process::exit(1);
        });
        (n, v[1].clone())
    } else {
        (0, String::new())
    };

    let (n_negative, neg_name) = if let Some(ref v) = args.negative {
        let n: usize = v[0].parse().unwrap_or_else(|_| {
            eprintln!("Error: invalid negative ion count '{}'", v[0]);
            process::exit(1);
        });
        (n, v[1].clone())
    } else {
        (0, String::new())
    };

    if n_positive == 0 && n_negative == 0 {
        eprintln!("Error: specify @positive and/or @negative");
        process::exit(1);
    }

    let cutoff = match args.potential {
        Some(c) => c,
        None => {
            eprintln!("Error: @potential cutoff required");
            process::exit(1);
        }
    };

    // Read topology for charges
    let parsed_topo = match read_topology_file(&args.topo) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error reading topology '{}': {:?}", args.topo, e);
            process::exit(1);
        }
    };

    let n_solute = parsed_topo.n_atoms;
    eprintln!("Topology: {} solute atoms", n_solute);

    // Read solvated coordinate file
    let data = match read_g96_labeled(&args.pos) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error reading '{}': {}", args.pos, e);
            process::exit(1);
        }
    };
    let atoms = data.atoms;
    let box_dim = [data.box_dims.x, data.box_dims.y, data.box_dims.z];

    let n_total = atoms.len();
    let n_solvent_atoms = n_total - n_solute;
    eprintln!("Coordinates: {} total atoms ({} solute, {} solvent)",
        n_total, n_solute, n_solvent_atoms);
    eprintln!("Box: ({:.3}, {:.3}, {:.3}) nm", box_dim[0], box_dim[1], box_dim[2]);

    // Identify solvent molecules (groups of 3 atoms for water)
    let atoms_per_solvent = 3; // SPC water
    if n_solvent_atoms % atoms_per_solvent != 0 {
        eprintln!("Warning: solvent atoms ({}) not divisible by {}", n_solvent_atoms, atoms_per_solvent);
    }
    let n_solvent_mols = n_solvent_atoms / atoms_per_solvent;
    eprintln!("Solvent molecules: {}", n_solvent_mols);

    if n_solvent_mols < n_positive + n_negative {
        eprintln!("Error: not enough solvent molecules ({}) for {} ions",
            n_solvent_mols, n_positive + n_negative);
        process::exit(1);
    }

    // Calculate electrostatic potential at each solvent molecule's first atom
    // from all solute charges within the cutoff
    let cutoff_sq = cutoff * cutoff;
    let mut potentials: Vec<f64> = Vec::with_capacity(n_solvent_mols);

    // Coulomb constant: 1/(4*pi*eps0) in GROMOS units (kJ/mol * nm / e^2) = 138.9354
    let coulomb_k = 138.9354;

    for mol_i in 0..n_solvent_mols {
        let ow_idx = n_solute + mol_i * atoms_per_solvent; // first atom (oxygen)
        let ow_pos = (atoms[ow_idx].pos.x, atoms[ow_idx].pos.y, atoms[ow_idx].pos.z);

        let mut phi = 0.0;
        for j in 0..n_solute {
            let q_j = parsed_topo.charges[j];
            if q_j.abs() < 1e-10 {
                continue;
            }
            let atom_pos = (atoms[j].pos.x, atoms[j].pos.y, atoms[j].pos.z);
            let d2 = dist_sq_pbc(ow_pos, atom_pos, &box_dim);
            if d2 < cutoff_sq && d2 > 1e-10 {
                phi += coulomb_k * q_j / d2.sqrt();
            }
        }
        potentials.push(phi);
    }

    // Select solvent molecules to replace
    let mut selected: Vec<usize> = Vec::new(); // indices into solvent molecule list
    let mut excluded: HashSet<usize> = HashSet::new();
    let mindist_sq = args.mindist * args.mindist;

    // For positive ions: pick molecules with LOWEST potential (most negative)
    for _ in 0..n_positive {
        let mut best_idx = None;
        let mut best_pot = f64::MAX;

        for j in 0..n_solvent_mols {
            if excluded.contains(&j) {
                continue;
            }
            if potentials[j] < best_pot {
                best_pot = potentials[j];
                best_idx = Some(j);
            }
        }

        let idx = best_idx.unwrap_or_else(|| {
            eprintln!("Error: not enough suitable solvent positions for positive ions");
            process::exit(1);
        });

        selected.push(idx);
        let ow_i = n_solute + idx * atoms_per_solvent;
        let pos_i = (atoms[ow_i].pos.x, atoms[ow_i].pos.y, atoms[ow_i].pos.z);

        // Exclude nearby solvent molecules
        for j in 0..n_solvent_mols {
            let ow_j = n_solute + j * atoms_per_solvent;
            let pos_j = (atoms[ow_j].pos.x, atoms[ow_j].pos.y, atoms[ow_j].pos.z);
            if dist_sq_pbc(pos_i, pos_j, &box_dim) < mindist_sq {
                excluded.insert(j);
            }
        }
    }

    // For negative ions: pick molecules with HIGHEST potential (most positive)
    for _ in 0..n_negative {
        let mut best_idx = None;
        let mut best_pot = f64::MIN;

        for j in 0..n_solvent_mols {
            if excluded.contains(&j) {
                continue;
            }
            if potentials[j] > best_pot {
                best_pot = potentials[j];
                best_idx = Some(j);
            }
        }

        let idx = best_idx.unwrap_or_else(|| {
            eprintln!("Error: not enough suitable solvent positions for negative ions");
            process::exit(1);
        });

        selected.push(idx);
        let ow_i = n_solute + idx * atoms_per_solvent;
        let pos_i = (atoms[ow_i].pos.x, atoms[ow_i].pos.y, atoms[ow_i].pos.z);

        for j in 0..n_solvent_mols {
            let ow_j = n_solute + j * atoms_per_solvent;
            let pos_j = (atoms[ow_j].pos.x, atoms[ow_j].pos.y, atoms[ow_j].pos.z);
            if dist_sq_pbc(pos_i, pos_j, &box_dim) < mindist_sq {
                excluded.insert(j);
            }
        }
    }

    // Report selected positions
    let replaced_set: HashSet<usize> = selected.iter().copied().collect();
    for (i, &mol_idx) in selected.iter().enumerate() {
        let ow = n_solute + mol_idx * atoms_per_solvent;
        let ion_type = if i < n_positive { &pos_name } else { &neg_name };
        eprintln!("Ion {}: {} at solvent mol {} ({:.3}, {:.3}, {:.3}), potential = {:.3} kJ/(mol*e)",
            i + 1, ion_type, mol_idx + 1,
            atoms[ow].pos.x, atoms[ow].pos.y, atoms[ow].pos.z, potentials[mol_idx]);
    }

    // Write output: TITLE + POSITION (solute + ions + remaining solvent) + BOX
    println!("TITLE");
    println!("    ion: replaced {} solvent molecules in {}", n_positive + n_negative, args.pos);
    if n_positive > 0 {
        println!("    {} positive ions ({})", n_positive, pos_name);
    }
    if n_negative > 0 {
        println!("    {} negative ions ({})", n_negative, neg_name);
    }
    println!("END");

    println!("POSITION");
    let mut atom_counter = 1;

    // Write solute atoms as-is
    for i in 0..n_solute {
        let a = &atoms[i];
        println!("{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
            a.res_num, a.res_name, a.atom_name, atom_counter,
            a.pos.x, a.pos.y, a.pos.z);
        atom_counter += 1;
    }

    // Write ion atoms (at the oxygen position of replaced water molecules)
    // Positive ions first, then negative
    let solute_res_max = if n_solute > 0 { atoms[n_solute - 1].res_num } else { 0 };
    let mut ion_res = solute_res_max;

    for (i, &mol_idx) in selected.iter().enumerate() {
        let ow = n_solute + mol_idx * atoms_per_solvent;
        let ion_name = if i < n_positive { &pos_name } else { &neg_name };
        ion_res += 1;
        println!("{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
            ion_res, ion_name, ion_name, atom_counter,
            atoms[ow].pos.x, atoms[ow].pos.y, atoms[ow].pos.z);
        atom_counter += 1;
    }

    // Write remaining solvent
    for mol_i in 0..n_solvent_mols {
        if replaced_set.contains(&mol_i) {
            continue;
        }
        for k in 0..atoms_per_solvent {
            let a = &atoms[n_solute + mol_i * atoms_per_solvent + k];
            println!("{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
                a.res_num, a.res_name, a.atom_name, atom_counter,
                a.pos.x, a.pos.y, a.pos.z);
            atom_counter += 1;
        }
    }

    println!("END");
    println!("BOX");
    println!("{:15.9}{:15.9}{:15.9}", box_dim[0], box_dim[1], box_dim[2]);
    println!("END");

    let final_solvent = n_solvent_mols - n_positive - n_negative;
    eprintln!("Output: {} solute + {} ions + {} solvent molecules ({} atoms)",
        n_solute, n_positive + n_negative, final_solvent,
        n_solute + n_positive + n_negative + final_solvent * atoms_per_solvent);
}
