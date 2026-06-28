//! build_box - Generate a condensed-phase system on a regular grid
//!
//! Usage: build_box @topo <topology> @pos <molecule.g96> @nsm <N> @dens <density>
//!
//! Replicates a single molecule on an N x N x N grid into a cubic box whose
//! size is chosen so that the resulting system matches the target density.

use clap::Parser;
use gromos::math::Vec3;
use gromos_io::topology::read_topology_file;
use gromos_io::{gromos_args, read_g96_labeled, G96Atom};
use std::process;

/// Conversion factor from atomic mass units to kg per nm^3 at unit density
/// (gromos-rs `build_box.cc`: `vtot = nsm * (weight * 1.66056) / density`).
const AMU_TO_KG_PER_NM3: f64 = 1.66056;

#[derive(Parser)]
#[command(
    name = "build_box",
    version,
    about = "Generate a condensed-phase system on a grid"
)]
struct Args {
    /// Molecular topology file for a single molecule
    #[arg(long)]
    topo: String,

    /// Input coordinate file for a single molecule
    #[arg(long)]
    pos: String,

    /// Number of molecules per dimension (total = nsm^3)
    #[arg(long)]
    nsm: usize,

    /// Density of the liquid (kg/m^3)
    #[arg(long)]
    dens: f64,
}

/// Computes (box_side, cell_side) for `nsm` copies per dimension of a molecule
/// of the given `weight` (u) so that the system matches the target `density`
/// (kg/m^3). Mirrors gromos-rs `build_box.cc`: `vtot = nsm^3 * (weight * 1.66056) / density`.
fn box_dimensions(nsm: usize, weight: f64, density: f64) -> (f64, f64) {
    let nsm_total = nsm * nsm * nsm;
    let vtot = nsm_total as f64 * (weight * AMU_TO_KG_PER_NM3) / density;
    let box_side = vtot.cbrt();
    let cell = box_side / nsm as f64;
    (box_side, cell)
}

/// Mass-weighted center of geometry (center of mass)
fn center_of_mass(atoms: &[G96Atom], masses: &[f64]) -> Vec3 {
    let mut com = Vec3::ZERO;
    let mut total_mass = 0.0;
    for (atom, &mass) in atoms.iter().zip(masses) {
        com = com + atom.pos * mass;
        total_mass += mass;
    }
    com / total_mass
}

fn main() {
    let args = Args::parse_from(gromos_args());

    if args.nsm == 0 {
        eprintln!("Error: @nsm must be a positive integer");
        process::exit(1);
    }

    // Read topology to get the molecule's total mass
    let parsed_topo = match read_topology_file(&args.topo) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("Error reading topology '{}': {:?}", args.topo, e);
            process::exit(1);
        },
    };
    let weight: f64 = parsed_topo.masses.iter().sum();

    // Read the single-molecule coordinates
    let mol_atoms = match read_g96_labeled(&args.pos) {
        Ok(data) => data.atoms,
        Err(e) => {
            eprintln!("Error reading '{}': {}", args.pos, e);
            process::exit(1);
        },
    };
    if mol_atoms.len() != parsed_topo.masses.len() {
        eprintln!(
            "Error: topology has {} atoms but '{}' has {}",
            parsed_topo.masses.len(),
            args.pos,
            mol_atoms.len()
        );
        process::exit(1);
    }

    // Compute box size from the requested number of copies and target density
    let nsm_total = args.nsm * args.nsm * args.nsm;
    let (box_side, cell) = box_dimensions(args.nsm, weight, args.dens);

    eprintln!("build_box - Building condensed-phase system");
    eprintln!("  Topology:        {}", args.topo);
    eprintln!("  Molecule:        {}", args.pos);
    eprintln!("  Copies per dim:  {}", args.nsm);
    eprintln!("  Total copies:    {}", nsm_total);
    eprintln!("  Molecular mass:  {:.4} u", weight);
    eprintln!("  Target density:  {} kg/m^3", args.dens);
    eprintln!("  Box side:        {:.6} nm", box_side);
    eprintln!("  Cell side:       {:.6} nm", cell);

    // Center the molecule's center of mass in the first grid cell
    let cell_center = Vec3::new(cell, cell, cell) * 0.5;
    let com = center_of_mass(&mol_atoms, &parsed_topo.masses);
    let mut centered_atoms = mol_atoms.clone();
    for atom in &mut centered_atoms {
        atom.pos = atom.pos - com + cell_center;
    }

    // Replicate on an nsm x nsm x nsm grid
    let mut all_atoms: Vec<G96Atom> = Vec::with_capacity(centered_atoms.len() * nsm_total);
    let mut atom_num = 1usize;
    let mut res_num = 1usize;
    for i in 0..args.nsm {
        for j in 0..args.nsm {
            for k in 0..args.nsm {
                let shift = Vec3::new(i as f64 * cell, j as f64 * cell, k as f64 * cell);
                for atom in &centered_atoms {
                    all_atoms.push(G96Atom {
                        res_num,
                        res_name: atom.res_name.clone(),
                        atom_name: atom.atom_name.clone(),
                        atom_num,
                        pos: atom.pos + shift,
                    });
                    atom_num += 1;
                }
                res_num += 1;
            }
        }
    }

    eprintln!("  Total atoms:     {}", all_atoms.len());
    eprintln!();

    // Write to stdout in GROMOS96 format
    println!("TITLE");
    println!(
        "build_box: {} copies of {}\nDensity : {} kg/m^3\tMolecular weight : {:.4} u",
        nsm_total, args.pos, args.dens, weight
    );
    println!("END");

    println!("POSITION");
    for atom in &all_atoms {
        println!(
            "{:>5} {:5} {:>5}{:7}{:15.9}{:15.9}{:15.9}",
            atom.res_num,
            atom.res_name,
            atom.atom_name,
            atom.atom_num,
            atom.pos.x,
            atom.pos.y,
            atom.pos.z
        );
    }
    println!("END");

    println!("BOX");
    println!("{:15.9}{:15.9}{:15.9}", box_side, box_side, box_side);
    println!("END");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_dimensions_match_target_density() {
        // 27 SPC waters (M = 18.0154 u) at 1000 kg/m^3
        let (box_side, cell) = box_dimensions(3, 18.0154, 1000.0);
        let nsm_total = 27.0;
        let vtot = nsm_total * (18.0154 * AMU_TO_KG_PER_NM3) / 1000.0;
        assert!((box_side - vtot.cbrt()).abs() < 1e-12);
        assert!((cell * 3.0 - box_side).abs() < 1e-12);
    }

    #[test]
    fn center_of_mass_is_mass_weighted() {
        let masses = vec![16.0, 1.0, 1.0];
        let atoms = vec![
            G96Atom {
                res_num: 1,
                res_name: "SOL".into(),
                atom_name: "OW".into(),
                atom_num: 1,
                pos: Vec3::new(0.0, 0.0, 0.0),
            },
            G96Atom {
                res_num: 1,
                res_name: "SOL".into(),
                atom_name: "HW1".into(),
                atom_num: 2,
                pos: Vec3::new(1.0, 0.0, 0.0),
            },
            G96Atom {
                res_num: 1,
                res_name: "SOL".into(),
                atom_name: "HW2".into(),
                atom_num: 3,
                pos: Vec3::new(0.0, 1.0, 0.0),
            },
        ];
        let com = center_of_mass(&atoms, &masses);
        let total_mass: f64 = masses.iter().sum();
        assert!((com.x - 1.0 / total_mass).abs() < 1e-12);
        assert!((com.y - 1.0 / total_mass).abs() < 1e-12);
        assert!((com.z - 0.0).abs() < 1e-12);
    }
}
