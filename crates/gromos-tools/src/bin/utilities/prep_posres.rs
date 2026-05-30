//! prep_posres - Generate position restraint files (.por and .rpr)
//!
//! Usage: prep_posres @topo <topology> @pos <coordinates.cnf>
//!                    [@por <output.por>] [@rpr <output.rpr>]
//!
//! Creates two files from a solvated coordinate file:
//!   .por (POSRESSPEC) — lists solute atoms to be restrained
//!   .rpr (REFPOSITION) — reference positions for all atoms
//!
//! The number of solute atoms is determined from the topology.

use clap::Parser;
use gromos_io::gromos_args;
use gromos_io::topology::read_topology_file;
use gromos_io::{read_g96_labeled, write_por, write_rpr};
use std::process;

#[derive(Parser)]
#[command(name = "prep_posres", version, about = "Generate position restraint files (.por/.rpr)")]
struct Args {
    /// Molecular topology file
    #[arg(long)]
    topo: String,

    /// Input coordinate file (solvated system)
    #[arg(long)]
    pos: String,

    /// Output .por file (POSRESSPEC, solute atoms only)
    #[arg(long)]
    por: Option<String>,

    /// Output .rpr file (REFPOSITION, all atoms)
    #[arg(long)]
    rpr: Option<String>,
}

fn main() {
    let args = Args::parse_from(gromos_args());

    // Read topology to get number of solute atoms
    let topo = match read_topology_file(&args.topo) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        }
    };
    let n_solute = topo.n_atoms;

    // Read coordinate file
    let coord = match read_g96_labeled(&args.pos) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error reading coordinates: {}", e);
            process::exit(1);
        }
    };

    eprintln!("prep_posres:");
    eprintln!("  Topology:     {}", args.topo);
    eprintln!("  Coordinates:  {}", args.pos);
    eprintln!("  Total atoms:  {}", coord.atoms.len());
    eprintln!("  Solute atoms: {}", n_solute);

    if n_solute > coord.atoms.len() {
        eprintln!(
            "Error: topology has {} solute atoms but coordinate file has only {} atoms",
            n_solute,
            coord.atoms.len()
        );
        process::exit(1);
    }

    // Derive default names from input
    let base = args.pos.trim_end_matches(".cnf").trim_end_matches(".g96");

    // Write .por
    let por_path = args.por.unwrap_or_else(|| format!("{}.por", base));
    if let Err(e) = write_por(&por_path, &coord.atoms, n_solute) {
        eprintln!("Error writing .por: {}", e);
        process::exit(1);
    }
    eprintln!("  Wrote .por:   {} ({} solute atoms)", por_path, n_solute);

    // Write .rpr
    let rpr_path = args.rpr.unwrap_or_else(|| format!("{}.rpr", base));
    if let Err(e) = write_rpr(&rpr_path, &coord.atoms) {
        eprintln!("Error writing .rpr: {}", e);
        process::exit(1);
    }
    eprintln!("  Wrote .rpr:   {} ({} atoms)", rpr_path, coord.atoms.len());
}
