//! sim_box - Solvate a solute in a box of pre-equilibrated solvent
//!
//! Usage: sim_box @topo <topology> @pbc <r|t> @pos <solute.cnf> @solvent <solvent.g96>
//!               [@minwall <distance>] [@thresh <distance>] [@boxsize]
//!
//! Solvates a molecular system the way gromos++ `sim_box` does (`programs/sim_box.cc`):
//! 1. box from the longest solute atom–atom distance plus twice `@minwall` (cubic with one
//!    value; with three values the x/xy/xyz maximum distances, which assumes a solute rotated
//!    with its longest axis along z), or the solute file's box with `@boxsize`;
//! 2. solute shifted to its centre of geometry, the solvent template centred likewise and
//!    replicated `int(box / template) + 1` times per dimension around the origin;
//! 3. a solvent molecule is kept when its centre of geometry lies inside the box and is farther
//!    than `@thresh` from every non-hydrogen solute atom (hydrogens: mass 1.008 in `@topo`);
//! 4. the solvated system is written to stdout as GROMOS96 with a GENBOX block.
//!
//! Reference: `tests/sim_box_reference.rs` reproduces gromos++'s output for methanol in SPC.

use gromos::math::Vec3;
use gromos_io::{read_g96_labeled, G96Atom};
use std::env;
use std::process;

fn print_usage() {
    eprintln!("sim_box - Solvate a solute in a box of pre-equilibrated solvent");
    eprintln!();
    eprintln!(
        "Usage: sim_box @topo <topology> @pbc <r|t> @pos <solute.cnf> @solvent <solvent.g96> \\"
    );
    eprintln!("               [@minwall <distance>] [@thresh <distance>] [@boxsize]");
    eprintln!();
    eprintln!("Arguments:");
    eprintln!("  @topo      Molecular topology file");
    eprintln!("  @pbc       Periodic boundary conditions:");
    eprintln!("             r = rectangular box");
    eprintln!("             t = truncated octahedron (not yet implemented)");
    eprintln!("  @pos       Input coordinate file for the solute");
    eprintln!("  @solvent   Input coordinate file for pre-equilibrated solvent");
    eprintln!("  @minwall   Minimum solute-to-wall distance (nm)");
    eprintln!("             - 1 value: cubic box");
    eprintln!("             - 3 values: rectangular box (x y z)");
    eprintln!("  @thresh    Minimum solvent-to-solute distance (default: 0.23 nm)");
    eprintln!("  @boxsize   Use box dimensions from solute coordinate file");
    eprintln!();
    eprintln!("Description:");
    eprintln!("  Solvates a solute molecule in a pre-equilibrated solvent box.");
    eprintln!("  The solvent configuration should contain a BOX block.");
    eprintln!();
    eprintln!("  Box sizing options:");
    eprintln!("  1. @boxsize: Use box from solute file");
    eprintln!("  2. @minwall: Calculate box from solute extent + minwall distance");
    eprintln!();
    eprintln!("  Solvent removal criteria:");
    eprintln!("  - Solvent molecule center of geometry must be:");
    eprintln!("    a) Inside the target box");
    eprintln!("    b) At least @thresh distance from any solute atom");
    eprintln!();
    eprintln!("Examples:");
    eprintln!("  # Cubic box with 1.4 nm minimum wall distance");
    eprintln!("  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @minwall 1.4");
    eprintln!();
    eprintln!("  # Use box size from solute file");
    eprintln!("  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @boxsize");
    eprintln!();
    eprintln!("  # Rectangular box with different wall distances");
    eprintln!(
        "  sim_box @topo system.top @pbc r @pos solute.cnf @solvent h2o.g96 @minwall 1.0 1.5 2.0"
    );
    eprintln!();
    eprintln!("Output:");
    eprintln!("  GROMOS96 format to stdout containing solvated system");
}

#[derive(Debug)]
struct SimBoxArgs {
    topo: String,
    pbc: char,
    pos: String,
    solvent: String,
    minwall: Vec<f64>,
    thresh: f64,
    use_boxsize: bool,
}

impl Default for SimBoxArgs {
    fn default() -> Self {
        Self {
            topo: String::new(),
            pbc: 'r',
            pos: String::new(),
            solvent: String::new(),
            minwall: Vec::new(),
            thresh: 0.23, // Default from GROMOS++
            use_boxsize: false,
        }
    }
}

fn parse_args(args: Vec<String>) -> Result<SimBoxArgs, String> {
    let mut sb_args = SimBoxArgs::default();

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing topology file for @topo".to_string());
                }
                sb_args.topo = args[i].clone();
            },
            "@pbc" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing boundary type for @pbc".to_string());
                }
                let pbc_str = args[i].clone();
                if pbc_str.len() != 1 {
                    return Err(format!("Invalid @pbc: {}", pbc_str));
                }
                sb_args.pbc = pbc_str.chars().next().unwrap();
                if sb_args.pbc != 'r' && sb_args.pbc != 't' {
                    return Err(format!(
                        "Invalid @pbc: {} (must be 'r' or 't')",
                        sb_args.pbc
                    ));
                }
            },
            "@pos" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing solute file for @pos".to_string());
                }
                sb_args.pos = args[i].clone();
            },
            "@solvent" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing solvent file for @solvent".to_string());
                }
                sb_args.solvent = args[i].clone();
            },
            "@minwall" => {
                // Read 1 or 3 values
                loop {
                    i += 1;
                    if i >= args.len() || args[i].starts_with('@') {
                        i -= 1;
                        break;
                    }
                    let val: f64 = args[i]
                        .parse()
                        .map_err(|_| format!("Invalid @minwall value: {}", args[i]))?;
                    sb_args.minwall.push(val);
                }
                if sb_args.minwall.len() != 1 && sb_args.minwall.len() != 3 {
                    return Err(format!(
                        "@minwall requires 1 or 3 values, got {}",
                        sb_args.minwall.len()
                    ));
                }
            },
            "@thresh" => {
                i += 1;
                if i >= args.len() {
                    return Err("Missing threshold for @thresh".to_string());
                }
                sb_args.thresh = args[i]
                    .parse()
                    .map_err(|_| format!("Invalid @thresh: {}", args[i]))?;
            },
            "@boxsize" => {
                sb_args.use_boxsize = true;
            },
            _ => {
                return Err(format!("Unknown argument: {}", args[i]));
            },
        }
        i += 1;
    }

    // Validate
    if sb_args.topo.is_empty() {
        return Err("Missing required argument: @topo".to_string());
    }
    if sb_args.pos.is_empty() {
        return Err("Missing required argument: @pos".to_string());
    }
    if sb_args.solvent.is_empty() {
        return Err("Missing required argument: @solvent".to_string());
    }
    if !sb_args.use_boxsize && sb_args.minwall.is_empty() {
        return Err("Either @boxsize or @minwall must be specified".to_string());
    }
    if sb_args.use_boxsize && !sb_args.minwall.is_empty() {
        return Err("Cannot specify both @boxsize and @minwall".to_string());
    }

    Ok(sb_args)
}

/// Calculate center of geometry
fn calc_cog(atoms: &[G96Atom]) -> Vec3 {
    let mut cog = Vec3::ZERO;
    for a in atoms {
        cog += a.pos;
    }
    cog / atoms.len() as f64
}

/// Shift all positions by a vector
fn shift_atoms(atoms: &mut [G96Atom], shift: Vec3) {
    for a in atoms {
        a.pos += shift;
    }
}

/// gromos++ `calc_max_size`: the longest distance between two solute atoms considering the
/// first `dim` coordinates (3 = full distance, 2 = in the xy plane, 1 = along x).
fn calc_max_size(atoms: &[G96Atom], dim: usize) -> (f64, usize, usize) {
    let mut max2 = 0.0;
    let mut pair = (0, 0);
    for (i, a) in atoms.iter().enumerate() {
        for (j, b) in atoms.iter().enumerate().skip(i) {
            let d = [a.pos.x - b.pos.x, a.pos.y - b.pos.y, a.pos.z - b.pos.z];
            let d2: f64 = d.iter().take(dim).map(|x| x * x).sum();
            if d2 > max2 {
                max2 = d2;
                pair = (i, j);
            }
        }
    }
    (max2.sqrt(), pair.0, pair.1)
}

/// gromos++ hydrogen rule (`setHmass(1.008)` + `isH()`): the topology mass is exactly 1.008.
fn is_hydrogen(mass: f64) -> bool {
    (mass - 1.008).abs() < 1e-9
}

/// Rectangular nearest image of `r` to the origin, as gromos++ `nearestImage` computes it
/// (`r - box * rint(r / box)`); a molecule is inside the box when this returns `r` unchanged.
fn nearest_image(r: Vec3, box_dims: Vec3) -> Vec3 {
    Vec3::new(
        r.x - box_dims.x * (r.x / box_dims.x).round(),
        r.y - box_dims.y * (r.y / box_dims.y).round(),
        r.z - box_dims.z * (r.z / box_dims.z).round(),
    )
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 || args[1] == "-h" || args[1] == "--help" {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let sb_args = match parse_args(args) {
        Ok(args) => args,
        Err(e) => {
            eprintln!("Error: {}", e);
            eprintln!();
            print_usage();
            process::exit(1);
        },
    };

    if sb_args.pbc == 't' {
        eprintln!("Error: Truncated octahedron (@pbc t) not yet implemented");
        eprintln!("       Use rectangular boxes (@pbc r) for now");
        process::exit(1);
    }

    eprintln!("sim_box - Solvating system");
    eprintln!("  Topology:   {}", sb_args.topo);
    eprintln!("  PBC:        {}", sb_args.pbc);
    eprintln!("  Solute:     {}", sb_args.pos);
    eprintln!("  Solvent:    {}", sb_args.solvent);
    eprintln!("  Threshold:  {} nm", sb_args.thresh);
    eprintln!();

    // Read solute coordinates
    eprintln!("Reading solute coordinates...");
    let (mut solute_atoms, solute_box) = match read_g96_labeled(&sb_args.pos) {
        Ok(data) => {
            let box_opt = if data.box_dims == Vec3::ZERO {
                None
            } else {
                Some(data.box_dims)
            };
            (data.atoms, box_opt)
        },
        Err(e) => {
            eprintln!("Error reading solute: {}", e);
            process::exit(1);
        },
    };
    eprintln!("  Solute: {} atoms", solute_atoms.len());
    // Solute masses: the clash filter ignores hydrogens, as gromos++ does.
    let solute_is_h: Vec<bool> = match gromos_io::topology::read_topology_file(&sb_args.topo) {
        Ok(parsed) => {
            if parsed.masses.len() < solute_atoms.len() {
                eprintln!(
                    "Error: topology has {} solute atoms, coordinate file {}",
                    parsed.masses.len(),
                    solute_atoms.len()
                );
                process::exit(1);
            }
            parsed.masses[..solute_atoms.len()]
                .iter()
                .map(|&m| is_hydrogen(m))
                .collect()
        },
        Err(e) => {
            eprintln!("Error reading topology: {}", e);
            process::exit(1);
        },
    };

    // Read solvent coordinates
    eprintln!("Reading solvent coordinates...");
    let (solvent_atoms_orig, solvent_box_opt) = match read_g96_labeled(&sb_args.solvent) {
        Ok(data) => {
            let box_opt = if data.box_dims == Vec3::ZERO {
                None
            } else {
                Some(data.box_dims)
            };
            (data.atoms, box_opt)
        },
        Err(e) => {
            eprintln!("Error reading solvent: {}", e);
            process::exit(1);
        },
    };

    let solvent_box = match solvent_box_opt {
        Some(b) => b,
        None => {
            eprintln!("Error: Solvent file must contain a BOX block");
            process::exit(1);
        },
    };
    eprintln!(
        "  Solvent: {} atoms in box ({:.3}, {:.3}, {:.3}) nm",
        solvent_atoms_orig.len(),
        solvent_box.x,
        solvent_box.y,
        solvent_box.z
    );

    // Determine target box size
    let target_box: Vec3 = if sb_args.use_boxsize {
        // Use box from solute file
        match solute_box {
            Some(b) => {
                eprintln!(
                    "Using box size from solute file: ({:.3}, {:.3}, {:.3}) nm",
                    b.x, b.y, b.z
                );
                b
            },
            None => {
                eprintln!("Error: @boxsize specified but solute file has no BOX block");
                process::exit(1);
            },
        }
    } else {
        // gromos++: box = longest solute atom-atom distance + 2 * minwall
        if sb_args.minwall.len() == 1 {
            let (max_dist, a, b) = calc_max_size(&solute_atoms, 3);
            let size = max_dist + 2.0 * sb_args.minwall[0];
            eprintln!(
                "  Cubic box {:.6} nm: maximum solute atom-atom distance {:.6} nm (atoms {} and {}) + 2 x minwall {:.3} nm",
                size,
                max_dist,
                a + 1,
                b + 1,
                sb_args.minwall[0]
            );
            Vec3::new(size, size, size)
        } else {
            // three minwall values: gromos++ expects the solute rotated (longest axis along z)
            // and takes the x, xy and xyz maximum distances for K, L and M respectively
            let (dx, _, _) = calc_max_size(&solute_atoms, 1);
            let (dxy, _, _) = calc_max_size(&solute_atoms, 2);
            let (dxyz, _, _) = calc_max_size(&solute_atoms, 3);
            let size = Vec3::new(
                dx + 2.0 * sb_args.minwall[0],
                dxy + 2.0 * sb_args.minwall[1],
                dxyz + 2.0 * sb_args.minwall[2],
            );
            eprintln!(
                "  Rectangular box ({:.6}, {:.6}, {:.6}) nm from solute distances ({:.6}, {:.6}, {:.6}) nm along x / in xy / in xyz",
                size.x, size.y, size.z, dx, dxy, dxyz
            );
            eprintln!("  (gromos++ assumes a solute rotated with its longest axis along z)");
            size
        }
    };
    // Move solute to center
    let solute_cog = calc_cog(&solute_atoms);
    let shift_to_origin = Vec3::ZERO - solute_cog;
    shift_atoms(&mut solute_atoms, shift_to_origin);
    eprintln!("  Shifted solute to origin");

    // Centre the solvent template on its centre of geometry (gromos++ does the same).
    let mut solvent_atoms_orig = solvent_atoms_orig;
    let solvent_cog = calc_cog(&solvent_atoms_orig);
    shift_atoms(&mut solvent_atoms_orig, Vec3::ZERO - solvent_cog);

    // Number of template copies per dimension: gromos++ `int(box / template) + 1`, laid out
    // symmetrically around the origin.
    let needed = [
        (target_box.x / solvent_box.x) as i32 + 1,
        (target_box.y / solvent_box.y) as i32 + 1,
        (target_box.z / solvent_box.z) as i32 + 1,
    ];
    let start = Vec3::new(
        -0.5 * (needed[0] - 1) as f64 * solvent_box.x,
        -0.5 * (needed[1] - 1) as f64 * solvent_box.y,
        -0.5 * (needed[2] - 1) as f64 * solvent_box.z,
    );
    eprintln!();
    eprintln!("Replicating solvent box...");
    eprintln!(
        "  {} x {} x {} = {} copies of the template",
        needed[0],
        needed[1],
        needed[2],
        needed[0] * needed[1] * needed[2]
    );
    let mut all_solvent: Vec<G96Atom> = Vec::new();
    for ix in 0..needed[0] {
        for iy in 0..needed[1] {
            for iz in 0..needed[2] {
                let shift = start
                    + Vec3::new(
                        ix as f64 * solvent_box.x,
                        iy as f64 * solvent_box.y,
                        iz as f64 * solvent_box.z,
                    );
                for atom in &solvent_atoms_orig {
                    let mut shifted = atom.clone();
                    shifted.pos = atom.pos + shift;
                    all_solvent.push(shifted);
                }
            }
        }
    }
    eprintln!(
        "  Total solvent atoms before filtering: {}",
        all_solvent.len()
    );

    // Atoms per solvent molecule from the template's residue numbering (3 for water).
    let atoms_per_molecule = {
        let first = solvent_atoms_orig.first().map(|a| a.res_num);
        let n = solvent_atoms_orig
            .iter()
            .take_while(|a| Some(a.res_num) == first)
            .count();
        n.max(1)
    };
    let num_molecules = all_solvent.len() / atoms_per_molecule;

    // Keep a molecule when its centre of geometry is inside the box (nearest image to the
    // origin is itself) and farther than @thresh from every non-hydrogen solute atom.
    eprintln!("Filtering solvent molecules...");
    let mut kept_solvent: Vec<G96Atom> = Vec::new();
    let thresh2 = sb_args.thresh * sb_args.thresh;
    let min_init = target_box.length_squared();
    for mol_idx in 0..num_molecules {
        let mol_atoms =
            &all_solvent[mol_idx * atoms_per_molecule..(mol_idx + 1) * atoms_per_molecule];
        let cog = calc_cog(mol_atoms);
        let check = nearest_image(cog, target_box);
        if check != cog {
            continue; // outside the box
        }
        let mut min2 = min_init;
        for (atom, &is_h) in solute_atoms.iter().zip(&solute_is_h) {
            if is_h {
                continue;
            }
            let d2 = (check - atom.pos).length_squared();
            if d2 < min2 {
                min2 = d2;
            }
        }
        if min2 > thresh2 {
            kept_solvent.extend(mol_atoms.iter().cloned());
        }
    }
    let kept_molecules = kept_solvent.len() / atoms_per_molecule;
    eprintln!(
        "  Kept {} solvent molecules ({} atoms), removed {}",
        kept_molecules,
        kept_solvent.len(),
        num_molecules - kept_molecules
    );

    eprintln!();
    eprintln!("Final system:");
    eprintln!(
        "  {} solute atoms + {} solvent atoms = {}",
        solute_atoms.len(),
        kept_solvent.len(),
        solute_atoms.len() + kept_solvent.len()
    );
    eprintln!(
        "  Box: ({:.6}, {:.6}, {:.6}) nm",
        target_box.x, target_box.y, target_box.z
    );
    eprintln!();

    // Write to stdout (gromos++ OutG96S layout: solute as read, solvent residues SOLV numbered
    // per molecule from 1, GENBOX with the box angles, Euler angles and origin).
    println!("TITLE");
    println!("Solvating {} in {}", sb_args.pos, sb_args.solvent);
    println!("Added {} solvent molecules", kept_molecules);
    println!("END");
    println!("POSITION");
    let mut serial = 0usize;
    for atom in &solute_atoms {
        serial += 1;
        println!(
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}",
            atom.res_num, atom.res_name, atom.atom_name, serial, atom.pos.x, atom.pos.y, atom.pos.z
        );
    }
    for (i, atom) in kept_solvent.iter().enumerate() {
        serial += 1;
        println!(
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}",
            i / atoms_per_molecule + 1,
            "SOLV",
            atom.atom_name,
            serial,
            atom.pos.x,
            atom.pos.y,
            atom.pos.z
        );
    }
    println!("END");
    println!("GENBOX");
    println!("{:>8}", 1);
    println!(
        "{:15.9}{:15.9}{:15.9}",
        target_box.x, target_box.y, target_box.z
    );
    println!("{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0);
    println!("{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0);
    println!("{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0);
    println!("END");
}
