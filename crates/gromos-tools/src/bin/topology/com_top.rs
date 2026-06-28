//! com_top - Combine multiple GROMOS topology files
//!
//! Usage: com_top @topo <file1> [N:file2] ... @param <n> @solv <n>
//!    or: com_top @f <argfile>
//!
//! Combines multiple topology files into a single topology.
//! Use prefix 'N:' before a file to repeat that topology N times.
//! Parameters and solvent are taken from the topology indicated by @param/@solv.

use clap::Parser;
use gromos_core::topology::{
    AngleParameters, BondParameters, DihedralParameters, ImproperDihedralParameters, LJParameters,
};
use gromos_io::gromos_args;
use gromos_io::topology::{
    read_topology_file, ParsedSolventAtom, ParsedSolventConstraint, ParsedTopology,
};
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
    write_combined(&combined, &title);
}

/// Combined topology ready for output
struct CombinedTopology {
    atom_names: Vec<String>,
    residue_numbers: Vec<usize>,
    residue_names: Vec<String>,
    atom_type_names: Vec<String>,
    masses: Vec<f64>,
    charges: Vec<f64>,
    iac: Vec<usize>,
    chargegroup_codes: Vec<usize>,
    exclusions: Vec<Vec<usize>>,
    one_four_pairs: Vec<Vec<usize>>,
    bonds: Vec<(usize, usize, usize)>,
    bond_parameters: Vec<BondParameters>,
    angles: Vec<(usize, usize, usize, usize)>,
    angle_parameters: Vec<AngleParameters>,
    proper_dihedrals: Vec<(usize, usize, usize, usize, usize)>,
    dihedral_parameters: Vec<DihedralParameters>,
    improper_dihedrals: Vec<(usize, usize, usize, usize, usize)>,
    improper_dihedral_parameters: Vec<ImproperDihedralParameters>,
    lj_parameters: HashMap<(usize, usize), LJParameters>,
    temperature_groups: Vec<usize>,
    pressure_groups: Vec<usize>,
    solvent_atoms: Vec<ParsedSolventAtom>,
    solvent_constraints: Vec<ParsedSolventConstraint>,
    solute_molecules: Vec<usize>,
}

fn combine_topologies(
    specs: &[(usize, String)],
    file_cache: &HashMap<String, ParsedTopology>,
    param_file: &str,
    solv_file: &str,
) -> CombinedTopology {
    let param_topo = &file_cache[param_file];
    let solv_topo = &file_cache[solv_file];

    let mut combined = CombinedTopology {
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

    eprintln!(
        "Combined topology: {} atoms, {} residues, {} molecules",
        combined.masses.len(),
        combined.residue_names.len(),
        combined.solute_molecules.len()
    );

    combined
}

fn write_combined(topo: &CombinedTopology, title: &str) {
    let n_atoms = topo.masses.len();

    // TITLE
    println!("TITLE");
    print!("{}", title);
    println!("END");

    // PHYSICALCONSTANTS
    println!("PHYSICALCONSTANTS");
    println!("# FPEPSI");
    println!("{:15.7}", gromos_core::units::four_pi_eps_i);
    println!("# HBAR");
    println!("{:15.7}", gromos_core::units::hBar);
    println!("# SPDL");
    println!("{:15.3}", gromos_core::units::spd_l);
    println!("# BOLTZ");
    println!("{:15.8}", gromos_core::units::kB);
    println!("END");

    // TOPVERSION
    println!("TOPVERSION");
    println!("2.0");
    println!("END");

    // ATOMTYPENAME
    println!("ATOMTYPENAME");
    println!("# NRATT: number of atom types");
    println!("{:5}", topo.atom_type_names.len());
    println!("# TYPE: atom type names");
    for name in &topo.atom_type_names {
        println!("{}", name);
    }
    println!("END");

    // RESNAME
    println!("RESNAME");
    println!("# NRAA2: number of residue names");
    println!("{:5}", topo.residue_names.len());
    println!("# AANM: residue names");
    for name in &topo.residue_names {
        println!("{}", name);
    }
    println!("END");

    // SOLUTEATOM
    println!("SOLUTEATOM");
    println!("#   NRP: number of solute atoms");
    println!("{:6}", n_atoms);
    println!("#  ATNM: atom number");
    println!("#  MRES: residue number");
    println!("#  PANM: atom name");
    println!("#   IAC: integer atom code");
    println!("#  MASS: mass");
    println!("#    CG: charge");
    println!("#   CGC: charge group code");
    println!("#   INE: number of excluded atoms");
    println!("# INE14: number of 1-4 interactions");
    println!("#  ATNM  MRES  PANM   IAC      MASS        CG  CGC   INE");
    println!("#                                                     INE14");

    for i in 0..n_atoms {
        let ne = topo.exclusions.get(i).map_or(0, |e| e.len());
        println!(
            "{:6}{:5} {:>5}{:5}{:10.5}{:12.5}{:3}{:4}",
            i + 1,
            topo.residue_numbers[i],
            topo.atom_names[i],
            topo.iac[i],
            topo.masses[i],
            topo.charges[i],
            topo.chargegroup_codes[i],
            ne
        );
        // Print exclusion list
        if let Some(excls) = topo.exclusions.get(i) {
            if !excls.is_empty() {
                let mut sorted: Vec<usize> = excls.clone();
                sorted.sort();
                for (j, &e) in sorted.iter().enumerate() {
                    print!("{:5}", e + 1);
                    if (j + 1) % 20 == 0 {
                        println!();
                    }
                }
                if sorted.len() % 20 != 0 {
                    println!();
                }
            }
        }
        // Print 1-4 pairs
        let n14 = topo.one_four_pairs.get(i).map_or(0, |p| p.len());
        println!("   {:4}", n14);
        if let Some(pairs) = topo.one_four_pairs.get(i) {
            if !pairs.is_empty() {
                let mut sorted: Vec<usize> = pairs.clone();
                sorted.sort();
                for (j, &p) in sorted.iter().enumerate() {
                    print!("{:5}", p + 1);
                    if (j + 1) % 20 == 0 {
                        println!();
                    }
                }
                if sorted.len() % 20 != 0 {
                    println!();
                }
            }
        }
    }
    println!("END");

    // BONDSTRETCHTYPE
    println!("BONDSTRETCHTYPE");
    println!("# NBTY: number of bond types");
    println!("{:5}", topo.bond_parameters.len());
    println!("#  CB     HB     B0");
    for bp in &topo.bond_parameters {
        println!(
            "{:15.7e}{:15.7e}{:13.7e}",
            bp.k_quartic, bp.k_harmonic, bp.r0
        );
    }
    println!("END");

    // Separate bonds into BONDH/BOND by H involvement
    let (bondh, bond): (Vec<_>, Vec<_>) = topo
        .bonds
        .iter()
        .partition(|&&(i, j, _)| topo.masses[i] < 2.0 || topo.masses[j] < 2.0);

    println!("BONDH");
    println!("# NBAH: number of bonds involving H atoms");
    println!("{:5}", bondh.len());
    println!("#  IB   JB  ICB");
    for &&(i, j, t) in &bondh {
        println!("{:5}{:5}{:5}", i + 1, j + 1, t + 1);
    }
    println!("END");

    println!("BOND");
    println!("# NBA: number of bonds not involving H atoms");
    println!("{:5}", bond.len());
    println!("#  IB   JB  ICB");
    for &&(i, j, t) in &bond {
        println!("{:5}{:5}{:5}", i + 1, j + 1, t + 1);
    }
    println!("END");

    // BONDANGLEBENDTYPE
    println!("BONDANGLEBENDTYPE");
    println!("# NTTY: number of angle types");
    println!("{:5}", topo.angle_parameters.len());
    println!("#   CT     CHT     T0");
    for ap in &topo.angle_parameters {
        println!(
            "{:15.7e}{:15.7e}{:13.7e}",
            ap.k_cosine, ap.k_harmonic, ap.theta0
        );
    }
    println!("END");

    // Separate angles
    let (angleh, angle): (Vec<_>, Vec<_>) = topo.angles.iter().partition(|&&(i, j, k, _)| {
        topo.masses[i] < 2.0 || topo.masses[j] < 2.0 || topo.masses[k] < 2.0
    });

    println!("BONDANGLEH");
    println!("# NTHEH: number of angles involving H atoms");
    println!("{:5}", angleh.len());
    println!("#  IT   JT   KT  ICT");
    for &&(i, j, k, t) in &angleh {
        println!("{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, t + 1);
    }
    println!("END");

    println!("BONDANGLE");
    println!("# NTHE: number of angles not involving H atoms");
    println!("{:5}", angle.len());
    println!("#  IT   JT   KT  ICT");
    for &&(i, j, k, t) in &angle {
        println!("{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, t + 1);
    }
    println!("END");

    // IMPDIHEDRALTYPECODE
    println!("IMPDIHEDRALTYPECODE");
    println!("# NQTY: number of improper dihedral types");
    println!("{:5}", topo.improper_dihedral_parameters.len());
    println!("#  CQ     Q0");
    for ip in &topo.improper_dihedral_parameters {
        println!("{:15.7e}{:13.4}", ip.k, ip.q0);
    }
    println!("END");

    // Separate impropers
    let (impdihedralh, impdihedral): (Vec<_>, Vec<_>) =
        topo.improper_dihedrals
            .iter()
            .partition(|&&(i, j, k, l, _)| {
                topo.masses[i] < 2.0
                    || topo.masses[j] < 2.0
                    || topo.masses[k] < 2.0
                    || topo.masses[l] < 2.0
            });

    println!("IMPDIHEDRALH");
    println!("# NQHIH: number of improper dihedrals involving H");
    println!("{:5}", impdihedralh.len());
    println!("#  IQ   JQ   KQ   LQ  ICQ");
    for &&(i, j, k, l, t) in &impdihedralh {
        println!("{:5}{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, l + 1, t + 1);
    }
    println!("END");

    println!("IMPDIHEDRAL");
    println!("# NQHI: number of improper dihedrals not involving H");
    println!("{:5}", impdihedral.len());
    println!("#  IQ   JQ   KQ   LQ  ICQ");
    for &&(i, j, k, l, t) in &impdihedral {
        println!("{:5}{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, l + 1, t + 1);
    }
    println!("END");

    // TORSDIHEDRALTYPECODE
    println!("TORSDIHEDRALTYPECODE");
    println!("# NPTY: number of dihedral types");
    println!("{:5}", topo.dihedral_parameters.len());
    println!("#   CP     PD    NP");
    for dp in &topo.dihedral_parameters {
        println!("{:15.7e}{:10.1}{:5}", dp.k, dp.pd, dp.m);
    }
    println!("END");

    // Separate dihedrals
    let (dihedralh, dihedral): (Vec<_>, Vec<_>) =
        topo.proper_dihedrals.iter().partition(|&&(i, j, k, l, _)| {
            topo.masses[i] < 2.0
                || topo.masses[j] < 2.0
                || topo.masses[k] < 2.0
                || topo.masses[l] < 2.0
        });

    println!("DIHEDRALH");
    println!("# NPHIH: number of dihedrals involving H");
    println!("{:5}", dihedralh.len());
    println!("#  IP   JP   KP   LP  ICP");
    for &&(i, j, k, l, t) in &dihedralh {
        println!("{:5}{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, l + 1, t + 1);
    }
    println!("END");

    println!("DIHEDRAL");
    println!("# NPHI: number of dihedrals not involving H");
    println!("{:5}", dihedral.len());
    println!("#  IP   JP   KP   LP  ICP");
    for &&(i, j, k, l, t) in &dihedral {
        println!("{:5}{:5}{:5}{:5}{:5}", i + 1, j + 1, k + 1, l + 1, t + 1);
    }
    println!("END");

    // LJPARAMETERS
    if !topo.lj_parameters.is_empty() {
        let n_types = topo.atom_type_names.len();
        let n_pairs = n_types * (n_types + 1) / 2;
        println!("LJPARAMETERS");
        println!("# NRATT2: number of LJ parameter pairs");
        println!("{:5}", n_pairs);
        println!("#  IAC  JAC    C12          C6           CS12         CS6");
        for i in 1..=n_types {
            for j in i..=n_types {
                if let Some(lj) = topo
                    .lj_parameters
                    .get(&(i - 1, j - 1))
                    .or_else(|| topo.lj_parameters.get(&(j - 1, i - 1)))
                {
                    println!(
                        "{:5}{:5}{:15.7e}{:15.7e}{:15.7e}{:15.7e}",
                        i, j, lj.c12, lj.c6, lj.cs12, lj.cs6
                    );
                }
            }
        }
        println!("END");
    }

    // SOLUTEMOLECULES
    println!("SOLUTEMOLECULES");
    println!("# NSPM: number of separate molecules");
    println!("{:5}", topo.solute_molecules.len());
    println!("# NSP: atom sequence number of last atom");
    for &mol_last in &topo.solute_molecules {
        println!("{:5}", mol_last);
    }
    println!("END");

    // TEMPERATUREGROUPS
    println!("TEMPERATUREGROUPS");
    println!("# NSTM: number of temperature atom groups");
    println!("{:5}", topo.temperature_groups.len());
    println!("# NST: atom sequence number of last atom");
    for &tg in &topo.temperature_groups {
        println!("{:5}", tg + 1); // stored 0-based, output 1-based
    }
    println!("END");

    // PRESSUREGROUPS
    println!("PRESSUREGROUPS");
    println!("# NSTM: number of pressure atom groups");
    println!("{:5}", topo.pressure_groups.len());
    println!("# NST: atom sequence number of last atom");
    for &pg in &topo.pressure_groups {
        println!("{:5}", pg);
    }
    println!("END");

    // SOLVENTATOM
    if !topo.solvent_atoms.is_empty() {
        println!("SOLVENTATOM");
        println!("# NRAM: number of atoms per solvent molecule");
        println!("{:5}", topo.solvent_atoms.len());
        println!("#     I  ANMS   IACS  MASS     CGS");
        for (idx, sa) in topo.solvent_atoms.iter().enumerate() {
            println!(
                "{:5} {:>5}{:5}{:10.5}{:12.5}",
                idx + 1,
                sa.name,
                sa.iac,
                sa.mass,
                sa.charge
            );
        }
        println!("END");
    }

    // SOLVENTCONSTR
    if !topo.solvent_constraints.is_empty() {
        println!("SOLVENTCONSTR");
        println!("# NCONS: number of constraints");
        println!("{:5}", topo.solvent_constraints.len());
        println!("#  IC   JC       CC");
        for sc in &topo.solvent_constraints {
            println!("{:5}{:5}{:15.7e}", sc.i + 1, sc.j + 1, sc.length);
        }
        println!("END");
    }
}
