//! duplicate — find atoms that sit exactly on top of each other (gromos++ `duplicate`).
//!
//! Usage: duplicate @topo <top> @pos <cnf> @pbc <r|v> [@write]
//!
//! Reports coincident solute atoms (across molecules) and solvent atoms (across molecules),
//! molecule by molecule as gromos++ does; with `@write` the system without the duplicated
//! molecules is written instead of the report.

use gromos_io::args::{fail, Arguments};
use gromos_io::coordinate::format_g96;
use gromos_io::pbc::Pbc;
use gromos_io::read_g96_labeled;
use gromos_io::topology::read_topology_file;

const USAGE: &str = "# duplicate
\t@topo       <input topology file>
\t@pos         <input coordinate file>
\t@pbc         <periodic boudary conditions>
\t[@write      <write out duplicate-filtered coordinate file, flag>]";

fn main() {
    if let Err(e) = run() {
        fail(&e);
    }
}

fn run() -> Result<(), String> {
    let args = Arguments::parse(&["topo", "pos", "pbc", "write"], USAGE)?;
    let topo = read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?;
    let data = read_g96_labeled(args.value("pos")?).map_err(|e| format!("@pos: {e}"))?;
    let pbc = Pbc::from_args(&args)?;
    let periodicity = pbc.periodicity(data.box_dims);
    let write = args.count("write") >= 0;
    let atoms = &data.atoms;
    let n_solute = topo.n_atoms.min(atoms.len());
    let mols: Vec<std::ops::Range<usize>> = {
        let mut v = Vec::new();
        let mut start = 0;
        for &end in &topo.solute_molecules {
            v.push(start..end.min(n_solute));
            start = end;
        }
        if v.is_empty() && n_solute > 0 {
            v.push(0..n_solute);
        }
        v
    };
    let per = topo.solvent_atoms.len().max(1);
    let n_solvent_mols = (atoms.len() - n_solute) / per;
    let mut report = String::new();
    let header = "Duplicated atoms (mol1:atom1 <-> mol2:atom2)\n";
    if write {
        report.push_str(header);
        report.push_str("-------------------------------------------\n");
    } else {
        report.push('\n');
        report.push_str(header);
        report.push_str("--------------------------------------------\n");
    }
    // gromos++: pos_a − nearestImage(pos_a, pos_b) == 0 ⇔ the minimum-image vector is zero
    let same = |a: usize, b: usize| {
        periodicity.nearest_image(atoms[a].pos, atoms[b].pos) == gromos_core::Vec3::ZERO
    };
    let mut delmol = Vec::new();
    let mut keep_solute = Vec::new();
    for (i, mi) in mols.iter().enumerate() {
        let mut addmol = true;
        for (j, a) in mi.clone().enumerate() {
            for (k, mk) in mols.iter().enumerate().skip(i) {
                // gromos++ starts the inner atom loop at j even in another molecule
                for (l, b) in mk.clone().enumerate().skip(j) {
                    if same(a, b) && (j != l || i != k) {
                        addmol = false;
                        report.push_str(&format!(
                            "solute: {}:{}  <->  {}:{}\n",
                            i + 1,
                            j + 1,
                            k + 1,
                            l + 1
                        ));
                    }
                }
            }
        }
        if addmol {
            keep_solute.push(i);
        } else {
            delmol.push(i + 1);
        }
    }
    let mut numdelsol = 0;
    let mut keep_solvent = Vec::new();
    for j in 0..n_solvent_mols {
        let mut addsol = true;
        for k in 0..per {
            let a = n_solute + j * per + k;
            for m in j..n_solvent_mols {
                for n in 0..per {
                    let b = n_solute + m * per + n;
                    if a != b && same(a, b) {
                        addsol = false;
                        report.push_str(&format!(
                            "solvent: {}:{}  <->  {}:{}\n",
                            j,
                            j * per + k,
                            m,
                            m * per + n
                        ));
                    }
                }
            }
        }
        if addsol {
            keep_solvent.push(j);
        } else {
            numdelsol += 1;
        }
    }
    if write {
        report.push_str(&format!(
            "\nDeleted Solute-Molecules:  {:>6}\n",
            delmol.len()
        ));
        for m in &delmol {
            report.push_str(&format!("  - Molecule {m}\n"));
        }
        report.push_str(&format!("Deleted Solvent-Molecules: {:>6}", numdelsol));
        let mut out = Vec::new();
        let mut offres = 0usize;
        for &i in &keep_solute {
            let first_res = atoms[mols[i].start].res_num;
            let mut n_res = 0;
            for a in mols[i].clone() {
                let mut atom = atoms[a].clone();
                atom.res_num = atom.res_num - first_res + 1 + offres;
                n_res = n_res.max(atom.res_num - offres);
                out.push(atom);
            }
            offres += n_res;
        }
        for (r, &j) in keep_solvent.iter().enumerate() {
            for k in 0..per {
                let mut a = atoms[n_solute + j * per + k].clone();
                a.res_num = r + 1;
                out.push(a);
            }
        }
        let box_dims = (data.box_dims != gromos_core::Vec3::ZERO).then_some(data.box_dims);
        print!("{}", format_g96(&report, &out, box_dims));
    } else {
        report.push_str(&format!(
            "\nDuplicated Solute-Molecules:  {:>6}\n",
            delmol.len()
        ));
        for m in &delmol {
            report.push_str(&format!("  - Molecule {m}\n"));
        }
        report.push_str(&format!("Duplicated Solvent-Molecules: {:>6}\n", numdelsol));
        print!("{report}");
    }
    Ok(())
}
