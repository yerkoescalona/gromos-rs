//! copy_box — duplicate a periodic box along its axes (gromos++ `copy_box`).
//!
//! Usage: copy_box @topo <file> @pos <file> @dir <x[:N],y[:N],z[:N]>
//!    or: copy_box @topo <file> @traj <file> @replicate <nx> <ny> <nz>
//!
//! As gromos++ does: the images are built direction by direction (each existing image is
//! shifted by 1..N box lengths), the solute molecules of every image come first, then the
//! solvent of every image, and the box grows by N lengths along each duplicated axis.

use gromos::math::Vec3;
use gromos_io::topology::read_topology_file;
use gromos_io::{read_g96_labeled, G96Atom};
use std::env;
use std::process;

fn usage() -> ! {
    eprintln!("Usage: copy_box @topo <file> @pos <file> @dir <x[:N],y[:N],z[:N]>");
    eprintln!("   or: copy_box @topo <file> @traj <file> @replicate <nx> <ny> <nz>");
    process::exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut topo_file = None;
    let mut pos_file = None;
    let mut dirs: Vec<(char, usize)> = Vec::new();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => {
                i += 1;
                topo_file = args.get(i).cloned();
            },
            "@pos" | "@traj" => {
                i += 1;
                pos_file = args.get(i).cloned();
            },
            "@dir" => {
                i += 1;
                for item in args.get(i).cloned().unwrap_or_default().split(',') {
                    let (axis, n) = item.split_once(':').unwrap_or((item, "1"));
                    let axis = axis.trim().chars().next().unwrap_or_else(|| usage());
                    if !"xyz".contains(axis) {
                        eprintln!(
                            "Error: only x, y and z directions are supported (rectangular boxes)"
                        );
                        process::exit(1);
                    }
                    let n: usize = n.trim().parse().unwrap_or_else(|_| usage());
                    if n == 0 {
                        eprintln!("Error: multiplicity in @dir must be >= 1");
                        process::exit(1);
                    }
                    dirs.push((axis, n));
                }
            },
            "@replicate" => {
                for (k, axis) in ['x', 'y', 'z'].into_iter().enumerate() {
                    let n: usize = args
                        .get(i + 1 + k)
                        .and_then(|s| s.parse().ok())
                        .unwrap_or_else(|| usage());
                    if n > 1 {
                        dirs.push((axis, n - 1));
                    }
                }
                i += 3;
            },
            _ => {},
        }
        i += 1;
    }
    let (Some(topo_file), Some(pos_file)) = (topo_file, pos_file) else {
        usage()
    };
    let topo = match read_topology_file(&topo_file) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Error reading topology '{}': {}", topo_file, e);
            process::exit(1);
        },
    };
    let data = match read_g96_labeled(&pos_file) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("Error reading '{}': {}", pos_file, e);
            process::exit(1);
        },
    };
    if data.box_dims == Vec3::ZERO {
        eprintln!("Error: '{}' has no box", pos_file);
        process::exit(1);
    }
    let n_solute = topo.n_atoms;
    if data.atoms.len() < n_solute {
        eprintln!(
            "Error: {} atoms in the coordinates, {} in the topology",
            data.atoms.len(),
            n_solute
        );
        process::exit(1);
    }
    let per_solvent = topo.solvent_atoms.len().max(1);
    let base = data.box_dims;
    let axis_shift = |axis: char| match axis {
        'x' => Vec3::new(base.x, 0.0, 0.0),
        'y' => Vec3::new(0.0, base.y, 0.0),
        _ => Vec3::new(0.0, 0.0, base.z),
    };
    // images, as gromos++ builds them: for each direction, every existing image shifted 1..N times
    let mut images: Vec<Vec3> = vec![Vec3::ZERO];
    for &(axis, n) in &dirs {
        let mut new_images = Vec::new();
        for k in 1..=n {
            let shift = axis_shift(axis) * k as f64;
            for &im in &images {
                new_images.push(im + shift);
            }
        }
        images.extend(new_images);
    }
    let mut new_box = base;
    for &(axis, n) in &dirs {
        match axis {
            'x' => new_box.x += base.x * n as f64,
            'y' => new_box.y += base.y * n as f64,
            _ => new_box.z += base.z * n as f64,
        }
    }
    let n_res = topo.residue_names.len();
    let (solute, solvent) = data.atoms.split_at(n_solute);
    let n_solvent_mols = solvent.len() / per_solvent;
    eprintln!(
        "# {} images of {} solute atoms and {} solvent molecules; box ({:.6}, {:.6}, {:.6})",
        images.len(),
        n_solute,
        n_solvent_mols,
        new_box.x,
        new_box.y,
        new_box.z
    );
    println!("TITLE");
    println!(
        "Copy_box: {} duplicated with dir={} ({} copies)",
        pos_file,
        dirs.iter()
            .map(|(a, n)| format!("{a}:{n}"))
            .collect::<Vec<_>>()
            .join(","),
        images.len()
    );
    println!("END");
    println!("POSITION");
    let line = |res: usize, res_name: &str, atom: &G96Atom, serial: usize, p: Vec3| {
        println!(
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}",
            res, res_name, atom.atom_name, serial, p.x, p.y, p.z
        );
    };
    let mut serial = 0;
    for (k, &shift) in images.iter().enumerate() {
        for a in solute {
            serial += 1;
            line(a.res_num + k * n_res, &a.res_name, a, serial, a.pos + shift);
        }
    }
    let mut mol = 0;
    for &shift in &images {
        for (j, a) in solvent.iter().enumerate() {
            if j % per_solvent == 0 {
                mol += 1;
            }
            serial += 1;
            line(mol, "SOLV", a, serial, a.pos + shift);
        }
    }
    println!("END");
    println!("GENBOX");
    println!("{:>8}", 1);
    println!("{:15.9}{:15.9}{:15.9}", new_box.x, new_box.y, new_box.z);
    println!("{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0);
    println!("{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0);
    println!("{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0);
    println!("END");
}
