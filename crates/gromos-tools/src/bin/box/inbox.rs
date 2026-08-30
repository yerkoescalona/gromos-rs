//! inbox — put every atom into the periodic box (gromos++ `inbox`).
//!
//! Usage: inbox @topo <file> @pos <file> @pbc <r|v> [@shift <x> <y> <z>]
//!    or: inbox @topo <file> @traj <file> @pbc <r|v>
//!
//! gromos++: each atom is replaced by its periodic image nearest to the box centre
//! (K/2, L/2, M/2), after an optional shift; vacuum leaves the coordinates alone.

use gromos::math::Vec3;
use gromos_io::read_g96_labeled;
use std::env;
use std::process;

fn usage() -> ! {
    eprintln!("Usage: inbox @topo <file> @pos <file> @pbc <r|v> [@shift <x> <y> <z>]");
    process::exit(1);
}

/// gromos++ rectangular `nearestImage(r, centre)`: the image of `r` nearest to `centre`.
fn nearest_image_to(r: Vec3, centre: Vec3, b: Vec3) -> Vec3 {
    let d = r - centre;
    Vec3::new(
        r.x - b.x * (d.x / b.x).round(),
        r.y - b.y * (d.y / b.y).round(),
        r.z - b.z * (d.z / b.z).round(),
    )
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut pos_file = None;
    let mut pbc = 'r';
    let mut shift = Vec3::ZERO;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "@topo" => i += 1, // the topology is not needed: every atom is treated alike
            "@pos" | "@traj" => {
                i += 1;
                pos_file = args.get(i).cloned();
            },
            "@pbc" => {
                i += 1;
                pbc = args
                    .get(i)
                    .and_then(|s| s.chars().next())
                    .unwrap_or_else(|| usage());
            },
            "@shift" => {
                let v: Vec<f64> = (1..=3)
                    .filter_map(|k| args.get(i + k).and_then(|s| s.parse().ok()))
                    .collect();
                if v.len() != 3 {
                    usage();
                }
                shift = Vec3::new(v[0], v[1], v[2]);
                i += 3;
            },
            _ => {},
        }
        i += 1;
    }
    let Some(pos_file) = pos_file else { usage() };
    let data = match read_g96_labeled(&pos_file) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("Error reading '{}': {}", pos_file, e);
            process::exit(1);
        },
    };
    let b = data.box_dims;
    if pbc != 'v' && b == Vec3::ZERO {
        eprintln!("Error: no box in '{}'", pos_file);
        process::exit(1);
    }
    let centre = b * 0.5;
    println!("TITLE");
    println!("inbox: {}", pos_file);
    println!("END");
    println!("POSITION");
    for (k, a) in data.atoms.iter().enumerate() {
        let p = a.pos + shift;
        let p = if pbc == 'v' {
            p
        } else {
            nearest_image_to(p, centre, b)
        };
        println!(
            "{:>5} {:<5} {:<6}{:>6}{:15.9}{:15.9}{:15.9}",
            a.res_num,
            a.res_name,
            a.atom_name,
            k + 1,
            p.x,
            p.y,
            p.z
        );
    }
    println!("END");
    println!("GENBOX");
    println!("{:>8}", if pbc == 'v' { 0 } else { 1 });
    println!("{:15.9}{:15.9}{:15.9}", b.x, b.y, b.z);
    println!("{:15.9}{:15.9}{:15.9}", 90.0, 90.0, 90.0);
    println!("{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0);
    println!("{:15.9}{:15.9}{:15.9}", 0.0, 0.0, 0.0);
    println!("END");
}
