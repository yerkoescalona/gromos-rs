//! rmsd — Root Mean Square Deviation (gromos-rs-compatible).
//!
//! Usage:
//!   rmsd @topo <topology> @traj <trajectory>
//!        [@ref <frame_spec>] [@atomspec <spec>] [@nofit] [@pbc]
//!
//!   @topo      Topology file (.top/.topo)
//!   @traj      Trajectory file (.trc)
//!   @ref       Reference frame: "first" (default), "last", or 0-based integer
//!   @atomspec  gromos-rs AtomSpecifier for fit+RMSD atoms (default: all)
//!              e.g. "1:a", "1:CA,N,C", "1:res(1-10:a)"
//!   @nofit     Skip rotational fit (translation-only centering still applied)
//!   @pbc       Gather molecules before computing RMSD (default: off)
//!
//! Output:
//!   # time(ps)  RMSD(nm)
//!     0.000     0.000000
//!     0.002     0.012345

use gromos_analysis::fit::{rmsd, superimpose};
use gromos_core::{
    gather::gather_molecules,
    math::Periodicity,
    selection::AtomSelection,
};
use gromos_io::{
    gromos_args,
    topology::{build_topology, read_topology_file},
    trajectory::TrajectoryReader,
};
use std::process;

fn print_usage() {
    eprintln!("rmsd — Root Mean Square Deviation");
    eprintln!();
    eprintln!("Usage: rmsd @topo <top> @traj <trc> [@ref <frame>] [@atomspec <spec>] [@nofit] [@pbc]");
    eprintln!();
    eprintln!("  @topo      Topology (.top/.topo)");
    eprintln!("  @traj      Trajectory (.trc)");
    eprintln!("  @ref       Reference: 'first' (default), 'last', or 0-based integer");
    eprintln!("  @atomspec  AtomSpecifier for fit + RMSD (default: 1:a = all solute)");
    eprintln!("  @nofit     Skip rotational fit");
    eprintln!("  @pbc       Gather molecules before RMSD");
}

fn main() {
    let args = gromos_args();

    if args.len() < 2 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut ref_spec  = "first".to_string();
    let mut atomspec  = "1:a".to_string();
    let mut do_fit    = true;
    let mut do_pbc    = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--topo"     => { i += 1; topo_file = Some(args[i].clone()); }
            "--traj"     => { i += 1; traj_file = Some(args[i].clone()); }
            "--ref"      => { i += 1; ref_spec  = args[i].clone(); }
            "--atomspec" => { i += 1; atomspec  = args[i].clone(); }
            "--nofit"    => { do_fit = false; }
            "--pbc"      => { do_pbc = true; }
            other if other.starts_with("--") => {
                eprintln!("Unknown argument: {other}"); process::exit(1);
            }
            _ => {}
        }
        i += 1;
    }

    let topo_file = topo_file.unwrap_or_else(|| { eprintln!("Error: @topo required"); process::exit(1); });
    let traj_file = traj_file.unwrap_or_else(|| { eprintln!("Error: @traj required"); process::exit(1); });

    // Load topology
    let topo_data = read_topology_file(&topo_file).unwrap_or_else(|e| {
        eprintln!("Error reading topology: {e}"); process::exit(1);
    });
    let topo = build_topology(topo_data);

    // Resolve atom selection
    let sel = AtomSelection::from_string(&atomspec, &topo).unwrap_or_else(|e| {
        eprintln!("Error in @atomspec '{atomspec}': {e}"); process::exit(1);
    });
    let fit_indices: Vec<usize> = sel.indices().to_vec();
    eprintln!("# topology: {} atoms, {} selected by '{atomspec}'",
        topo.num_atoms(), fit_indices.len());

    // Read all trajectory frames
    let mut reader = TrajectoryReader::new(&traj_file).unwrap_or_else(|e| {
        eprintln!("Error opening trajectory: {e}"); process::exit(1);
    });
    let frames = reader.read_all_frames().unwrap_or_else(|e| {
        eprintln!("Error reading trajectory: {e}"); process::exit(1);
    });
    if frames.is_empty() {
        eprintln!("Error: no frames in trajectory"); process::exit(1);
    }
    eprintln!("# frames: {}", frames.len());

    // Choose reference frame
    let ref_idx = match ref_spec.to_lowercase().as_str() {
        "first" => 0,
        "last"  => frames.len() - 1,
        s => s.parse::<usize>().unwrap_or_else(|_| {
            eprintln!("Error: invalid @ref '{s}' — use 'first', 'last', or integer"); process::exit(1);
        }),
    };
    if ref_idx >= frames.len() {
        eprintln!("Error: @ref {ref_idx} out of range (0..{})", frames.len()-1);
        process::exit(1);
    }
    eprintln!("# reference frame: {ref_idx} (t={:.4} ps)", frames[ref_idx].time);

    // Reference positions (gather if PBC, then use as-is — reference is fixed)
    let mut reference = frames[ref_idx].positions.clone();
    if do_pbc {
        let periodicity = {
            let b = frames[ref_idx].box_dims;
            if b.x > 0.0 {
                Periodicity::Rectangular(gromos_core::math::Rectangular::new(b))
            } else {
                Periodicity::Vacuum(gromos_core::math::Vacuum)
            }
        };
        gather_molecules(&mut reference, &topo.molecules, &periodicity, gromos_core::math::Vec3::ZERO);
    }

    // Output header
    println!("# rmsd: topo={topo_file} traj={traj_file} ref={ref_idx} atoms={} fit={do_fit}",
        fit_indices.len());
    println!("# {:>12}  {:>14}", "time(ps)", "RMSD(nm)");

    // Compute RMSD per frame
    for frame in &frames {
        let mut pos = frame.positions.clone();

        if do_pbc {
            let periodicity = {
                let b = frame.box_dims;
                if b.x > 0.0 {
                    Periodicity::Rectangular(gromos_core::math::Rectangular::new(b))
                } else {
                    Periodicity::Vacuum(gromos_core::math::Vacuum)
                }
            };
            gather_molecules(&mut pos, &topo.molecules, &periodicity, gromos_core::math::Vec3::ZERO);
        }

        let d = if do_fit {
            superimpose(&mut pos, &reference, &fit_indices, None);
            rmsd(&pos, &reference, &fit_indices)
        } else {
            rmsd(&pos, &reference, &fit_indices)
        };

        println!("{:14.6}  {:14.9}", frame.time, d);
    }
}
