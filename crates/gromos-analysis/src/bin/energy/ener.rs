//! ener — Recalculate interaction energies from a trajectory.
//!
//! Usage:
//!   ener @topo <top> @traj <trc> [@cutoff <nm>] [@epsrf <eps>] [@pbc]

use gromos_forces::energy::{single_point_energy, EnergyParams};
use gromos_core::{
    gather::gather_molecules,
    math::{Periodicity, Rectangular, Vec3},
};
use gromos_io::{
    gromos_args,
    topology::{build_topology, read_topology_file},
    trajectory::TrajectoryReader,
};
use std::process;

fn print_usage() {
    eprintln!("ener — recalculate interaction energies from trajectory");
    eprintln!();
    eprintln!("Usage: ener @topo <top> @traj <trc> [@cutoff <nm>] [@epsrf <eps>] [@pbc]");
    eprintln!();
    eprintln!("  @topo    Topology");
    eprintln!("  @traj    Trajectory (.trc)");
    eprintln!("  @cutoff  Nonbonded cutoff in nm (default: 1.4)");
    eprintln!("  @epsrf   Reaction-field dielectric (default: 78.5)");
    eprintln!("  @pbc     Gather molecules before computing energy");
}

fn main() {
    let args = gromos_args();

    if args.len() < 2 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut params    = EnergyParams::default();
    let mut do_pbc    = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--topo"   => { i += 1; topo_file = Some(args[i].clone()); }
            "--traj"   => { i += 1; traj_file = Some(args[i].clone()); }
            "--cutoff" => { i += 1; params.cutoff = args[i].parse().unwrap_or(1.4); }
            "--epsrf"  => { i += 1; params.rf_epsilon = args[i].parse().unwrap_or(78.5); }
            "--pbc"    => { do_pbc = true; }
            other if other.starts_with("--") => {
                eprintln!("Unknown argument: {other}"); process::exit(1);
            }
            _ => {}
        }
        i += 1;
    }

    let topo_file = topo_file.unwrap_or_else(|| { eprintln!("Error: @topo required"); process::exit(1); });
    let traj_file = traj_file.unwrap_or_else(|| { eprintln!("Error: @traj required"); process::exit(1); });

    let topo_data = read_topology_file(&topo_file).unwrap_or_else(|e| {
        eprintln!("Error reading topology: {e}"); process::exit(1);
    });
    let topo = build_topology(topo_data);
    params.atoms_per_solvent = topo.solvent_atom_template.len().max(3);

    let mut reader = TrajectoryReader::new(&traj_file).unwrap_or_else(|e| {
        eprintln!("Error opening trajectory: {e}"); process::exit(1);
    });

    println!("# ener: topo={topo_file} traj={traj_file}");
    println!("# cutoff={:.3} nm  epsrf={:.1}", params.cutoff, params.rf_epsilon);
    println!("# {:>12}  {:>14}  {:>14}  {:>14}  {:>14}",
        "time(ps)", "E_bond(kJ/mol)", "E_lj(kJ/mol)", "E_crf(kJ/mol)", "E_pot(kJ/mol)");

    while let Ok(Some(frame)) = reader.read_frame() {
        let mut pos = frame.positions.clone();

        if do_pbc && frame.box_dims.x > 0.0 {
            let periodicity = Periodicity::Rectangular(Rectangular::new(frame.box_dims));
            gather_molecules(&mut pos, &topo.molecules, &periodicity, Vec3::ZERO);
        }

        let e = single_point_energy(&topo, &pos, frame.box_dims, &params);

        println!("{:14.6}  {:14.6}  {:14.6}  {:14.6}  {:14.6}",
            frame.time, e.bond, e.lj, e.crf, e.potential);
    }
}
