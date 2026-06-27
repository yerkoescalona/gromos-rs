//! nhoparam — N-H bond order parameters for NMR analysis.
//!
//! Calculates the Lipari-Szabo order parameter S² for N-H bonds:
//!
//!   S² = ½ [ 3·Σᵢⱼ ⟨μᵢ(t)μⱼ(t)⟩² − 1 ]
//!
//! where μ is the unit vector along the N-H bond after rotational fitting.
//!
//! Reference: Lipari & Szabo (1982) J. Am. Chem. Soc. 104:4546.
//! Algorithm: faithful port of gromos-rs nhoparam (GROMOS analysis tools).
//!
//! Usage:
//!   nhoparam @topo <top> @traj <trc> @atoms <spec>
//!            [@ref <cnf>] [@atomsfit <spec>] [@winframe <n>] [@pbc]

use gromos_analysis::fit::superimpose;
use gromos_core::{
    gather::gather_molecules,
    math::{Periodicity, Rectangular, Vacuum, Vec3},
    selection::AtomSelection,
    stat::Stat,
};
use gromos_io::{
    coordinate::read_coordinates,
    gromos_args,
    topology::{build_topology, read_topology_file},
    trajectory::TrajectoryReader,
};
use std::process;

const H_MASS: f64 = 1.10; // upper bound for hydrogen mass (Da)

fn print_usage() {
    eprintln!("nhoparam — N-H bond order parameters (S²)");
    eprintln!();
    eprintln!("Usage: nhoparam @topo <top> @traj <trc> @atoms <spec>");
    eprintln!("                [@ref <cnf>] [@atomsfit <spec>] [@winframe <n>] [@pbc]");
    eprintln!();
    eprintln!("  @topo      Topology");
    eprintln!("  @traj      Trajectory (.trc)");
    eprintln!("  @atoms     Nitrogen atoms to analyse (e.g. '1:N')");
    eprintln!("  @ref       Reference coordinates for fit (default: first frame)");
    eprintln!("  @atomsfit  Atoms used for rotational fit (default: @atoms)");
    eprintln!("  @winframe  Window size for running average (default: 10)");
    eprintln!("  @pbc       Gather molecules before computing vectors");
}

/// S² from the accumulated tensor sum and number of frames.
fn s2_from_tensor(tensor: &[[f64; 3]; 3], n: usize) -> f64 {
    if n == 0 { return 0.0; }
    let mut sum = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            let ave = tensor[a][b] / n as f64;
            sum += ave * ave;
        }
    }
    0.5 * (3.0 * sum - 1.0)
}

fn main() {
    let args = gromos_args();

    if args.len() < 2 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut topo_file  = None;
    let mut traj_file  = None;
    let mut ref_file   = None;
    let mut atom_spec  = None;
    let mut fit_spec   = None;
    let mut winframe   = 10usize;
    let mut do_pbc     = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--topo"     => { i += 1; topo_file  = Some(args[i].clone()); }
            "--traj"     => { i += 1; traj_file  = Some(args[i].clone()); }
            "--ref"      => { i += 1; ref_file   = Some(args[i].clone()); }
            "--atoms"    => { i += 1; atom_spec  = Some(args[i].clone()); }
            "--atomsfit" => { i += 1; fit_spec   = Some(args[i].clone()); }
            "--winframe" => { i += 1; winframe   = args[i].parse().unwrap_or(10); }
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
    let atom_spec = atom_spec.unwrap_or_else(|| { eprintln!("Error: @atoms required"); process::exit(1); });

    // Load topology
    let topo_data = read_topology_file(&topo_file).unwrap_or_else(|e| {
        eprintln!("Error reading topology: {e}"); process::exit(1);
    });
    let topo = build_topology(topo_data);

    // Resolve N atoms selection
    let n_sel = AtomSelection::from_string(&atom_spec, &topo).unwrap_or_else(|e| {
        eprintln!("Error in @atoms '{atom_spec}': {e}"); process::exit(1);
    });
    let n_indices: Vec<usize> = n_sel.indices().to_vec();
    eprintln!("# N atoms selected: {}", n_indices.len());

    // Verify all selected atoms are nitrogen (mass > H_MASS)
    for &ni in &n_indices {
        let mass = topo.mass.get(ni).copied().unwrap_or(0.0);
        if mass < 12.0 {
            eprintln!("Warning: atom {ni} has mass {mass:.3} — may not be nitrogen");
        }
    }

    // For each N, find bonded H atoms (mass < H_MASS)
    let mut nh_pairs: Vec<(usize, Vec<usize>)> = Vec::new(); // (N_idx, [H_idx, ...])
    for &ni in &n_indices {
        let mut h_atoms: Vec<usize> = Vec::new();
        for bond in &topo.moltypes[0].bonds {
            let (a, b) = (bond.i, bond.j);
            let partner = if a == ni { b } else if b == ni { a } else { continue };
            let h_mass = topo.mass.get(partner).copied().unwrap_or(99.0);
            if h_mass < H_MASS {
                h_atoms.push(partner);
            }
        }
        if h_atoms.is_empty() {
            eprintln!("Warning: N atom {ni} has no bonded H — skipping");
        } else {
            nh_pairs.push((ni, h_atoms));
        }
    }
    eprintln!("# N-H pairs found: {}", nh_pairs.len());

    // Resolve fit atoms
    let fit_spec_str = fit_spec.as_deref().unwrap_or(&atom_spec);
    let fit_sel = AtomSelection::from_string(fit_spec_str, &topo).unwrap_or_else(|e| {
        eprintln!("Error in @atomsfit '{fit_spec_str}': {e}"); process::exit(1);
    });
    let fit_indices: Vec<usize> = fit_sel.indices().to_vec();

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
    eprintln!("# Frames: {}", frames.len());

    // Reference positions: from @ref file or first frame
    let reference: Vec<Vec3> = if let Some(ref path) = ref_file {
        read_coordinates(path).unwrap_or_else(|e| {
            eprintln!("Error reading @ref: {e}"); process::exit(1);
        }).positions
    } else {
        frames[0].positions.clone()
    };

    // Accumulate order parameter tensors (one 3×3 per N-H pair)
    let np = nh_pairs.len();
    let mut tensors: Vec<[[f64; 3]; 3]> = vec![[[0.0; 3]; 3]; np];
    let mut win_tensors: Vec<[[f64; 3]; 3]> = vec![[[0.0; 3]; 3]; np];
    let mut win_stats: Vec<Stat> = (0..np).map(|_| Stat::new()).collect();
    let mut n_frames = 0usize;
    let mut win_count = 0usize;

    for frame in &frames {
        let mut pos = frame.positions.clone();

        // PBC gathering
        if do_pbc {
            let periodicity = if frame.box_dims.x > 0.0 {
                Periodicity::Rectangular(Rectangular::new(frame.box_dims))
            } else {
                Periodicity::Vacuum(Vacuum)
            };
            gather_molecules(&mut pos, &topo.molecules, &periodicity, Vec3::ZERO);
        }

        // Rotational fit
        superimpose(&mut pos, &reference, &fit_indices, None);

        // Accumulate N-H tensors
        for (k, (ni, h_list)) in nh_pairs.iter().enumerate() {
            // Average unit vector over bonded H atoms
            let mut mu = Vec3::ZERO;
            for &hi in h_list {
                let v = pos[hi] - pos[*ni];
                let len = v.length();
                if len > 1e-10 { mu += v / len; }
            }
            // Normalize average (handles both single and multiple H)
            let mu_len = mu.length();
            if mu_len < 1e-10 { continue; }
            let mu = mu / mu_len;

            let mv = [mu.x, mu.y, mu.z];
            for a in 0..3 {
                for b in 0..3 {
                    tensors[k][a][b]     += mv[a] * mv[b];
                    win_tensors[k][a][b] += mv[a] * mv[b];
                }
            }
        }

        n_frames += 1;
        win_count += 1;

        // Window averaging
        if win_count >= winframe {
            for (k, ws) in win_stats.iter_mut().enumerate() {
                ws.add(s2_from_tensor(&win_tensors[k], win_count));
                win_tensors[k] = [[0.0; 3]; 3];
            }
            win_count = 0;
        }
    }

    // Flush incomplete window
    if win_count > 0 {
        for (k, ws) in win_stats.iter_mut().enumerate() {
            ws.add(s2_from_tensor(&win_tensors[k], win_count));
        }
    }

    // Output
    println!("# nhoparam — N-H order parameters (S²)");
    println!("# topo={topo_file}  traj={traj_file}  frames={n_frames}  winframe={winframe}");
    println!("# {:>10} {:>6} {:>8} {:>10} {:>10} {:>10} {:>10}",
        "ATOM", "RESNUM", "RESNAME", "S2", "WINAV", "STDDEV", "EE");

    for (k, (ni, _h_list)) in nh_pairs.iter().enumerate() {
        let s2 = s2_from_tensor(&tensors[k], n_frames);
        let res_nr   = topo.residue_nr(*ni).unwrap_or(0);
        let res_name = topo.residue_name(*ni).unwrap_or("UNK");
        let atom_name = topo.atom_name(*ni).unwrap_or("N");
        let label = format!("{}:{}", 1, ni + 1); // mol:atom (1-based)

        let (winav, stddev, ee) = if win_stats[k].n() >= 2 {
            (win_stats[k].ave(), win_stats[k].rmsd(), win_stats[k].ee())
        } else if win_stats[k].n() == 1 {
            (win_stats[k].ave(), 0.0, 0.0)
        } else {
            (s2, 0.0, 0.0)
        };

        println!("  {:>10} {:>6} {:>8} {:>10.6} {:>10.6} {:>10.6} {:>10.6}",
            label, res_nr, format!("{res_name}/{atom_name}"), s2, winav, stddev, ee);
    }
}
