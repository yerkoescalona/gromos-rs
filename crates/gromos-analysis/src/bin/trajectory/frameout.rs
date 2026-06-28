//! frameout — Write individual snapshots from a trajectory.
//!
//! Full port of GROMOS frameout. All I/O goes through gromos-io.
//!
//! Usage:
//!   frameout @topo <top> @traj <trc>
//!            [@pbc <r|t|v> [gather]]
//!            [@spec <ALL|EVERY <n>|SPEC>] [@frames <list>]
//!            [@include <SOLUTE|SOLVENT|ALL>]
//!            [@ref <cnf>] [@atomsfit <spec>]
//!            [@out <cnf|pdb|trc>] [@name <prefix>] [@single]
//!            [@notimeblock] [@time <t_start> [t_end] [dt]]

use gromos_analysis::fit::superimpose;
use gromos_core::{
    gather::gather_molecules,
    math::{Mat3, Periodicity, Rectangular, Triclinic, Vacuum, Vec3},
    selection::AtomSelection,
};
use gromos_io::{
    coordinate::read_coordinates,
    g96::write_g96,
    gromos_args,
    pdb::write_pdb_positions,
    topology::{build_topology, read_topology_file},
    trajectory::{TrajectoryReader, TrajectoryWriter},
};
use std::process;

fn print_usage() {
    eprintln!("frameout — write individual snapshots from a trajectory");
    eprintln!();
    eprintln!("Usage: frameout @topo <top> @traj <trc> [@pbc <r|t|v>]");
    eprintln!("                [@spec <ALL|EVERY <n>|SPEC>] [@frames <list>]");
    eprintln!("                [@include <SOLUTE|SOLVENT|ALL>]");
    eprintln!("                [@ref <cnf>] [@atomsfit <spec>]");
    eprintln!("                [@out <cnf|pdb|trc>] [@name <prefix>] [@single]");
    eprintln!("                [@notimeblock] [@time <t_start> [t_end] [dt]]");
}

enum FrameSpec {
    All,
    Every(usize),
    Spec(Vec<usize>),
}

impl FrameSpec {
    fn matches(&self, idx: usize) -> bool {
        match self {
            FrameSpec::All => true,
            FrameSpec::Every(n) => idx % n == 0,
            FrameSpec::Spec(list) => list.contains(&idx),
        }
    }
}

fn parse_index_list(s: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let p = part.trim();
        if let Some((a, b)) = p.split_once('-') {
            if let (Ok(a), Ok(b)) = (a.parse::<usize>(), b.parse::<usize>()) {
                out.extend(a..=b);
                continue;
            }
        }
        if let Ok(n) = p.parse::<usize>() {
            out.push(n);
        }
    }
    out
}

fn main() {
    let args = gromos_args();
    if args.len() < 2 || args.contains(&"--help".to_string()) {
        print_usage();
        process::exit(if args.len() < 2 { 1 } else { 0 });
    }

    let mut topo_file = None;
    let mut traj_file = None;
    let mut ref_file = None;
    let mut pbc_type = "v".to_string();
    let mut do_gather = false;
    let mut frame_spec = FrameSpec::All;
    let mut extra_frames: Vec<usize> = Vec::new();
    let mut time_range = (f64::NEG_INFINITY, f64::INFINITY, 0.0_f64);
    let mut include_str = "ALL".to_string();
    let mut fit_spec = None::<String>;
    let mut outfmt = "cnf".to_string();
    let mut name_pfx = "frame".to_string();
    let mut single = false;
    let mut notimeblock = false;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--topo" => {
                i += 1;
                topo_file = Some(args[i].clone());
            },
            "--traj" => {
                i += 1;
                traj_file = Some(args[i].clone());
            },
            "--ref" => {
                i += 1;
                ref_file = Some(args[i].clone());
            },
            "--atomsfit" => {
                i += 1;
                fit_spec = Some(args[i].clone());
            },
            "--out" => {
                i += 1;
                outfmt = args[i].clone();
            },
            "--name" => {
                i += 1;
                name_pfx = args[i].clone();
            },
            "--include" => {
                i += 1;
                include_str = args[i].to_uppercase();
            },
            "--single" => {
                single = true;
            },
            "--notimeblock" => {
                notimeblock = true;
            },
            "--pbc" => {
                i += 1;
                pbc_type = args[i].to_lowercase();
                if i + 1 < args.len() && !args[i + 1].starts_with("--") {
                    i += 1;
                    do_gather = !args[i].eq_ignore_ascii_case("nog");
                } else {
                    do_gather = pbc_type != "v";
                }
            },
            "--spec" => {
                i += 1;
                match args[i].to_uppercase().as_str() {
                    "ALL" => frame_spec = FrameSpec::All,
                    "EVERY" => {
                        i += 1;
                        frame_spec = FrameSpec::Every(args[i].parse().unwrap_or(1));
                    },
                    "SPEC" => { /* will use extra_frames below */ },
                    _ => {},
                }
            },
            "--frames" => {
                i += 1;
                extra_frames = parse_index_list(&args[i]);
            },
            "--time" => {
                let mut tv = Vec::new();
                i += 1;
                while i < args.len() && !args[i].starts_with("--") {
                    if let Ok(v) = args[i].parse::<f64>() {
                        tv.push(v);
                        i += 1;
                    } else {
                        break;
                    }
                }
                time_range = (
                    tv.get(0).copied().unwrap_or(f64::NEG_INFINITY),
                    tv.get(1).copied().unwrap_or(f64::INFINITY),
                    tv.get(2).copied().unwrap_or(0.0),
                );
                continue;
            },
            other if other.starts_with("--") => {
                eprintln!("Unknown argument: {other}");
                process::exit(1);
            },
            _ => {},
        }
        i += 1;
    }

    if !extra_frames.is_empty() {
        if matches!(frame_spec, FrameSpec::All) {
            frame_spec = FrameSpec::Spec(extra_frames);
        }
    }

    let traj_file = traj_file.unwrap_or_else(|| {
        eprintln!("Error: @traj required");
        process::exit(1);
    });

    // Load topology
    let topo = topo_file.as_ref().map(|p| {
        let data = read_topology_file(p).unwrap_or_else(|e| {
            eprintln!("Error reading topology: {e}");
            process::exit(1);
        });
        build_topology(data)
    });
    let topo_ref = topo.as_ref();

    // Reference positions for fit
    let reference = ref_file.as_ref().map(|p| {
        read_coordinates(p)
            .unwrap_or_else(|e| {
                eprintln!("Error reading @ref: {e}");
                process::exit(1);
            })
            .positions
    });

    // Fit atom indices
    let fit_indices: Option<Vec<usize>> = match (&topo_ref, &fit_spec, &reference) {
        (Some(t), Some(spec), _) => Some(
            AtomSelection::from_string(spec, t)
                .unwrap_or_else(|e| {
                    eprintln!("Error in @atomsfit '{spec}': {e}");
                    process::exit(1);
                })
                .indices()
                .to_vec(),
        ),
        (Some(t), None, Some(_)) => Some((0..t.num_solute_atoms()).collect()),
        _ => None,
    };

    // Atom index filter for @include
    let write_indices: Option<Vec<usize>> = topo_ref.map(|t| match include_str.as_str() {
        "SOLUTE" => (0..t.num_solute_atoms()).collect(),
        "SOLVENT" => (t.num_solute_atoms()..t.num_atoms()).collect(),
        _ => (0..t.num_atoms()).collect(),
    });

    // Periodicity builder
    let make_periodicity = |b: Vec3| -> Periodicity {
        match pbc_type.as_str() {
            "r" if b.x > 0.0 => Periodicity::Rectangular(Rectangular::new(b)),
            "t" if b.x > 0.0 => {
                let m = Mat3::from_cols(
                    Vec3::new(b.x, 0., 0.),
                    Vec3::new(0., b.y, 0.),
                    Vec3::new(0., 0., b.z),
                );
                Periodicity::Triclinic(Triclinic::new(m))
            },
            _ => Periodicity::Vacuum(Vacuum),
        }
    };

    // Open input
    let mut reader = TrajectoryReader::new(&traj_file).unwrap_or_else(|e| {
        eprintln!("Error opening trajectory: {e}");
        process::exit(1);
    });

    // Single-file .trc writer
    let mut trc_writer: Option<TrajectoryWriter> = if single && outfmt == "trc" {
        let p = format!("{name_pfx}.trc");
        Some(
            TrajectoryWriter::new(&p, "frameout sub-trajectory", false, false).unwrap_or_else(
                |e| {
                    eprintln!("Error creating {p}: {e}");
                    process::exit(1);
                },
            ),
        )
    } else {
        None
    };

    let (t_start, t_end, t_dt) = time_range;
    let mut frame_idx = 0usize;
    let mut written = 0usize;

    loop {
        match reader.read_frame() {
            Ok(None) => break,
            Err(e) => {
                eprintln!("Read error: {e}");
                break;
            },
            Ok(Some(frame)) => {
                let t = frame.time;
                let time_ok = t >= t_start - 1e-9
                    && t <= t_end + 1e-9
                    && (t_dt <= 0.0 || ((t - t_start) / t_dt).round() * t_dt + t_start - t < 1e-6);

                if frame_spec.matches(frame_idx) && time_ok {
                    // 1. Gather
                    let mut pos = frame.positions.clone();
                    if do_gather && frame.box_dims.x > 0.0 {
                        if let Some(ref t) = topo {
                            let per = make_periodicity(frame.box_dims);
                            gather_molecules(&mut pos, &t.molecules, &per, Vec3::ZERO);
                        }
                    }

                    // 2. Rotational fit
                    if let (Some(ref refpos), Some(ref fidx)) = (&reference, &fit_indices) {
                        superimpose(&mut pos, refpos, fidx, None);
                    }

                    // 3. Filter atoms
                    let out_pos: Vec<Vec3> = if let Some(ref idx) = write_indices {
                        idx.iter().map(|&i| pos[i]).collect()
                    } else {
                        pos
                    };

                    let box_opt = if frame.box_dims.x > 0.0 {
                        Some(frame.box_dims)
                    } else {
                        None
                    };

                    // Topology for labelling (only when writing all atoms)
                    let topo_for_write = if write_indices.as_ref().map_or(false, |idx| {
                        idx.len() < topo_ref.map_or(0, |t| t.num_atoms())
                    }) {
                        None
                    } else {
                        topo_ref
                    };

                    // 4. Write via gromos-io
                    let tag = format!("{frame_idx:06}_{:.6}", frame.time).replace('.', "_");

                    if single {
                        if let Some(ref mut tw) = trc_writer {
                            tw.write_trc_frame(frame.step, frame.time, &out_pos, box_opt)
                                .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                        } else {
                            // cnf single: last frame wins
                            let title = if notimeblock {
                                "frameout".to_string()
                            } else {
                                format!("t={:.9} ps step={}", frame.time, frame.step)
                            };
                            write_g96(
                                format!("{name_pfx}.cnf"),
                                &title,
                                &out_pos,
                                None,
                                box_opt,
                                topo_for_write,
                            )
                            .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                        }
                    } else {
                        let ext = match outfmt.as_str() {
                            "pdb" => "pdb",
                            "trc" => "trc",
                            _ => "cnf",
                        };
                        let path = format!("{name_pfx}_{tag}.{ext}");
                        match outfmt.as_str() {
                            "pdb" => {
                                let title = format!("t={:.6} ps step={}", frame.time, frame.step);
                                write_pdb_positions(
                                    &path,
                                    &title,
                                    &out_pos,
                                    box_opt,
                                    topo_for_write,
                                )
                                .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                            },
                            "trc" => {
                                let mut tw = TrajectoryWriter::new(&path, "frameout", false, false)
                                    .unwrap_or_else(|e| {
                                        eprintln!("Error: {e}");
                                        process::exit(1);
                                    });
                                tw.write_trc_frame(frame.step, frame.time, &out_pos, box_opt)
                                    .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                                tw.flush().ok();
                            },
                            _ => {
                                let title = if notimeblock {
                                    "frameout".to_string()
                                } else {
                                    format!("t={:.9} ps step={}", frame.time, frame.step)
                                };
                                write_g96(&path, &title, &out_pos, None, box_opt, topo_for_write)
                                    .unwrap_or_else(|e| eprintln!("Write error: {e}"));
                            },
                        }
                        eprintln!("# frame {frame_idx} t={:.4} ps → {path}", frame.time);
                    }
                    written += 1;
                }
                frame_idx += 1;
            },
        }
    }

    if let Some(ref mut tw) = trc_writer {
        tw.flush().ok();
    }
    eprintln!("# frameout: {written} frames written ({frame_idx} total read)");
}
