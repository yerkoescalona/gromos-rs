//! ene_ana — averages and time series of properties from energy trajectories
//! (gromos++ `ene_ana`, `programs/ene_ana.cc`).
//!
//! The layout of a `.tre`/`.trg` is not built in: `@library` declares it and names the
//! properties, which is what lets one library file span GROMOS versions. Without `@library`
//! the one `gromos-io` ships is used. Before a file is read its layout is checked against the
//! library — from the file's own self-description, or from the built-in layout for its
//! `ENEVERSION` — and a mismatch is refused with a diff (`gromos_io::energy_traj`, and
//! `docs/src/reference/energy-library.md`).
//!
//! Each `@prop` is written as a `<prop>.dat` time series and summarised on stdout as
//! average / rmsd / error estimate. The error estimate is gromos++'s block-averaging one,
//! which needs more than 200 frames — with fewer, gromos++ prints `-nan` and so do we.

use gromos_analysis::args::{fail, Arguments};
use gromos_core::stat::Stat;
use gromos_io::energy_traj::{
    builtin_schema, library_text, official_variables, read_preamble, EnergyTraj, OFFICIAL_LIBRARY,
};
use gromos_io::gz::TextReader;
use gromos_io::table::cpp_g;
use gromos_io::topology::{build_topology, read_topology_file};
use std::fmt::Write as _;

const USAGE: &str = "# ene_ana
\t@en_files    <energy files> (and/or)
\t@fr_files    <free energy files>
\t@prop        <properties to monitor>
\t[@library    <library for property names> [print]] (default: the library gromos-io ships)
\t[@topo       <molecular topology file> (for MASS and NUMMOL)]
\t[@time       <t and dt> (overwrites TIME in the trajectory files)]
\t[@print_library [<trajectory>]] (write the library to stdout and stop: the built-in one,
\t             or one generated from the layout the trajectory describes itself with)";

fn main() {
    if let Err(e) = run() {
        fail(&e);
    }
}

/// gromos++ writes these with `precision(9)` and no `fixed`, i.e. C++ default `%g`.
fn num(x: f64) -> String {
    if x.is_nan() {
        return format!("{:>15}", "-nan");
    }
    format!("{:>15}", cpp_g(x, 9))
}

/// The `@en_files`/`@fr_files` list as one stream of frames.
struct Sequence {
    files: Vec<String>,
    next: usize,
    reader: Option<TextReader>,
    kind: &'static str,
    lib_name: String,
}

impl Sequence {
    fn new(files: Vec<String>, kind: &'static str, lib_name: &str) -> Self {
        Sequence {
            files,
            next: 0,
            reader: None,
            kind,
            lib_name: lib_name.to_string(),
        }
    }

    fn is_empty(&self) -> bool {
        self.files.is_empty()
    }

    /// Read the next frame, opening the next file when the current one runs out.
    fn read_frame(&mut self, etrj: &mut EnergyTraj) -> Result<bool, String> {
        loop {
            if self.reader.is_none() {
                let Some(path) = self.files.get(self.next) else {
                    return Ok(false);
                };
                self.next += 1;
                let mut reader = TextReader::open(path)
                    .map_err(|e| format!("{}: cannot open {path}: {e}", self.kind))?;
                let preamble = read_preamble(&mut reader).map_err(|e| format!("{path}: {e}"))?;
                etrj.bind(self.kind, &preamble, &self.lib_name, path)?;
                for w in etrj.drain_warnings() {
                    eprintln!("{w}");
                }
                self.reader = Some(reader);
            }
            let reader = self.reader.as_mut().unwrap();
            let file = &self.files[self.next - 1];
            let more = etrj
                .read_frame(reader, self.kind)
                .map_err(|e| format!("{file}: {e}"))?;
            for w in etrj.drain_warnings() {
                eprintln!("{w}");
            }
            if more {
                return Ok(true);
            }
            self.reader = None;
        }
    }
}

/// `@print_library [trajectory]`: the library for the built-in layout, or for the layout a
/// trajectory declares for itself (its `VARIABLES` are md++'s either way).
fn print_library(trajectory: Option<&str>) -> Result<String, String> {
    let Some(path) = trajectory else {
        return Ok(OFFICIAL_LIBRARY.to_string());
    };
    let mut reader = TextReader::open(path).map_err(|e| format!("{path}: {e}"))?;
    let pre = read_preamble(&mut reader).map_err(|e| format!("{path}: {e}"))?;
    let version = pre.version.clone().unwrap_or_default();
    let (schemas, origin) = match (&pre.schema, pre.version.as_deref()) {
        (Some(schema), _) => (vec![schema.clone()], "the self-description of"),
        (None, Some(v)) if builtin_schema(v, "ENERTRJ").is_some() => (
            ["ENERTRJ", "FRENERTRJ"]
                .iter()
                .filter_map(|t| builtin_schema(v, t))
                .collect(),
            "the built-in layout for the ENEVERSION of",
        ),
        _ => {
            return Err(format!(
                "{path}: no self-description, and ENEVERSION {} has no built-in layout",
                pre.version.as_deref().unwrap_or("NONE")
            ))
        },
    };
    let title = format!(
        "  ene_ana library generated from {origin} {path}\n  \
         schema from the trajectory; VARIABLES from md++"
    );
    Ok(library_text(
        &title,
        &version,
        &schemas,
        &official_variables(),
    ))
}

fn run() -> Result<(), String> {
    let args = Arguments::parse(
        &[
            "topo",
            "time",
            "en_files",
            "fr_files",
            "prop",
            "library",
            "print_library",
        ],
        USAGE,
    )?;

    if args.has("print_library") {
        print!(
            "{}",
            print_library(args.values("print_library").first().map(String::as_str))?
        );
        return Ok(());
    }

    // gromos++ only switches to a user time when both t0 and dt are given.
    let tv = args.values("time");
    let mut t0: f64 = tv
        .first()
        .map(|s| s.parse().map_err(|_| "@time: t is not numeric".to_string()))
        .transpose()?
        .unwrap_or(0.0);
    let dt: f64 = tv
        .get(1)
        .map(|s| {
            s.parse::<f64>()
                .map_err(|_| "@time: dt is not numeric".to_string())
        })
        .transpose()?
        .unwrap_or(1.0);
    let usertime = tv.len() >= 2;
    if usertime {
        t0 -= dt;
    }

    let en_files: Vec<String> = args.values("en_files").to_vec();
    let fr_files: Vec<String> = args.values("fr_files").to_vec();
    if en_files.is_empty() && fr_files.is_empty() {
        return Err(format!("no data specified:\n{USAGE}"));
    }
    let props: Vec<String> = args.values("prop").to_vec();
    if props.is_empty() {
        return Err(format!("no properties to follow:\n{USAGE}"));
    }
    let library = args.values("library");
    let print_library = library.iter().any(|s| s == "print");
    let lib_path = library.iter().find(|s| *s != "print");
    let (lib_name, lib_text) = match lib_path {
        Some(path) => {
            let text = std::io::read_to_string(
                gromos_io::gz::open_text(path)
                    .map_err(|e| format!("failed to open library file {path}: {e}"))?,
            )
            .map_err(|e| format!("{path}: {e}"))?;
            (path.clone(), text)
        },
        None => (
            "<built-in library>".to_string(),
            OFFICIAL_LIBRARY.to_string(),
        ),
    };
    let mut etrj = EnergyTraj::from_library(&lib_text).map_err(|e| format!("{lib_name}: {e}"))?;

    // gromos++ takes MASS and NUMMOL from the topology, for densit/Hvap and friends.
    if args.has("topo") {
        let path = args.value("topo")?;
        let topo = build_topology(read_topology_file(path).map_err(|e| format!("@topo: {e}"))?);
        etrj.add_constant("NUMMOL", topo.molecules.len() as f64);
        let mass: f64 = (0..topo.num_solute_atoms())
            .filter_map(|i| topo.mass.get(i))
            .sum();
        etrj.add_constant("MASS", mass);
    } else {
        etrj.add_constant("MASS", 0.0);
    }

    if print_library {
        for name in etrj.variable_names() {
            println!("{name}");
        }
    }

    let mut en = Sequence::new(en_files, "ENERTRJ", &lib_name);
    let mut fr = Sequence::new(fr_files, "FRENERTRJ", &lib_name);
    let do_en = !en.is_empty();
    let do_fr = !fr.is_empty();

    let mut stats: Vec<Stat> = vec![Stat::new(); props.len()];
    let mut times: Vec<f64> = Vec::new();

    loop {
        if do_en && !en.read_frame(&mut etrj)? {
            if do_fr && fr.read_frame(&mut etrj)? {
                eprintln!("# WARNING: frames left over in free energy trajectory,");
                eprintln!("#   they will be disregarded.");
            }
            break;
        }
        if do_fr && !fr.read_frame(&mut etrj)? {
            if do_en {
                eprintln!("# WARNING: frames left over in energy trajectory,");
                eprintln!("#   they will be disregarded.");
            }
            break;
        }
        for (s, p) in stats.iter_mut().zip(&props) {
            s.add(etrj.value(p)?);
        }
        if usertime {
            t0 += dt;
        } else {
            t0 = etrj.value("TIME[2]")?;
        }
        times.push(t0);
    }

    if stats.iter().any(|s| s.ee_strict().is_nan()) {
        eprintln!("# WARNING: One of the values is a NaN,");
        eprintln!("#   the data provided are not enough to ");
        eprintln!("#   give a sensible error estimate");
    }

    let mut out = String::new();
    let _ = writeln!(
        out,
        "{:>10} {:>15} {:>15} {:>15}",
        "property", "average", "rmsd", "error est."
    );
    for (s, p) in stats.iter().zip(&props) {
        let mut series = String::new();
        let _ = writeln!(series, "#{:>14} {:>15}", "time", p);
        for (t, v) in times.iter().zip(s.values()) {
            let _ = writeln!(series, "{} {}", num(*t), num(*v));
        }
        std::fs::write(format!("{p}.dat"), series).map_err(|e| format!("{p}.dat: {e}"))?;

        let _ = writeln!(
            out,
            "{:>10} {} {} {}",
            p,
            num(s.ave()),
            num(s.rmsd_strict()),
            num(s.ee_strict())
        );
    }
    print!("{out}");
    Ok(())
}
