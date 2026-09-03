//! ene_ana — averages and time series of properties from energy trajectories
//! (gromos++ `ene_ana`, `programs/ene_ana.cc`).
//!
//! The layout of a `.tre`/`.trg` is not built in: `@library` declares it and names the
//! properties, which is what lets one library file span GROMOS versions. See
//! `gromos_io::energy_traj`.
//!
//! Each `@prop` is written as a `<prop>.dat` time series and summarised on stdout as
//! average / rmsd / error estimate. The error estimate is gromos++'s block-averaging one,
//! which needs more than 200 frames — with fewer, gromos++ prints `-nan` and so do we.

use gromos_analysis::args::{fail, Arguments};
use gromos_core::stat::Stat;
use gromos_io::energy_traj::EnergyTraj;
use gromos_io::gz::TextReader;
use gromos_io::table::cpp_g;
use gromos_io::topology::{build_topology, read_topology_file};
use std::fmt::Write as _;

const USAGE: &str = "# ene_ana
\t@en_files    <energy files> (and/or)
\t@fr_files    <free energy files>
\t@prop        <properties to monitor>
\t@library     <library for property names> [print]
\t[@topo       <molecular topology file> (for MASS and NUMMOL)]
\t[@time       <t and dt> (overwrites TIME in the trajectory files)]";

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
}

impl Sequence {
    fn new(files: Vec<String>, kind: &'static str) -> Self {
        Sequence {
            files,
            next: 0,
            reader: None,
            kind,
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
                skip_preamble(&mut reader, path)?;
                self.reader = Some(reader);
            }
            let reader = self.reader.as_mut().unwrap();
            let file = &self.files[self.next - 1];
            if etrj
                .read_frame(reader, self.kind)
                .map_err(|e| format!("{file}: {e}"))?
            {
                return Ok(true);
            }
            self.reader = None;
        }
    }
}

/// A trajectory opens with a `TITLE` and may carry an `ENEVERSION`; frames start after them.
/// gromos++ consumes both in `Ginstream::open` before the first `read_frame`.
fn skip_preamble(reader: &mut TextReader, path: &str) -> Result<Option<String>, String> {
    let mut version = None;
    let mut line = String::new();
    loop {
        let mut peeked = String::new();
        if reader.read_line(&mut peeked).map_err(|e| e.to_string())? == 0 {
            return Ok(version);
        }
        let name = peeked.split('#').next().unwrap_or("").trim().to_string();
        if name.is_empty() {
            continue;
        }
        if name != "TITLE" && name != "ENEVERSION" {
            reader.unread_line(&peeked);
            return Ok(version);
        }
        let mut body = String::new();
        loop {
            if reader.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
                return Err(format!("{path}: no END in {name} block"));
            }
            let l = line.split('#').next().unwrap_or("").trim();
            if l == "END" {
                break;
            }
            body.push_str(l);
        }
        if name == "ENEVERSION" {
            version = Some(body.split_whitespace().collect());
        }
    }
}

fn run() -> Result<(), String> {
    let args = Arguments::parse(
        &["topo", "time", "en_files", "fr_files", "prop", "library"],
        USAGE,
    )?;

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
    let Some(lib_path) = library.first() else {
        return Err(format!("no library file specified:\n{USAGE}"));
    };
    let print_library = library.get(1).is_some_and(|s| s == "print");

    let lib_text = std::io::read_to_string(
        gromos_io::gz::open_text(lib_path)
            .map_err(|e| format!("failed to open library file {lib_path}: {e}"))?,
    )
    .map_err(|e| format!("{lib_path}: {e}"))?;
    let mut etrj = EnergyTraj::from_library(&lib_text)?;

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

    let mut en = Sequence::new(en_files, "ENERTRJ");
    let mut fr = Sequence::new(fr_files, "FRENERTRJ");
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
