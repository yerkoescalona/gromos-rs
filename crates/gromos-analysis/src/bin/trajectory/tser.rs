//! tser — time series of properties (gromos++ `tser`).
//!
//! Usage: tser @topo <top> @pbc <r|v> [cog|nog] [@time <t0 dt>] @prop <properties> @traj <trc…>
//!             [@nots] [@dist <steps [min max] [periodic]>] [@norm] [@solv] [@skip <n>] [@stride <n>]
//!
//! Properties are gromos++ property specifiers (`gromos_analysis::property`): `d%1:1,2`,
//! `a%1:1,2,3`, `t%1:1,2,3,4`, `tp%…`, `o%atom(…)%cart(…)`, `op%…`, `pr%…`, `pa%…`.

use gromos_analysis::args::{fail, Arguments};
use gromos_analysis::distribution::{trim_float, Distribution};
use gromos_analysis::pbc::Pbc;
use gromos_analysis::property::PropertyContainer;
use gromos_analysis::time::Time;
use gromos_io::topology::{build_topology, read_topology_file};
use gromos_io::trajectory::TrajectoryReader;

const USAGE: &str = "# tser
\t@topo      <molecular topology file>
\t@pbc       <boundary type> [<gathermethod>]
\t@time      <time and dt>
\t@prop      <properties>
\t[@nots     (do not write time series)]
\t[@dist     <steps [min max] [periodic]>]
\t[@norm     (normalise distribution)]
\t[@solv     (read in solvent)]
\t@traj      <trajectory files>
\t[@skip     <skip n first frames>]
\t[@stride   <take every n-th frame>]";

fn main() {
    if let Err(e) = run() {
        fail(&e);
    }
}

fn run() -> Result<(), String> {
    let args = Arguments::parse(
        &[
            "topo", "pbc", "time", "prop", "traj", "skip", "stride", "nots", "dist", "norm", "solv",
        ],
        USAGE,
    )?;
    let mut time = Time::from_args(&args)?;
    let do_tser = args.count("nots") < 0;
    let normalize = args.count("norm") >= 0;
    let dist = args.values("dist");
    let do_dist = !dist.is_empty();
    let mut dist_steps = 0usize;
    let mut bounds: Option<(f64, f64)> = None;
    let mut periodic = false;
    if do_dist {
        dist_steps = dist[0]
            .parse::<usize>()
            .map_err(|_| format!("distribution: wrong number of steps specified {}", dist[0]))?;
        if dist_steps == 0 {
            return Err(format!(
                "distribution: wrong number of steps specified {}",
                dist[0]
            ));
        }
        if dist.len() >= 3 {
            let lo: f64 = dist[1]
                .parse()
                .map_err(|_| "distribution: bad minimum".to_string())?;
            let hi: f64 = dist[2]
                .parse()
                .map_err(|_| "distribution: bad maximum".to_string())?;
            if hi < lo {
                return Err("distribution: maximum value smaller than minimum value".into());
            }
            bounds = Some((lo, hi));
        }
        if let Some(kw) = dist.get(3) {
            if kw == "periodic" {
                periodic = true;
            } else {
                return Err(format!("distribution: keyword {kw} not known"));
            }
        }
    }
    let topo =
        build_topology(read_topology_file(args.value("topo")?).map_err(|e| format!("@topo: {e}"))?);
    let pbc = Pbc::from_args(&args)?;
    let mut props = PropertyContainer::parse(args.values("prop"), &topo)?;
    let skip: usize = args.get("skip", 0)?;
    let stride: usize = args.get("stride", 1)?.max(1);
    let mut out = String::new();
    if do_tser {
        out.push_str(&format!("#\n#{:>9}\t\t{}\n", "time", props.title(&topo)));
    }
    let mut frame_no = 0usize;
    for traj in args.values("traj") {
        let mut reader = TrajectoryReader::new(traj).map_err(|e| format!("@traj {traj}: {e}"))?;
        while let Some(frame) = reader.read_frame().map_err(|e| format!("{traj}: {e}"))? {
            let use_frame = frame_no >= skip && (frame_no - skip).is_multiple_of(stride);
            frame_no += 1;
            if !use_frame {
                continue;
            }
            let t = time.next(frame.time);
            let periodicity = pbc.periodicity(frame.box_dims);
            let mut pos = frame.positions;
            pbc.gather(&topo, &mut pos, &periodicity);
            let values = props.calc(&pos, &periodicity);
            if do_tser {
                out.push_str(&Time::fmt(t));
                out.push_str("\t\t");
                for v in &values {
                    out.push_str(&trim_float(*v));
                    out.push_str("\t\t");
                }
                out.push('\n');
            }
        }
    }
    if do_tser {
        out.push_str("# Averages over run: (<average> <rmsd> <error estimate>)\n");
        for p in &props.props {
            out.push_str(&format!("# {}\t{}\n", p.title(&topo), p.average()));
        }
    }
    if do_dist {
        for p in &props.props {
            let values = p.stat.values();
            let d = match bounds {
                Some((lo, hi)) => {
                    Distribution::over_values_in(values, lo, hi, dist_steps, periodic)?
                },
                None => Distribution::over_values(values, dist_steps)?,
            };
            out.push_str("\n#\n");
            out.push_str(&format!("# Distribution of      {}\n", p.title(&topo)));
            out.push_str(&format!("# values:              {}\n", p.stat.n()));
            out.push_str(&format!("# average value:       {:.6}\n", p.stat.ave()));
            out.push_str(&format!("# rmsd:                {:.6}\n", p.stat.rmsd()));
            let ee = p.stat.ee_strict();
            out.push_str(&format!(
                "# error estimate:      {}\n",
                if ee.is_nan() {
                    "-nan".to_string()
                } else {
                    format!("{:.6}", ee)
                }
            ));
            out.push_str(&format!("# maximum value at:    {:.6}\n", d.max_val_at()));
            let fixed6 = |x: f64| format!("{:.6}", x);
            out.push_str(&if normalize {
                d.write_normalized(&fixed6)
            } else {
                d.write(&fixed6)
            });
        }
    }
    print!("{out}");
    Ok(())
}
