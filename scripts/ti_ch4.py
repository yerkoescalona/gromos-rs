#!/usr/bin/env python3
"""
Thermodynamic integration of CH4 -> dummy in SPC water (the vol. 7 tutorial system), through
gromosXX and gromos-rs, from the reference inputs that are proven bit-for-bit at five λ values
(`crates/gromos-md/tests/gromosXX_references/ch4_water_fep*`).

For every λ window the same `.imd` is run by both engines (NVT 300 K, velocities generated,
dH/dλ written every `--ntwg` steps); the windows are integrated with the project's `ext_ti_ana`
(trapezoid over ⟨∂H/∂λ⟩_λ, error from the per-window standard error) and the two ΔG values are
compared. ΔG(CH4 -> dummy) is the negative hydration free energy: experiment ΔG_hyd(CH4) ≈
+8.4 kJ/mol, so ΔG here is expected around −8 kJ/mol for 54a7 (the tutorial's own result).

    python3 scripts/ti_ch4.py --steps 5000                # 10 ps per window, 11 windows, both engines
    python3 scripts/ti_ch4.py --engine gromos-rs --steps 20000 --lambdas 0,0.25,0.5,0.75,1

Needs `cargo build --release --bin md --bin ext_ti_ana` and (for gromosXX) the native build at
`.local/gromosXX/md++/BUILD/program/md`. Re-running skips windows whose `.trg` already exists.
Writes `<workdir>/README.md` (the table) and `<workdir>/results.json`.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
REF = REPO / "crates/gromos-md/tests/gromosXX_references"
TEMPLATE = REF / "ch4_water_fep/ch4_water_fep.in"
TOPO = REF / "shared/ch4_spc.top"
CONF = REF / "shared/ch4_spc.cnf"
PTTOPO = REF / "shared/ch4_spc_dummy.ptp"
ENGINES = {
    "gromosxx": REPO / ".local/gromosXX/md++/BUILD/program/md",
    "gromos-rs": REPO / "target/release/md",
}
EXT_TI_ANA = REPO / "target/release/ext_ti_ana"


def make_imd(lam: float, steps: int, ntwe: int, ntwg: int, seed: int) -> str:
    s = TEMPLATE.read_text()
    subs = [
        (r"(#   NSTLIM         T        DT\n)\s+\d+\s+0\.0\s+0\.002", rf"\g<1>{steps:>9}       0.0     0.002"),
        (r"(#  NTIVEL  NTISHK  NTINHT  NTINHB\n)\s+\d", r"\g<1>        1"),
        (r"(#      IG   TEMPI\n)\s+\d+\s+[\d.]+", rf"\g<1>{seed:>9}   300.0"),
        (r"(#   TEMP0  TAU\n    300\.0)\s+-?[\d.]+", r"\g<1>  0.1"),
        (r"(#     NTPR      NTPP\n)\s+\d+", r"\g<1>      1000"),
        (r"(#     NTWX     NTWSE      NTWV      NTWF      NTWE      NTWG      NTWB\n)[ \d]+",
         rf"\g<1>         0         0         0         0{ntwe:>10}{ntwg:>10}         0"),
        (r"(standard\s+)1(\s+0\.8\s+1\.4)", r"\g<1>5\g<2>"),
        (r"(#     NTG   NRDGL    RLAM   DLAMT\n        1       0)\s+[\d.]+", rf"\g<1>   {lam:.3f}"),
    ]
    for pat, rep in subs:
        s, n = re.subn(pat, rep, s)
        assert n == 1, f"template edit failed: {pat}"
    s = s.replace("One step for reference energy validation.", f"TI window lambda={lam:.3f}, {steps} steps, NVT 300 K")
    return s


def run_window(engine: str, lam: float, wdir: Path, imd: Path) -> tuple[Path, float]:
    tag = f"L_{lam:.3f}"
    out = wdir / engine / tag
    out.mkdir(parents=True, exist_ok=True)
    trg = out / "window.trg"
    if trg.exists() and trg.stat().st_size > 1000:
        return trg, 0.0
    cmd = [str(ENGINES[engine]), "@topo", str(TOPO), "@pttopo", str(PTTOPO), "@conf", str(CONF), "@input", str(imd),
           "@fin", str(out / "final.cnf"), "@tre", str(out / "window.tre"), "@trg", str(trg)]
    t0 = time.perf_counter()
    with open(out / "run.log", "w") as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    dt = time.perf_counter() - t0
    if rc != 0 or not trg.exists():
        raise SystemExit(f"{engine} failed on {tag}: see {out / 'run.log'}")
    return trg, dt


def analyse(trgs: list[Path], lambdas: list[float], skip: int) -> dict:
    cmd = [str(EXT_TI_ANA), "@trg", *map(str, trgs), "@lambda", *[f"{l:.6f}" for l in lambdas], "@skip", str(skip)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise SystemExit(f"ext_ti_ana failed:\n{res.stderr}")
    windows = []
    for line in res.stdout.splitlines():
        m = re.match(r"\s+([\d.]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+(\d+)\s*$", line)
        if m:
            windows.append({"lambda": float(m[1]), "mean": float(m[2]), "ee": float(m[3]), "n": int(m[4])})
    dg = re.search(r"ΔG = ([-\d.]+) kJ/mol", res.stdout)
    err = re.search(r"±\s*([\d.e+-]+)", res.stdout[res.stdout.find("TI integration"):] if dg else "")
    return {"windows": windows, "dG": float(dg[1]) if dg else None, "dG_err": float(err[1]) if err else None,
            "stdout": res.stdout}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", choices=["both", *ENGINES], default="both")
    ap.add_argument("--steps", type=int, default=5000, help="MD steps per window (dt = 2 fs)")
    ap.add_argument("--lambdas", default=",".join(f"{i / 10:.1f}" for i in range(11)))
    ap.add_argument("--ntwg", type=int, default=10, help="dH/dlambda output frequency (steps)")
    ap.add_argument("--ntwe", type=int, default=100, help="energy output frequency (steps)")
    ap.add_argument("--skip", type=int, default=None, help="dH/dlambda frames to drop per window (default: 20 %%)")
    ap.add_argument("--seed", type=int, default=210185)
    ap.add_argument("--workdir", type=Path, default=REPO / "bench/work/ti_ch4")
    args = ap.parse_args()

    lambdas = [float(x) for x in args.lambdas.split(",")]
    engines = list(ENGINES) if args.engine == "both" else [args.engine]
    for e in engines:
        if not ENGINES[e].exists():
            raise SystemExit(f"{e}: binary not found at {ENGINES[e]}")
    if not EXT_TI_ANA.exists():
        raise SystemExit(f"ext_ti_ana not built: cargo build --release -p gromos-analysis --bin ext_ti_ana")
    n_frames = args.steps // args.ntwg
    skip = args.skip if args.skip is not None else n_frames // 5
    wdir = args.workdir
    wdir.mkdir(parents=True, exist_ok=True)

    results: dict = {"steps": args.steps, "ntwg": args.ntwg, "skip_frames": skip, "lambdas": lambdas, "engines": {}}
    for engine in engines:
        trgs, timings = [], {}
        for i, lam in enumerate(lambdas):
            imd = wdir / f"ti_L_{lam:.3f}.imd"
            imd.write_text(make_imd(lam, args.steps, args.ntwe, args.ntwg, args.seed + i))
            trg, dt = run_window(engine, lam, wdir, imd)
            trgs.append(trg)
            timings[f"{lam:.3f}"] = dt
            print(f"{engine:10s} λ={lam:.3f}  {'cached' if dt == 0 else f'{dt:7.1f} s'}", flush=True)
        ana = analyse(trgs, lambdas, skip)
        results["engines"][engine] = {"timings_s": timings, **{k: v for k, v in ana.items() if k != "stdout"}}
        (wdir / f"ext_ti_ana_{engine}.txt").write_text(ana["stdout"])

    # ---- table ----
    lines = [f"# TI: CH4 -> dummy in 999 SPC (54a7), {args.steps} steps/window (dt 2 fs), "
             f"{len(lambdas)} windows, dH/dλ every {args.ntwg} steps, first {skip} frames dropped", "",
             "| λ | " + " | ".join(f"⟨∂H/∂λ⟩ {e} ± ee (n)" for e in engines) + (" | Δ (rs − XX) |" if len(engines) == 2 else " |"),
             "|---|" + "---|" * len(engines) + ("---|" if len(engines) == 2 else "")]
    per = {e: {round(w["lambda"], 6): w for w in results["engines"][e]["windows"]} for e in engines}
    for lam in lambdas:
        cells = []
        for e in engines:
            w = per[e].get(round(lam, 6))
            cells.append(f"{w['mean']:9.3f} ± {w['ee']:.3f} ({w['n']})" if w else "—")
        row = f"| {lam:.2f} | " + " | ".join(cells)
        if len(engines) == 2:
            a, b = per[engines[0]].get(round(lam, 6)), per[engines[1]].get(round(lam, 6))
            row += f" | {b['mean'] - a['mean']:+.3f} |" if a and b else " | — |"
        else:
            row += " |"
        lines.append(row)
    lines.append("")
    for e in engines:
        r = results["engines"][e]
        tot = sum(r["timings_s"].values())
        if r["dG"] is None:
            lines.append(f"- **{e}**: no ΔG (ext_ti_ana found no usable frames — see ext_ti_ana_{e}.txt)")
            continue
        err = f" ± {r['dG_err']:.2f}" if r.get("dG_err") is not None else ""
        lines.append(f"- **{e}**: ΔG = {r['dG']:.2f}{err} kJ/mol; wall time {tot:.0f} s for {len(lambdas)} windows"
                     + ("" if tot else " (all cached)"))
    if len(engines) == 2 and all(results["engines"][e]["dG"] is not None for e in engines):
        a, b = (results["engines"][e]["dG"] for e in engines)
        lines.append(f"- ΔG(gromos-rs) − ΔG(gromosXX) = {b - a:+.3f} kJ/mol")
    lines.append("- Reference: ΔG_hyd(CH4) experiment ≈ +8.4 kJ/mol, i.e. ΔG(CH4→dummy) ≈ −8 kJ/mol; the "
                 "tutorial's 54a7 result is of that size. Windows are independent short runs, not a converged study.")
    table = "\n".join(lines) + "\n"
    (wdir / "README.md").write_text(table)
    (wdir / "results.json").write_text(json.dumps(results, indent=2))
    print("\n" + table)


if __name__ == "__main__":
    main()
