#!/usr/bin/env python3
"""Regenerate gromosXX reference data using the patched CODATA binary.

Run from the repo root:
    python3 scripts/regen_gromosXX_references.py

The script runs the gromosXX md binary (patched to use CODATA 2018 constants)
on every system in crates/gromos-md/tests/gromosXX_references/ and overwrites
the expected/ output files in place.

Requires the patched binary at:
    .local/gromosXX/md++/build/program/md
"""

import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
REF_DIR = REPO_ROOT / "crates" / "gromos-md" / "tests" / "gromosXX_references"
MD_BIN = REPO_ROOT / ".local" / "gromosXX" / "md++" / "build" / "program" / "md"


def parse_toml_str(text: str, key: str) -> str | None:
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(f"{key} ="):
            parts = line.split("=", 1)[1].strip()
            if parts.startswith('"') and parts.endswith('"'):
                return parts[1:-1]
    return None


def run_system(sys_dir: Path) -> bool:
    name = sys_dir.name
    toml_path = sys_dir / "input.toml"
    if not toml_path.exists():
        print(f"  SKIP {name}: no input.toml")
        return True

    toml = toml_path.read_text()

    topo = sys_dir / parse_toml_str(toml, "topology")
    conf = sys_dir / parse_toml_str(toml, "configuration")
    params = sys_dir / parse_toml_str(toml, "parameters")

    with tempfile.TemporaryDirectory(prefix=f"regen_{name}_") as tmp:
        tmp = Path(tmp)
        tre = tmp / "energies.tre"
        trf = tmp / "forces.trf"
        trc = tmp / "trajectory.trc"
        trg = tmp / "free_energy.trg"
        fin = tmp / "final.conf"

        cmd = [
            str(MD_BIN),
            "@topo", str(topo),
            "@conf", str(conf),
            "@input", str(params),
            "@fin", str(fin),
            "@tre", str(tre),
            "@trf", str(trf),
            "@trc", str(trc),
            "@trg", str(trg),
        ]

        posresspec = parse_toml_str(toml, "posresspec")
        if posresspec:
            cmd += ["@posresspec", str(sys_dir / posresspec)]

        refpos = parse_toml_str(toml, "refpos")
        if refpos:
            cmd += ["@refpos", str(sys_dir / refpos)]

        distrest = parse_toml_str(toml, "distrest")
        if distrest:
            cmd += ["@distrest", str(sys_dir / distrest)]

        pttopo = parse_toml_str(toml, "pttopo")
        if pttopo:
            cmd += ["@pttopo", str(sys_dir / pttopo)]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"  FAIL {name}: md exited {result.returncode}")
            print(result.stderr[-2000:])
            return False

        # Copy outputs into expected/
        expected = sys_dir / "expected"
        expected.mkdir(exist_ok=True)

        copied = []
        for src, dst_name in [
            (tre, "energies.tre"),
            (trf, "forces.trf"),
            (trc, "trajectory.trc"),
            (fin, "final.conf"),
        ]:
            if src.exists():
                shutil.copy2(src, expected / dst_name)
                copied.append(dst_name)

        fe_path = parse_toml_str(toml, "free_energy")
        if fe_path and trg.exists():
            dst = sys_dir / fe_path
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(trg, dst)
            copied.append(fe_path)

        print(f"  OK  {name}: {', '.join(copied)}")
        return True


def main():
    if not MD_BIN.exists():
        print(f"ERROR: gromosXX binary not found at {MD_BIN}")
        sys.exit(1)

    print(f"Using binary: {MD_BIN}")
    print(f"Reference dir: {REF_DIR}")
    print()

    systems = sorted(p for p in REF_DIR.iterdir() if p.is_dir() and p.name != "shared")

    failed = []
    for sys_dir in systems:
        ok = run_system(sys_dir)
        if not ok:
            failed.append(sys_dir.name)

    print()
    print(f"Done: {len(systems) - len(failed)}/{len(systems)} systems regenerated.")
    if failed:
        print(f"FAILED: {', '.join(failed)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
