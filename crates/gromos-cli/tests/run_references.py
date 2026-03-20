#!/usr/bin/env python3
"""
GROMOS-RS Reference Data Generator
===================================
Runs gromosXX (md++) on the reference systems and saves raw output files
(.trf, .tre, .trc) alongside a TOML manifest for direct use in Rust tests.

Usage:
    python3 run_references.py --md-binary /path/to/md

Output structure per system (self-contained, copy to gromos-rs/tests/data/):
    pair_lj/
    ├── input.toml          # Test metadata + key reference values + tolerances
    ├── pair_lj.topo        # Input (already exists)
    ├── pair_lj.conf        # Input (already exists)
    ├── pair_lj.in          # Input (already exists)
    └── expected/
        ├── forces.trf      # Raw gromosXX force trajectory (full double precision)
        ├── energies.tre    # Raw gromosXX energy trajectory
        ├── trajectory.trc  # Raw gromosXX coordinate trajectory
        └── final.conf      # Final configuration after N steps

Precision:
    - Internal: double (64-bit IEEE 754)
    - Output positions: 9 decimal digits (fixed-point)
    - Output forces: 9 decimal digits (fixed-point)
    - Output energies: 10 significant digits (scientific notation)
    This preserves ~53 bits of mantissa — effectively full double.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
REFERENCES_DIR = SCRIPT_DIR / "gromosXX_references"

# Systems in pyramid order (simplest first)
SYSTEMS = [
    {
        "name": "pair_lj",
        "description": "Level 0: Two Ar atoms, pure LJ, no PBC",
        "level": 0,
        "dir": "pair_lj",
        "topo": "pair_lj.topo",
        "conf": "pair_lj.conf",
        "input": "pair_lj.in",
        "pttopo": None,
        "natoms": 2,
        "isolates": "Lennard-Jones pair force",
    },
    {
        "name": "water_single",
        "description": "Level 1: Single SPC water, bond+angle+coulomb",
        "level": 1,
        "dir": "water_single",
        "topo": "water_single.topo",
        "conf": "water_single.conf",
        "input": "water_single.in",
        "pttopo": None,
        "natoms": 3,
        "isolates": "bond + angle + intramolecular Coulomb",
    },
    {
        "name": "aladip_vacuum",
        "description": "Level 1: Aladip in vacuum, all bonded + nonbonded intra",
        "level": 1,
        "dir": "aladip_vacuum",
        "topo": "../../aladip.topo",
        "conf": "aladip_vacuum.conf",
        "input": "aladip_vacuum.in",
        "pttopo": None,
        "natoms": 12,
        "isolates": "all bonded terms + exclusions + 1-4 interactions",
    },
    {
        "name": "water_3_box",
        "description": "Level 2: 3 waters in PBC box, tests minimum image",
        "level": 2,
        "dir": "water_3_box",
        "topo": "water_3_box.topo",
        "conf": "water_3_box.conf",
        "input": "water_3_box.in",
        "pttopo": None,
        "natoms": 9,
        "isolates": "PBC + minimum image + pairlist + intermolecular CRF",
    },
    {
        "name": "water_216_box",
        "description": "Level 3: 216 waters bulk, scaling + long-range",
        "level": 3,
        "dir": "water_216_box",
        "topo": "water_216_box.topo",
        "conf": "water_216_box.conf",
        "input": "water_216_box.in",
        "pttopo": None,
        "natoms": 648,
        "isolates": "bulk nonbonded scaling + pairlist update + virial",
    },
    {
        "name": "aladip_solvated",
        "description": "Level 4: Aladip + 20 SPC, full system with SHAKE",
        "level": 4,
        "dir": "aladip_solvated",
        "topo": "../../aladip.topo",
        "conf": "../../aladip.conf",
        "input": "aladip_solvated.in",
        "pttopo": None,
        "natoms": 72,
        "isolates": "SHAKE constraints + solute-solvent + full MD",
    },
    # --- New systems: ions, unlike LJ, aromatics ---
    {
        "name": "nacl_pair",
        "description": "Level 0: Na+ Cl- pair in vacuum, full unit charges",
        "level": 0,
        "dir": "nacl_pair",
        "topo": "nacl_pair.topo",
        "conf": "nacl_pair.conf",
        "input": "nacl_pair.in",
        "pttopo": None,
        "natoms": 2,
        "isolates": "Coulomb with full unit charges + LJ ion params",
    },
    {
        "name": "pair_lj_mixed",
        "description": "Level 0: Ar-CH4 unlike LJ pair, combination rules",
        "level": 0,
        "dir": "pair_lj_mixed",
        "topo": "pair_lj_mixed.topo",
        "conf": "pair_lj_mixed.conf",
        "input": "pair_lj_mixed.in",
        "pttopo": None,
        "natoms": 2,
        "isolates": "LJ combination rules (unlike atom types)",
    },
    {
        "name": "nacl_water_box",
        "description": "Level 2: Na+ Cl- + 20 SPC waters in PBC",
        "level": 2,
        "dir": "nacl_water_box",
        "topo": "nacl_water_box.topo",
        "conf": "nacl_water_box.conf",
        "input": "nacl_water_box.in",
        "pttopo": None,
        "natoms": 62,
        "isolates": "ion-water RF + monatomic charge groups in PBC",
    },
    {
        "name": "benzene_vacuum",
        "description": "Level 1: Benzene all-atom in vacuum, aromatic ring",
        "level": 1,
        "dir": "benzene_vacuum",
        "topo": "benzene_vacuum.topo",
        "conf": "benzene_vacuum.conf",
        "input": "benzene_vacuum.in",
        "pttopo": None,
        "natoms": 12,
        "isolates": "aromatic ring + improper dihedral planarity + ring torsion",
    },
    # --- Ensemble tests: NVT, NPT ---
    {
        "name": "water_216_nvt",
        "description": "Level 3: 216 waters, weak-coupling thermostat (NVT)",
        "level": 3,
        "dir": "water_216_nvt",
        "topo": "../water_216_box/water_216_box.topo",
        "conf": "../water_216_box/water_216_box.conf",
        "input": "water_216_nvt.in",
        "pttopo": None,
        "natoms": 648,
        "isolates": "Berendsen thermostat + kinetic energy rescaling",
    },
    {
        "name": "water_216_npt",
        "description": "Level 3: 216 waters, thermostat + Berendsen barostat (NPT)",
        "level": 3,
        "dir": "water_216_npt",
        "topo": "../water_216_box/water_216_box.topo",
        "conf": "../water_216_box/water_216_box.conf",
        "input": "water_216_npt.in",
        "pttopo": None,
        "natoms": 648,
        "isolates": "Berendsen barostat + pressure coupling + box rescaling",
    },
    # --- Diagnostic systems for bug isolation ---
    {
        "name": "nacl_pair_box",
        "description": "Level 1: Na+ Cl- in PBC box, ion pair with reaction field",
        "level": 1,
        "dir": "nacl_pair_box",
        "topo": "nacl_pair_box.topo",
        "conf": "nacl_pair_box.conf",
        "input": "nacl_pair_box.in",
        "pttopo": None,
        "natoms": 2,
        "isolates": "Coulomb + LJ with PBC + reaction field (no solvent)",
    },
    {
        "name": "nacl_1water_box",
        "description": "Level 2: Na+ Cl- + 1 SPC water in PBC",
        "level": 2,
        "dir": "nacl_1water_box",
        "topo": "nacl_1water_box.topo",
        "conf": "nacl_1water_box.conf",
        "input": "nacl_1water_box.in",
        "pttopo": None,
        "natoms": 5,
        "isolates": "minimal solute-solvent with SHAKE solvent in PBC",
    },
    {
        "name": "nacl_3water_box",
        "description": "Level 2: Na+ Cl- + 3 SPC waters in PBC",
        "level": 2,
        "dir": "nacl_3water_box",
        "topo": "nacl_3water_box.topo",
        "conf": "nacl_3water_box.conf",
        "input": "nacl_3water_box.in",
        "pttopo": None,
        "natoms": 11,
        "isolates": "multiple solvent molecules + solute-solvent pairlist",
    },
    {
        "name": "water_3_box_twinrange",
        "description": "Level 2: 3 waters with twin-range cutoff (RCUTP < RCUTL)",
        "level": 2,
        "dir": "water_3_box_twinrange",
        "topo": "water_3_box_twinrange.topo",
        "conf": "water_3_box_twinrange.conf",
        "input": "water_3_box_twinrange.in",
        "pttopo": None,
        "natoms": 9,
        "isolates": "twin-range pairlist with NSNB=5 + short/long-range split",
    },
    {
        "name": "nacl_3water_cutoff",
        "description": "Level 2: Na+ Cl- + 3 SPC waters at cutoff boundaries",
        "level": 2,
        "dir": "nacl_3water_cutoff",
        "topo": "nacl_3water_cutoff.topo",
        "conf": "nacl_3water_cutoff.conf",
        "input": "nacl_3water_cutoff.in",
        "pttopo": None,
        "natoms": 11,
        "isolates": "CRF at short-range/long-range/outside cutoff boundaries",
    },
    {
        "name": "water_10_box",
        "description": "Level 2: 10 SPC waters + 2 ions, positions away from cutoff boundaries",
        "level": 2,
        "dir": "water_10_box",
        "topo": "water_10_box.topo",
        "conf": "water_10_box.conf",
        "input": "water_10_box.in",
        "pttopo": None,
        "natoms": 32,
        "isolates": "medium solvent system with no pairs at cutoff boundary",
    },
    {
        "name": "nacl_water_box_shifted",
        "description": "Level 2: nacl_water_box with perturbed positions off cutoff boundaries",
        "level": 2,
        "dir": "nacl_water_box_shifted",
        "topo": "nacl_water_box_shifted.topo",
        "conf": "nacl_water_box_shifted.conf",
        "input": "nacl_water_box_shifted.in",
        "pttopo": None,
        "natoms": 62,
        "isolates": "same topology as nacl_water_box but no pairs at cutoff boundary",
    },
]


def find_md_binary(hint=None):
    """Try to locate the md binary."""
    if hint and os.path.isfile(hint) and os.access(hint, os.X_OK):
        return hint

    gromosxx_root = SCRIPT_DIR.parents[3]
    candidates = [
        gromosxx_root / "build" / "program" / "md",
        gromosxx_root / "program" / "md",
        gromosxx_root / "build" / "md",
    ]

    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return str(c)

    result = shutil.which("md")
    if result:
        return result

    return None


def run_gromosxx(md_binary, system, work_dir, verbose=False):
    """Run gromosXX md on a system. Returns paths to output files or None on failure."""
    sys_dir = REFERENCES_DIR / system["dir"]

    topo_path = (sys_dir / system["topo"]).resolve()
    conf_path = (sys_dir / system["conf"]).resolve()
    input_path = (sys_dir / system["input"]).resolve()

    fin_path = work_dir / "final.conf"
    trc_path = work_dir / "trajectory.trc"
    trf_path = work_dir / "forces.trf"
    tre_path = work_dir / "energies.tre"

    cmd = [
        md_binary,
        "@topo",
        str(topo_path),
        "@conf",
        str(conf_path),
        "@input",
        str(input_path),
        "@fin",
        str(fin_path),
        "@trc",
        str(trc_path),
        "@trf",
        str(trf_path),
        "@tre",
        str(tre_path),
    ]

    if system["pttopo"]:
        pttopo_path = (sys_dir / system["pttopo"]).resolve()
        cmd.extend(["@pttopo", str(pttopo_path)])

    print(f"  Running: md @topo ...{system['topo']} @conf ...{system['conf']}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(work_dir),
        timeout=300,
    )

    # Save full md log (stdout + stderr) to work_dir for debugging
    log_path = work_dir / "md_output.log"
    with open(log_path, "w") as f:
        f.write(f"# gromosXX md output for {system['name']}\n")
        f.write(f"# Command: {' '.join(cmd)}\n")
        f.write(f"# Return code: {result.returncode}\n\n")
        f.write("=== STDOUT ===\n")
        f.write(result.stdout)
        f.write("\n=== STDERR ===\n")
        f.write(result.stderr)

    if result.returncode != 0:
        print(f"  ERROR: md returned {result.returncode}")
        print(f"  Full log saved to: {log_path}")
        if verbose:
            print(f"  --- STDOUT (last 80 lines) ---")
            for line in result.stdout.split("\n")[-80:]:
                print(f"    {line}")
            if result.stderr.strip():
                print(f"  --- STDERR ---")
                for line in result.stderr.split("\n")[-40:]:
                    print(f"    {line}")
        else:
            # Show relevant error lines
            for line in result.stdout.split("\n"):
                if "ERROR" in line or "error" in line or "WARNING" in line:
                    print(f"    {line.strip()}")
        return None

    return {
        "stdout": result.stdout,
        "stderr": result.stderr,
        "log": log_path,
        "fin": fin_path,
        "trc": trc_path,
        "trf": trf_path,
        "tre": tre_path,
    }


def extract_key_values(run_output):
    """Extract key numerical values from stdout for quick validation in TOML."""
    values = {}
    for line in run_output["stdout"].split("\n"):
        line = line.strip()
        # Match energy lines like "E_Total                :  -4.9387e-01"
        m = re.match(r"(E_\S+)\s+:\s+([-\d.eE+]+)", line)
        if m:
            key = m.group(1)
            val = m.group(2)
            # Only keep first occurrence (step 0)
            if key not in values:
                values[key] = val
    return values


def extract_step0_forces(trf_path):
    """Read first FREEFORCERED block from force trajectory."""
    forces = []
    with open(trf_path, "r") as f:
        in_first_frame = False
        in_forces = False
        past_first_timestep = False
        for line in f:
            stripped = line.strip()
            if stripped == "TIMESTEP":
                if past_first_timestep:
                    break  # only want first frame
                past_first_timestep = True
                in_first_frame = True
                continue
            if in_first_frame and stripped == "FREEFORCERED":
                in_forces = True
                continue
            if in_forces and stripped == "END":
                in_forces = False
                continue
            if in_forces:
                parts = stripped.split()
                if len(parts) == 3:
                    forces.append((parts[0], parts[1], parts[2]))
    return forces


def generate_toml(system, run_output, key_values, step0_forces):
    """Generate the input.toml manifest for a test system."""
    lines = []
    lines.append(f"# Auto-generated by run_references.py")
    lines.append(f"# gromosXX v1.6.0, double precision, {system['natoms']} atoms")
    lines.append(f"")
    lines.append(f"[system]")
    lines.append(f'name = "{system["name"]}"')
    lines.append(f'description = "{system["description"]}"')
    lines.append(f"level = {system['level']}")
    lines.append(f"natoms = {system['natoms']}")
    lines.append(f'isolates = "{system["isolates"]}"')
    lines.append(f"")
    lines.append(f"[input]")
    lines.append(f'topology = "{system["topo"]}"')
    lines.append(f'configuration = "{system["conf"]}"')
    lines.append(f'parameters = "{system["input"]}"')
    if system["pttopo"]:
        lines.append(f'perturbation_topology = "{system["pttopo"]}"')
    lines.append(f"")
    lines.append(f"[expected]")
    lines.append(f'forces = "expected/forces.trf"')
    lines.append(f'energies = "expected/energies.tre"')
    lines.append(f'trajectory = "expected/trajectory.trc"')
    lines.append(f'final_conf = "expected/final.conf"')
    lines.append(f"")
    lines.append(f"[precision]")
    lines.append(f"# gromosXX internal: double (64-bit)")
    lines.append(f"# Output positions: 9 decimal digits fixed-point")
    lines.append(f"# Output forces: 9 decimal digits fixed-point")
    lines.append(f"# Output energies: 10 significant digits scientific")
    lines.append(f'internal = "f64"')
    lines.append(f"position_digits = 9")
    lines.append(f"force_digits = 9")
    lines.append(f"energy_sigfigs = 10")
    lines.append(f"")
    lines.append(f"[tolerances]")
    lines.append(f"# Recommended tolerances for gromos-rs comparison")
    lines.append(f"force_abs = 1e-6        # kJ/(mol*nm)")
    lines.append(f"energy_rel = 1e-8       # relative")
    lines.append(f"position_abs = 1e-9     # nm")
    lines.append(f"")
    lines.append(f"[reference_values]")
    lines.append(f"# Key values from step 0 for quick validation")
    lines.append(f"# (parse the .tre/.trf files for full precision)")
    for key, val in key_values.items():
        # Sanitize key for TOML (replace special chars)
        toml_key = (
            key.replace("(", "_").replace(")", "").replace(" ", "_").replace("-", "_")
        )
        lines.append(f"{toml_key} = {val}")
    lines.append(f"")

    # Add first few forces as quick-check values
    if step0_forces:
        lines.append(f"[reference_values.forces_step0]")
        lines.append(f"# Force on atom 0 at step 0 [fx, fy, fz] in kJ/(mol*nm)")
        lines.append(
            f"atom_0 = [{step0_forces[0][0]}, {step0_forces[0][1]}, {step0_forces[0][2]}]"
        )
        if len(step0_forces) > 1:
            lines.append(
                f"atom_1 = [{step0_forces[1][0]}, {step0_forces[1][1]}, {step0_forces[1][2]}]"
            )
        if len(step0_forces) > 2:
            lines.append(
                f"atom_2 = [{step0_forces[2][0]}, {step0_forces[2][1]}, {step0_forces[2][2]}]"
            )

    return "\n".join(lines) + "\n"


def process_system(md_binary, system, verbose=False):
    """Run one system and save raw outputs + TOML manifest."""
    sys_dir = REFERENCES_DIR / system["dir"]
    expected_dir = sys_dir / "expected"
    expected_dir.mkdir(exist_ok=True)

    # Create work directory
    work_dir = Path(tempfile.mkdtemp(prefix=f"gromos_ref_{system['name']}_"))

    try:
        run_output = run_gromosxx(md_binary, system, work_dir, verbose=verbose)

        if run_output is None:
            # Even on failure, save the log to expected/ for inspection
            log_src = work_dir / "md_output.log"
            if log_src.is_file():
                shutil.copy2(log_src, expected_dir / "md_output.log")
                print(f"  Log saved to: {expected_dir / 'md_output.log'}")
            return False

        # Copy raw output files to expected/
        for src_name, dst_name in [
            ("trf", "forces.trf"),
            ("tre", "energies.tre"),
            ("trc", "trajectory.trc"),
            ("fin", "final.conf"),
            ("log", "md_output.log"),
        ]:
            src = run_output[src_name]
            dst = expected_dir / dst_name
            if src.is_file():
                shutil.copy2(src, dst)

        # Extract key values for TOML
        key_values = extract_key_values(run_output)
        step0_forces = extract_step0_forces(expected_dir / "forces.trf")

        # Generate TOML manifest
        toml_content = generate_toml(system, run_output, key_values, step0_forces)
        toml_path = sys_dir / "input.toml"
        with open(toml_path, "w") as f:
            f.write(toml_content)

        # Report sizes
        total_size = sum(
            (expected_dir / f).stat().st_size
            for f in ["forces.trf", "energies.tre", "trajectory.trc", "final.conf"]
            if (expected_dir / f).is_file()
        )
        print(f"  Saved to: {expected_dir}/")
        print(f"  Total size: {total_size / 1024:.1f} KB")
        print(f"  Manifest: {toml_path}")

        return True

    except Exception as e:
        print(f"  EXCEPTION: {e}")
        import traceback

        traceback.print_exc()
        return False

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description="Generate reference data for gromos-rs validation"
    )
    parser.add_argument(
        "--md-binary",
        "-m",
        help="Path to compiled gromosXX md binary",
        default=None,
    )
    parser.add_argument(
        "--systems",
        "-s",
        nargs="+",
        help="Run only specific systems (by name)",
        default=None,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Show full gromosXX output (especially on failure)",
    )

    args = parser.parse_args()

    md_binary = find_md_binary(args.md_binary)
    if not md_binary:
        print("ERROR: Cannot find md binary!")
        print("Build gromosXX first, then specify path with --md-binary")
        print("")
        print("Build with:")
        print("  cd md++ && mkdir build && cd build")
        print("  cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)")
        sys.exit(1)

    print(f"Using md binary: {md_binary}")

    # Filter systems
    systems_to_run = SYSTEMS
    if args.systems:
        systems_to_run = [s for s in SYSTEMS if s["name"] in args.systems]
        if not systems_to_run:
            print(
                f"ERROR: No matching systems. Available: {[s['name'] for s in SYSTEMS]}"
            )
            sys.exit(1)

    print(f"\nRunning {len(systems_to_run)} reference systems:\n")

    results = {}
    for system in systems_to_run:
        print(f"{'─' * 60}")
        print(f" {system['name']} (Level {system['level']}): {system['isolates']}")
        print(f"{'─' * 60}")

        ok = process_system(md_binary, system, verbose=args.verbose)
        results[system["name"]] = ok
        print()

    # Summary
    print(f"{'═' * 60}")
    print("SUMMARY")
    print(f"{'═' * 60}")
    for name, ok in results.items():
        icon = "✓" if ok else "✗"
        print(f"  {icon} {name}")

    failed = sum(1 for ok in results.values() if not ok)
    if failed == 0:
        print(f"\nAll {len(results)} systems passed.")
        print(
            "Copy each system directory to gromos-rs/tests/data/ for use in Rust tests."
        )
    else:
        print(f"\n{failed} system(s) failed.")

    return failed


if __name__ == "__main__":
    sys.exit(main())
