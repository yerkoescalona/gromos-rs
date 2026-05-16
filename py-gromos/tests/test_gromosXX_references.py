"""
Integration tests that exercise the GROMOS-RS engine through the Python
Simulation API (gromos.Simulation) and compare output against gromosXX
reference data.

This validates the full Rust simulation pipeline is accessible from Python
without shelling out to the md binary.

Usage style mirrors OpenMM but with gromosXX naming conventions:

    sim = Simulation("system.topo", "initial.cnf", "run.imd")
    sim.step(100)
    print(sim.energies, sim.positions, sim.forces)
"""

import re
from pathlib import Path

import numpy as np
import pytest

from gromos import Configuration, InputParameters, Simulation, Topology

# Paths
REPO_ROOT = Path(__file__).parent.parent.parent
REF_DIR = REPO_ROOT / "crates" / "gromos-cli" / "tests" / "gromosXX_references"

# Tolerances (same as Rust tests)
ENERGY_REL_TOL = 1e-8
FORCE_ABS_TOL = 1e-6
POSITION_ABS_TOL = 1e-9


# ============================================================================
# Helpers to parse reference (expected) data files
# ============================================================================


def _parse_input_toml(system_dir):
    """Parse input.toml to get file paths for this reference system."""
    toml_path = system_dir / "input.toml"
    text = toml_path.read_text()

    def get_val(key):
        m = re.search(rf'^{key}\s*=\s*"(.+)"', text, re.MULTILINE)
        return m.group(1) if m else None

    return {
        "topology": get_val("topology"),
        "configuration": get_val("configuration"),
        "parameters": get_val("parameters"),
    }


def _parse_expected_energies(tre_path):
    """Parse gromosXX ENERGY03 block format reference energies.

    Returns list of dicts: {total, kinetic, potential}
    """
    text = tre_path.read_text()
    frames = []
    blocks = text.split("TIMESTEP")
    for block in blocks[1:]:
        lines = block.strip().split("\n")
        energy_start = None
        for i, line in enumerate(lines):
            if "ENERGY03" in line:
                energy_start = i + 1
                break
        if energy_start is None:
            continue

        vals = []
        for i in range(energy_start, len(lines)):
            line = lines[i].strip()
            if line.startswith("#") and "totals" in line:
                continue
            if line.startswith("#") or line == "END":
                break
            try:
                vals.append(float(line))
            except ValueError:
                break

        if len(vals) >= 3:
            frames.append(
                {
                    "total": vals[0],
                    "kinetic": vals[1],
                    "potential": vals[2],
                }
            )
    return frames


def _parse_expected_positions(trc_path):
    """Parse gromosXX POSITIONRED trajectory (3 columns: x, y, z).

    Returns list of Nx3 numpy arrays, one per frame.
    """
    text = trc_path.read_text()
    frames = []
    blocks = text.split("POSITIONRED")
    for block in blocks[1:]:
        lines = block.strip().split("\n")
        positions = []
        for line in lines:
            line = line.strip()
            if line == "END" or line.startswith("TIMESTEP") or line.startswith("GENBOX"):
                break
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                try:
                    positions.append([float(x) for x in parts])
                except ValueError:
                    break
        if positions:
            frames.append(np.array(positions, dtype=np.float64))
    return frames


def _parse_expected_forces(trf_path):
    """Parse gromosXX FREEFORCERED trajectory (3 columns per atom).

    Returns list of Nx3 numpy arrays, one per frame.
    """
    text = trf_path.read_text()
    frames = []
    blocks = text.split("FREEFORCERED")
    for block in blocks[1:]:
        lines = block.strip().split("\n")
        forces = []
        for line in lines:
            line = line.strip()
            if line == "END" or line.startswith("CONSFORCERED"):
                break
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 3:
                try:
                    forces.append([float(x) for x in parts])
                except ValueError:
                    break
        if forces:
            frames.append(np.array(forces, dtype=np.float64))
    return frames


def _create_simulation(system_dir):
    """Create a Simulation from a reference system directory using compositional API."""
    inputs = _parse_input_toml(system_dir)
    topo = Topology(str((system_dir / inputs["topology"]).resolve()))
    conf = Configuration(str((system_dir / inputs["configuration"]).resolve()))
    params = InputParameters(str((system_dir / inputs["parameters"]).resolve()))
    return Simulation(topo, conf, params)


def _get_n_steps(system_dir):
    """Read NSTLIM from InputParameters."""
    inputs = _parse_input_toml(system_dir)
    params = InputParameters(str((system_dir / inputs["parameters"]).resolve()))
    return params.nstlim


def _get_write_frequencies(system_dir):
    """Read write frequencies from InputParameters.

    Returns (ntwe, ntwx, ntwf) — energy, trajectory, and force write intervals.
    For ntwf, falls back to parsing the WRITETRAJ block since InputParameters
    doesn't expose it directly yet.
    """
    inputs = _parse_input_toml(system_dir)
    params = InputParameters(str((system_dir / inputs["parameters"]).resolve()))
    ntwe = params.ntwe if params.ntwe > 0 else 1
    ntwx = params.ntwx if params.ntwx > 0 else 1

    # ntwf not yet exposed on InputParameters — parse from raw file
    param_file = system_dir / inputs["parameters"]
    text = param_file.read_text()
    m = re.search(
        r"WRITETRAJ\n#[^\n]*\n\s*(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)",
        text,
    )
    ntwf = int(m.group(2)) if m else ntwx
    if ntwf == 0:
        ntwf = 1
    return ntwe, ntwx, ntwf


# ============================================================================
# Reference test systems — same as test_gromosXX_references.rs
# ============================================================================

REFERENCE_SYSTEMS = [
    # Level 0: Pair interactions
    "pair_lj",
    "pair_lj_mixed",
    "nacl_pair",
    # Level 1: Single molecule
    "water_single",
    "benzene_vacuum",
    "nacl_pair_box",
    "butane_vacuum",
    "aladip_vacuum",
    # Level 2: Small solvated
    "water_3_box",
    "nacl_1water_box",
    "nacl_3water_box",
    "water_3_box_twinrange",
    "water_10_box",
    "nacl_3water_cutoff",
    "nacl_water_box",
    "nacl_water_box_shifted",
    # Level 3: Bulk
    "water_216_box",
    "water_216_box_com",
    "water_216_nvt",
    "water_216_npt",
    # Level 4: Full system
    "aladip_solvated",
]


# ============================================================================
# Energy tests — through Python Simulation API
# ============================================================================


@pytest.mark.parametrize("system_name", REFERENCE_SYSTEMS)
def test_reference_energies(system_name):
    """Compare energies from gromos.Simulation against gromosXX reference."""
    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")

    expected_tre = system_dir / "expected" / "energies.tre"
    if not expected_tre.exists():
        pytest.skip(f"No expected energies for {system_name}")

    expected = _parse_expected_energies(expected_tre)
    n_expected = len(expected)
    assert n_expected > 0, f"No expected energy frames for {system_name}"

    n_steps = _get_n_steps(system_dir)
    ntwe, _, _ = _get_write_frequencies(system_dir)

    # Create simulation (constructor initializes + runs step 0)
    sim = _create_simulation(system_dir)

    # Collect energies at write intervals (matching reference output)
    actual = [
        {
            "total": sim.total_energy,
            "kinetic": sim.kinetic_energy,
            "potential": sim.potential_energy,
        }
    ]

    for step in range(1, n_steps + 1):
        sim.step(1)
        if step % ntwe == 0:
            actual.append(
                {
                    "total": sim.total_energy,
                    "kinetic": sim.kinetic_energy,
                    "potential": sim.potential_energy,
                }
            )

    # Trim to match reference frame count
    if len(actual) > n_expected:
        actual = actual[:n_expected]

    assert len(actual) == n_expected, (
        f"{system_name}: frame count mismatch: {len(actual)} vs {n_expected}"
    )

    for i, (act, exp) in enumerate(zip(actual, expected)):
        for key in ["total", "kinetic", "potential"]:
            a, e = act[key], exp[key]
            if abs(e) > 1e-10:
                rel = abs(a - e) / abs(e)
                assert rel < ENERGY_REL_TOL, (
                    f"{system_name} frame {i} {key}: {a} vs {e} (rel={rel:.2e})"
                )
            else:
                assert abs(a - e) < 1e-10, f"{system_name} frame {i} {key}: {a} vs {e}"


# ============================================================================
# Position tests — through Python Simulation API
# ============================================================================


@pytest.mark.parametrize("system_name", REFERENCE_SYSTEMS)
def test_reference_positions(system_name):
    """Compare positions from gromos.Simulation against gromosXX reference."""
    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")

    expected_trc = system_dir / "expected" / "trajectory.trc"
    if not expected_trc.exists():
        pytest.skip(f"No expected trajectory for {system_name}")

    expected = _parse_expected_positions(expected_trc)
    n_expected = len(expected)
    assert n_expected > 0, f"No expected position frames for {system_name}"

    n_steps = _get_n_steps(system_dir)
    _, ntwx, _ = _get_write_frequencies(system_dir)

    sim = _create_simulation(system_dir)

    # Get box dimensions for minimum image convention
    inputs = _parse_input_toml(system_dir)
    conf = Configuration(str((system_dir / inputs["configuration"]).resolve()))
    box_dims = np.array(conf.box_dimensions)

    # Collect positions at write intervals
    actual = [sim.positions.copy()]
    for step in range(1, n_steps + 1):
        sim.step(1)
        if step % ntwx == 0:
            actual.append(sim.positions.copy())

    if len(actual) > n_expected:
        actual = actual[:n_expected]

    assert len(actual) == n_expected, (
        f"{system_name}: position frame count mismatch: {len(actual)} vs {n_expected}"
    )

    for i, (act, exp) in enumerate(zip(actual, expected)):
        diff = act - exp
        # Apply minimum image convention for periodic systems
        if box_dims.min() > 0:
            diff -= np.round(diff / box_dims) * box_dims
        max_diff = np.abs(diff).max()
        assert max_diff < POSITION_ABS_TOL, (
            f"{system_name} frame {i}: max position diff = {max_diff:.2e}"
        )


# ============================================================================
# Force tests — through Python Simulation API
# ============================================================================


@pytest.mark.parametrize("system_name", REFERENCE_SYSTEMS)
def test_reference_forces(system_name):
    """Compare forces from gromos.Simulation against gromosXX reference."""
    system_dir = REF_DIR / system_name
    if not system_dir.exists():
        pytest.skip(f"Reference system {system_name} not found")

    expected_trf = system_dir / "expected" / "forces.trf"
    if not expected_trf.exists():
        pytest.skip(f"No expected forces for {system_name}")

    expected = _parse_expected_forces(expected_trf)
    n_expected = len(expected)
    assert n_expected > 0, f"No expected force frames for {system_name}"

    n_steps = _get_n_steps(system_dir)
    _, _, ntwf = _get_write_frequencies(system_dir)

    sim = _create_simulation(system_dir)

    # Collect forces at write intervals
    actual = [sim.forces.copy()]
    for step in range(1, n_steps + 1):
        sim.step(1)
        if step % ntwf == 0:
            actual.append(sim.forces.copy())

    if len(actual) > n_expected:
        actual = actual[:n_expected]

    assert len(actual) == n_expected, (
        f"{system_name}: force frame count mismatch: {len(actual)} vs {n_expected}"
    )

    for i, (act, exp) in enumerate(zip(actual, expected)):
        diff = np.abs(act - exp)
        max_diff = diff.max()
        assert max_diff < FORCE_ABS_TOL, f"{system_name} frame {i}: max force diff = {max_diff:.2e}"
