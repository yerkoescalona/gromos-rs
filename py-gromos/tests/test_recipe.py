"""
The run recipe behind a Simulation — PLAN.md 3.9 step 2.

A `Simulation` built from parameters exposes the effective recipe (TOML), the algorithm plan
(JSON, one fully resolved entry per algorithm) and the parser diagnostics. Step 3 turns the
recipe into a first-class Python object; this file pins down what step 2 already guarantees.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from gromos import AlgorithmSequence, Configuration, InputParameters, Simulation, System, Topology

from test_gromosXX_references import REF_DIR, _parse_input_toml

WATER = REF_DIR / "water_216_box"
NOBATH = REF_DIR / "water_216_nve_nobath"


def _load(system_dir: Path):
    inputs = _parse_input_toml(system_dir)
    topo = Topology(str((system_dir / inputs["topology"]).resolve()))
    conf = Configuration(str((system_dir / inputs["configuration"]).resolve()))
    params = InputParameters(str((system_dir / inputs["parameters"]).resolve()))
    return topo, conf, params


def test_recipe_and_plan_are_exposed_and_ordered():
    topo, conf, params = _load(WATER)
    sim = Simulation(topo, conf, params)

    recipe = sim.recipe_toml
    assert recipe is not None
    assert "version = 1" in recipe
    for group in ("[control]", "[boundary]", "[forcefield", "[constraints]", "[outputs]"):
        assert group in recipe, f"missing {group} in\n{recipe}"

    plan = json.loads(sim.plan_json)
    kinds = [entry["kind"] for entry in plan]
    # GROMOS step order, one entry per algorithm, and it matches what actually runs.
    assert kinds[0] == "forcefield" or kinds[:2] == ["remove_com", "forcefield"]
    assert kinds[-1] == "energy_calculation"
    assert kinds.index("leap_frog_velocity") < kinds.index("leap_frog_position")
    assert len(kinds) == len(sim.algorithm_names)
    forcefield = next(e for e in plan if e["kind"] == "forcefield")
    # Fully resolved values, not "look it up in the parameter file".
    assert forcefield["cutoff_long"] == pytest.approx(params.cutoff)
    assert forcefield["four_pi_eps_i"] == pytest.approx(138.9354)
    assert forcefield["parallel"] is True  # 648 atoms > 100: the binary's policy


def test_diagnostics_report_absent_blocks():
    topo, conf, params = _load(NOBATH)
    sim = Simulation(topo, conf, params)
    notes = sim.diagnostics
    assert any("MULTIBATH block absent" in n and "no temperature coupling" in n for n in notes), notes
    # …and the plan indeed has no thermostat.
    kinds = [e["kind"] for e in json.loads(sim.plan_json)]
    assert "thermostat" not in kinds

    topo, conf, params = _load(WATER)
    sim = Simulation(topo, conf, params)
    # A reference input carries every block it needs; the notes only name the optional blocks
    # it deliberately omits (no pressure coupling, no restraints, no perturbation, no EM).
    assert not any("MULTIBATH" in n or "COMTRANSROT" in n for n in sim.diagnostics), sim.diagnostics
    assert all("block absent" in n for n in sim.diagnostics), sim.diagnostics


def test_descriptor_path_has_no_recipe():
    topo, conf, params = _load(WATER)
    seq = AlgorithmSequence.from_parameters(topo, params)
    sim = Simulation.from_sequence(topo, conf, params, seq)
    assert sim.recipe_toml is None
    assert sim.plan_json is None
    assert sim.diagnostics == []


def test_factory_parameters_default_to_what_gromosxx_does_without_the_blocks():
    """No COMTRANSROT in a factory-built run means no COM-motion removal (gromosXX:
    skip_step 0), and NVE means no bath — the recipe shows both."""
    system = System.from_files(
        str((WATER / "water_216_box.topo").resolve()), str((WATER / "water_216_box.conf").resolve())
    )
    sim = Simulation(system, InputParameters.nve(dt=0.002, steps=2))
    kinds = [e["kind"] for e in json.loads(sim.plan_json)]
    assert "remove_com" not in kinds
    assert "thermostat" not in kinds
