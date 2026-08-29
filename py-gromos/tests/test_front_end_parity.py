"""
Front-end parity — PLAN.md 3.9, step 0 (drift guard G2, v0).

The reference suite (`test_gromosXX_references.py`) guards each Python construction path
against gromosXX individually; nothing guarded the paths against *each other* — the trap
FUTURE.md named when it parked the builder unification. This file is that missing test.

v0 compares the two builders that exist today:

  path A — ``Simulation(topo, conf, params, **restraint_kwargs)``   (``build_simulation``)
  path C — ``Simulation.from_sequence(topo, conf, params,
                AlgorithmSequence.from_parameters(topo, params))``   (``resolve_algorithm_sequence``)

on every reference system, for NSTLIM steps, and asserts **exact** equality
(``np.array_equal``, never ``allclose``) of the per-step energy rows, final positions and final
forces. Same build, same kernels, same order: any difference is a defect or a missing feature,
never a tolerance question (PLAN.md 3.9 A1/A3).

Known divergences are ``xfail(strict=True)`` with a reason naming the feature, so that when the
shared builder lands (steps 1–2) they *must* be removed — a silent pass is drift in the other
direction. Rows are mirrored in PLAN.md 3.9's divergence table.

``test_same_path_twice`` is the determinism baseline (A3): path A run twice must be identical
before any cross-path comparison means anything. The Python path is serial
(``parallel_nonbonded`` is never set by ``build_simulation``); the multi-thread outcome of the
*binary* is measured separately by ``scripts/kernel_determinism.py`` and recorded in
BENCHMARKING.md.
"""

from __future__ import annotations

import numpy as np
import pytest

from gromos import AlgorithmSequence, Configuration, InputParameters, Simulation, Topology

from test_gromosXX_references import (
    REF_DIR,
    REFERENCE_SYSTEMS,
    _get_n_steps,
    _parse_input_toml,
)

# ---------------------------------------------------------------------------------------------
# Known divergences between path A and path C (filled from the first run of this file, step 0).
# Each entry names the *feature* path C lacks or gets wrong. strict=True: when the shared
# builder makes the two paths agree, the entry must be deleted, not left to rot.
# ---------------------------------------------------------------------------------------------
EXPECTED_DIVERGENCE: dict[str, str] = {
    # Constraint algorithms: the descriptor enum has only ShakeConstraints; path C runs SHAKE
    # where the IMD asks for SETTLE (NTCS=3) or LINCS (NTCP/NTCS=2).
    "nacl_1water_settle": "path C lacks SETTLE (ShakeConstraints is the only constraint descriptor)",
    "nacl_1water_lincs": "path C lacks LINCS for solvent (NTCS=2 resolved as SHAKE)",
    "aladip_vacuum_lincs": "path C lacks LINCS for solute (NTCP=2 resolved as SHAKE)",
    # Thermostats: BerendsenThermostat is the only thermostat descriptor.
    "water_216_nvt_nosehoover": "path C lacks Nosé-Hoover (MULTIBATH algorithm 1 resolved as Berendsen)",
    "water_216_nvt_nhc_chain": "path C lacks Nosé-Hoover chains (MULTIBATH algorithm >=2 resolved as Berendsen)",
    # COM motion: RemoveCOMMotion(initial: bool, nscm) collapses NTICOM=2 (translation +
    # rotation) to NTICOM=1 — the rotation removal is lost on path C.
    "water_216_box_com_rot": "path C collapses NTICOM=2 to a bool (rotation removal lost)",
    # Boundary: build_simulation_from_sequence hard-codes SimBox::rectangular; the NTB=-1
    # truncated-octahedron transform never happens on path C.
    "aladip_trunc_oct": "path C lacks the NTB=-1 truncated-octahedron box transform",
    # Restraints: constructor kwargs on path A, inexpressible on path C.
    "nacl_1water_distres": "restraints are Simulation kwargs; Simulation.from_sequence has none (distrest)",
    "aladip_solvated_em_posres": "restraints are Simulation kwargs; Simulation.from_sequence has none (posresspec/refpos)",
    "aladip_solvated_em": "restraints are Simulation kwargs; Simulation.from_sequence has none (posresspec/refpos)",
}


def _load(system_name: str):
    system_dir = REF_DIR / system_name
    inputs = _parse_input_toml(system_dir)
    topo = Topology(str((system_dir / inputs["topology"]).resolve()))
    conf = Configuration(str((system_dir / inputs["configuration"]).resolve()))
    params = InputParameters(str((system_dir / inputs["parameters"]).resolve()))
    kwargs = {}
    for key in ("distrest", "posresspec", "refpos"):
        if inputs.get(key):
            kwargs[key] = str((system_dir / inputs[key]).resolve())
    return topo, conf, params, kwargs


def _trace(sim: Simulation, n_steps: int):
    """Per-step energy rows (frame 0 = state after step 0), final positions and forces."""
    frames = sim.run(n_steps, ene_freq=1)
    return frames, sim.positions, sim.forces


def _run_path_a(system_name: str):
    topo, conf, params, kwargs = _load(system_name)
    sim = Simulation(topo, conf, params, **kwargs)
    return _trace(sim, _get_n_steps(REF_DIR / system_name))


def _run_path_c(system_name: str):
    topo, conf, params, kwargs = _load(system_name)
    if kwargs:
        pytest.fail(
            f"{system_name}: restraint inputs {sorted(kwargs)} cannot be expressed on the "
            "AlgorithmSequence path (Simulation.from_sequence has no restraint arguments)"
        )
    seq = AlgorithmSequence.from_parameters(topo, params)
    sim = Simulation.from_sequence(topo, conf, params, seq)
    return _trace(sim, _get_n_steps(REF_DIR / system_name))


def _first_difference(a: np.ndarray, b: np.ndarray) -> str:
    """Human-readable location of the first exact mismatch (for the divergence table)."""
    if a.shape != b.shape:
        return f"shape {a.shape} vs {b.shape}"
    idx = np.argwhere(a != b)
    if idx.size == 0:
        return "identical"
    first = tuple(int(i) for i in idx[0])
    return f"first mismatch at index {first}: {a[first]!r} vs {b[first]!r}"


def _assert_identical(label: str, a, b) -> None:
    fa, pa, xa = a
    fb, pb, xb = b
    problems = []
    if not np.array_equal(fa, fb):
        problems.append(f"energy rows — {_first_difference(fa, fb)}")
    if not np.array_equal(pa, pb):
        problems.append(f"final positions — {_first_difference(pa, pb)}")
    if not np.array_equal(xa, xb):
        problems.append(f"final forces — {_first_difference(xa, xb)}")
    assert not problems, f"{label}:\n  " + "\n  ".join(problems)


def _params_for(system_name: str):
    marks = []
    if system_name in EXPECTED_DIVERGENCE:
        marks.append(pytest.mark.xfail(strict=True, reason=EXPECTED_DIVERGENCE[system_name]))
    return pytest.param(system_name, marks=marks, id=system_name)


@pytest.mark.parametrize("system_name", [pytest.param(s, id=s) for s in REFERENCE_SYSTEMS])
def test_same_path_twice(system_name):
    """Determinism baseline (A3): the same construction path twice must be bit-identical."""
    _assert_identical(f"{system_name}: path A vs path A", _run_path_a(system_name), _run_path_a(system_name))


@pytest.mark.parametrize("system_name", [_params_for(s) for s in REFERENCE_SYSTEMS])
def test_front_end_parity(system_name):
    """Path A (build_simulation) vs path C (from_parameters + resolve_algorithm_sequence)."""
    _assert_identical(
        f"{system_name}: path A (Simulation) vs path C (AlgorithmSequence.from_parameters)",
        _run_path_a(system_name),
        _run_path_c(system_name),
    )
