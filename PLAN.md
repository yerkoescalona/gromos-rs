# gromos-rs — Status & Plan

Focus: `cargo build --release --bin md`
On commit: update CHANGELOG.md and Cargo.toml version.
DON'T modify `gromosXX_references/*/expected/` — those are ground truth.

References:
- gromosXX source (MD engine, "md++"): `.local/gromosXX/md++/src`
- gromosPlsPls source (analysis/tools, "gromos++"): `.local/gromosPlsPls/gromos++/src`
- Tutorials: `.local/gromos_tutorial_livecoms/tutorial_files`
- Theory: `.local/doc/gromos_book`
- Force fields: `.local/gromosXX/forcefields`

Doc style: Rust → KaTeX + `[^label]` footnotes; Python → NumPy docstrings + `.. math::`


## Architecture

| Crate | Role |
|-------|------|
| gromos-core | Types, math (Vec3), BoundaryCondition, Topology, Configuration, Algorithm trait, AtomSelection |
| gromos-forces | Bonded, nonbonded (LJ+CRF), pairlist, restraints, PME/QM-MM (stubs) |
| gromos-integrators | LeapFrog, SteepestDescent, SHAKE/SETTLE/LINCS, thermostats, barostats, EDS/GaMD/REMD/FEP |
| gromos-io | All file I/O: topology, coordinates, IMD, trajectories, energy, MTB/IFP, PDB, posres |
| gromos-md | Simulation engines: md, md_mpi, remd, eds, gamd, mdf |
| gromos-tools | System construction: make_top, sim_box, pdb2g96, mk_script, ion, prep_posres, ... |
| gromos-analysis | Trajectory analysis: rmsd, rdf, ene_ana, hbond, ... |
| pyo3-gromos | Python bindings (PyO3) |
| py-gromos | Python package (separate maturin build) |

**Design rules:**
- All file format parsing/writing lives exclusively in `gromos-io`. No duplication across crates.
- **Seamless gromosXX ↔ gromosPlsPls merge, no code duplication.** Shared primitives (math,
  boundary/PBC, topology, atom selection, rotational fit, statistics, single-point energy) live in
  the lower crates and are consumed by *both* the MD engine (gromosXX = forces/integrators/md) and
  the analysis/tools suite (gromosPlsPls = analysis/tools). The gromos++-style API is a **facade**
  over the gromosXX primitives: internally it uses the engine libraries; to the user it behaves like
  gromosPlsPls. In practice this means `gromos-analysis` should reuse engine code (energy eval,
  geometry, selection) rather than reimplement it — today it only depends on `gromos-core` + `gromos-io`.

## Decisions Taken

- f64 everywhere (not f32)
- gromosXX `@` CLI convention: `@topo @conf @input @fin @trc @tre @trf @trv @verb ...`
- All simulation parameters from `@input` .imd/.in file (not CLI flags)
- Bonded force vectors: `v = pos(i) - pos(j)` (gromosXX convention)
- Boundary condition from box_dims: vacuum if (0,0,0), rectangular otherwise
- GENBOX block parsed (box_type + dims) but box_type ignored — only dims used
- Energy output: full f64 scientific notation for exact comparison
- Tolerances: force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9
- CLI arg parsing: use clap `#[derive(Parser)]` with a `gromos_args()` pre-processor that translates `@key` → `--key` and expands `@f argfile`. No custom arg parsers.
- File format parsers (MTB, IFP, topology, coordinates, etc.) live exclusively in `gromos-io` — no duplication in tool binaries.
- Crate restructuring: gromos-md (simulation engines), gromos-tools (system construction), gromos-analysis (trajectory analysis)
- **AtomSpecifier strategy:** solve the gromosXX-level atom/topology primitives first (queryable
  per-atom metadata, including solvent), then build the gromos++ `AtomSpecifier` grammar as a facade
  on top of those primitives — zero duplication, facade behaves like gromosPlsPls.
- **Reference tests for analysis/tools:** hand-craft *minimal* systems inspired by the tutorials and
  compare against the C++ reference, mirroring the MD `gromosXX_references/` harness. (Full tutorial
  end-to-end runs are compute-limited and deferred to last.)

## Reference Test Status

Validation: `cargo test -p gromos-md --test test_gromosXX_references`
Generate refs: `python3 crates/gromos-md/tests/run_references.py --md-binary .local/gromosXX/md++/build/program/md`
Test file: `crates/gromos-md/tests/test_gromosXX_references.rs`
Ref data: `crates/gromos-md/tests/gromosXX_references/`

| Lvl | System           | Atoms | Isolates                              | Status   |
|-----|------------------|-------|---------------------------------------|----------|
| 0   | pair_lj          | 2     | Pure LJ, no PBC                      | **PASS** |
| 0   | pair_lj_mixed    | 2     | LJ combination rules                 | **PASS** |
| 0   | nacl_pair        | 2     | Coulomb + LJ ions                    | **PASS** |
| 1   | water_single     | 3     | bond + angle + intramolecular CRF    | **PASS** |
| 1   | water_single_genvel | 3  | NTIVEL=1 Maxwell-Boltzmann velocity generation | **PASS** |
| 1   | benzene_vacuum   | 12    | aromatic ring + improper + torsion   | **PASS** |
| 1   | nacl_pair_box    | 2     | Coulomb + LJ in PBC with RF (no solvent) | **PASS** |
| 1   | butane_vacuum    | 4     | dihedral + 1-4 LJ interaction        | **PASS** |
| 1   | aladip_vacuum    | 12    | all bonded + exclusions + 1-4        | **PASS** |
| 2   | water_3_box      | 9     | PBC + min image + pairlist + CRF     | **PASS** |
| 2   | nacl_1water_box  | 5     | minimal solute-solvent + SHAKE       | **PASS** |
| 2   | nacl_3water_box  | 11    | multiple solvent + solute-solvent pairlist | **PASS** |
| 2   | water_3_box_twinrange | 9 | twin-range pairlist (RCUTP<RCUTL, NSNB=5) | **PASS** |
| 2   | water_10_box     | 32    | 2 ions + 10 SPC, positions away from cutoff | **PASS** |
| 2   | nacl_3water_cutoff | 11  | nacl_3water near cutoff boundary     | **PASS** |
| 2   | nacl_water_box   | 62    | ion-water RF in PBC                  | **PASS** |
| 2   | nacl_water_box_shifted | 62 | nacl_water_box with perturbed positions | **PASS** |
| 3   | water_216_box    | 648   | bulk NVE, pairlist, virial           | **PASS** |
| 3   | water_216_box_com| 648   | bulk NVE + COM removal (NTICOM=1, NSCM=10) | **PASS** |
| 3   | water_216_nvt    | 648   | Berendsen thermostat                 | **PASS** |
| 3   | water_216_npt    | 648   | Berendsen barostat                   | **PASS** |
| 4   | aladip_solvated  | 72    | SHAKE + solute-solvent               | **PASS** |
| 4   | aladip_vacuum_em | 12    | steepest descent EM, vacuum          | **PASS** |
| 4   | aladip_vacuum_em_shake | 12 | SD EM + SHAKE, vacuum             | **PASS** |
| 4   | aladip_solvated_em_noshake | 72 | SD EM, solvated, no SHAKE      | **PASS** |
| 4   | aladip_solvated_em_shake | 72 | SD EM + SHAKE, solvated          | **PASS** |
| 4   | aladip_solvated_em_posres | 72 | SD EM + position restraints     | **PASS** |
| 4   | aladip_solvated_em | 72  | SD EM + SHAKE + posres, solvated    | **PASS** |

**28 of 28 tests pass.** All levels fully passing.

(No reference tests yet for `gromos-analysis` / `gromos-tools` — see Roadmap Priority 2 + the
cross-cutting minimal-reference-test theme.)

## What Works

### gromos-core
- Boundary conditions: Vacuum, Rectangular (minimum image), Triclinic (defined but not wired)
- Topology struct with solvent expansion: NSM auto-computed, chargegroups from CGC codes, intra-molecular exclusions
- AtomSelection: basic grammar only (atom/residue/molecule-1/all-solvent) — **underbuilt**, see Roadmap P2

### gromos-forces
- LJ + CRF nonbonded (vacuum and PBC), 1-4 interactions (cs6/cs12 + scaled CRF)
- All bonded types: quartic/harmonic bonds, cos-harmonic/harmonic angles, dihedrals, impropers, cross-dihedrals
- NTF flag control, RF excluded interactions (forces + energy + self-terms)
- Pairlist: chargegroup-based, atom-based, twin-range (RCUTP/RCUTL with force caching)
- Position restraints: harmonic F = -k·(r - r_ref), from .por + .rpr files

### gromos-integrators
- SHAKE: solute (NTC>1) + solvent (NTCS>0), virial contribution, skip optimization, NTISHK
- COM motion removal: NTICOM (initial) + NSCM (periodic), translational only
- SETTLE, LINCS (implemented but not wired)
- Berendsen thermostat & barostat (MULTIBATH/PRESSURESCALE)
- Pressure/virial: prepare_virial + atomic_to_molecular_virial, PressureCalculation, BerendsenBarostat
- Steepest Descent EM: adaptive step, FLIM, DELE, NMIN, works with SHAKE + posres

### gromos-io
- Topology parser: all blocks (SOLUTEATOM, SOLVENTCONSTR, LJPARAMETERS, PRESSUREGROUPS, ...)
- Coordinate reader: POSITION, POSITIONRED, VELOCITY, VELOCITYRED, GENBOX
- Energy/trajectory/force writers: ENERTRJ, POSITIONRED, FREEFORCERED/CONSFORCERED
- MTB parser: MTBUILDBLSOLUTE, MTBUILDBLEND, MTBUILDBLSOLVENT with continuation lines
- IFP parser: all type codes + SINGLEATOMLJPAIR (selection matrix) + MIXEDATOMLJPAIR
- Position restraints parser (.por/.rpr)

### gromos-md
- Force trajectory output (@trf): FREEFORCERED + CONSFORCERED, atom-by-atom tolerance 1e-6

### gromos-tools
- make_top: MTB + IFP → topology (end groups, exclusions, 1-4 pairs, LJ matrix, combining rules)
  - Tested: GB3 (56 res, 457 atoms) + Na+ with 54A7
- System building done: make_top, com_top, check_top, pdb2g96, sim_box, ion, mk_script, make_pt_top, prep_posres, build_box

---

## Roadmap (priority order)

Overarching principle: **gromosXX as the reference, no duplication between the engine (gromosXX) and
the analysis/tools facade (gromosPlsPls).** Every feature lands with a minimal gromosXX reference test.

### Priority 1 — MD engine physics (gromosXX-faithful, reference-tested)
Wire the already-coded-but-unwired physics; keep implementations in `gromos-forces`/`gromos-integrators`
(reusable), never duplicated into binaries. Each item gets a minimal reference test.

**1.1 — Reproducibility & correctness (small, do first)**
- [x] **NTIVEL=1 velocity generation** (Maxwell-Boltzmann) — `gromos-core/src/random.rs`
  (`GslMt19937` + `gsl_ran_gaussian` + `generate_velocities`), wired in `md.rs` velocity setup.
  Bit-for-bit match verified via `water_single_genvel` reference test (kinetic energy at step 0
  depends entirely on generated velocities and matches gromosXX to full precision).
  - [x] Read NTIVEL, IG (seed), TEMPI from INITIALISE block (already parsed in `imd.rs`)
  - [x] MT19937 RNG matching GSL's `gsl_rng_mt19937`; Gaussian matching `gsl_ran_gaussian` (polar Box-Muller)
  - [x] Per atom: σ = sqrt(k_B·T / m_i), v_i = gaussian(σ) for x,y,z; k_B = 0.00831441 kJ/(mol·K)
  - [x] Store in both current().vel and old().vel (gromosXX convention, via `copy_current_to_old`)
  - Refs: `util/generate_velocities.cc`, `math/random.h`
- [x] **Unit-conversion audit (topology parsing)** — audited every conversion in
  `gromos-io/src/topology.rs` line-by-line against `in_topology.cc` (lines 854, 928, 1055):
  `BONDANGLEBENDTYPE` (CHT ×(180/π)², T0 deg→rad), `TORSDIHEDRALTYPE` (PD deg→rad + cos(PD)),
  `IMPDIHEDRALTYPE` (CQ ×(180/π)², Q0 deg→rad), `BONDSTRETCHTYPE`/LJ/charges (no conversion) —
  **all bit-correct, no bugs found**. Locked in via `test_parse_dihedral_and_improper_type_conversions`
  and the extended `test_parse_cg16_topology` angle assertions.
  **Decision: `BONDANGLETYPE`/`HARMBONDANGLETYPE`/`DIHEDRALTYPE` intentionally unsupported** —
  these are GROMOS96-era split/legacy encodings of the exact same data `BONDANGLEBENDTYPE`/
  `TORSDIHEDRALTYPE` already carry in unified form (kept in gromosXX only as a back-compat shim).
  No `.top` in the repo's corpus uses them; gromos-rs treats them like `TITLE`/`TOPVERSION`/
  `PHYSICALCONSTANTS` — silently ignored, by design, with no fallback dispatch to maintain.

**1.2 — Constraints (code exists, wire + test)**
- [ ] **SETTLE** for rigid water — analytical 3-site solver, O(N_water), single-pass; code exists, not wired.
  - Ref: `algorithm/constraints/settle.cc` (Miyamoto & Kollman 1992)
- [ ] **LINCS** — linear constraint solver (recursion order param); better for long chains; code exists, not wired.
  - Ref: `algorithm/constraints/lincs.cc`
- [ ] **COM rotation removal** — currently only translational. Need L = Σ mᵢ·rᵢ×(vᵢ−v_com), inertia I (3×3),
  ω = I⁻¹·L, then vᵢ −= ω×rᵢ. Ref: `algorithm/constraints/remove_com_motion.cc`

**1.3 — Thermostat**
- [ ] **Nosé-Hoover** — single bath: ζ̇ = (1/τ²)(T/T₀−1), scale = 1−ζ·dt; chain variant: M coupled baths.
  Code exists, not wired/tested. Ref: `algorithm/temperature/nosehoover_thermostat.cc`

**1.4 — Boundary**
- [ ] **Triclinic box** — code exists in `math.rs` but md.rs never creates Triclinic periodicity.
  - [ ] Wire GENBOX box_type into periodicity selection
  - [ ] Test with truncated octahedron / non-rectangular boxes

**1.5 — Advanced sampling (bigger; code exists, untested)**
- [ ] **EDS** — V_mixed = −1/β·ln(Σ exp(−β(Eᵢ−eir_i))); per-state force eval + blending; AEDS emax/emin.
  Ref: `algorithm/integration/eds.cc`
- [ ] **GaMD** — V_boost = k·(V−E_threshold)² when V>E_threshold; Welford running stats; dihedral/total/dual.
  Ref: `algorithm/integration/gamd.cc`
- [ ] **FEP / TI** — K(λ) = (1−λ)K_A + λK_B; ∂V/∂λ; soft-core LJ. Perturbed bond forces exist, need testing.
  Ref: `interaction/bonded/perturbed_*.cc`
- [ ] **REMD** (large) — MPI parallel tempering; Δ = (β₁−β₂)(E₁−E₂), accept if rand < exp(−Δ); feature-gated MPI.
  Ref: `algorithm/integration/replicaExchange/`

### Priority 2 — Analysis foundations (correlated with P1; the no-duplication layer)
The gromosPlsPls facade built on gromosXX primitives. Order matters: foundations before consumers.

**2.1 — Atom selection (gromosXX primitives first, gromos++ facade on top)**
- [ ] Solidify queryable per-atom metadata in `gromos-core` (`topology.rs`, `selection.rs`):
  molecule-membership map, residue number/name **including solvent** (today solvent = residue_nr 0 /
  name "SOLV", `topology.rs:570-571`), atom name/type/mass as flat queryable data.
- [ ] Fix the known `AtomSelection` gaps: `parse_solvent` ignores its spec & grabs all solvent;
  `parse_molecule` rejects non-first molecules; name/residue search only covers solute.
- [ ] Build the gromos++ `AtomSpecifier` grammar as a facade over the primitives (molecule ranges,
  residue by number+name, atom by name/type/number, wildcards). Ref:
  `.local/gromosPlsPls/gromos++/src/utils/AtomSpecifier.{h,cc}`, `Neighbours.{h,cc}`.
  (Grammar-depth for the first slice TBD — sequence with first real consumer.)

**2.2 — Extract shared analysis infra into the lower crates (kills duplication)**
- [ ] Rotational fit (Kabsch/SVD) replacing `simple_rotation_fit` (rmsd.rs:161 only re-centers, no rotation).
  Ref: `.local/gromosPlsPls/gromos++/src/fit/RotationalFit.cc`, `Reference.cc`.
- [ ] Statistics + error-estimate (block averaging, `ee()`) per gromos++ `gmath/Stat`.
- [ ] PBC gathering, shared across programs.
- [ ] Single-point energy entry point so `ener.rs` calls `gromos-forces` instead of hardcoding LJ σ/ε
  (requires `gromos-analysis` to depend on the engine crates — the no-duplication change).

**2.3 — Make the blocked/real programs tractable & verifiable** (once 2.1/2.2 land)
- [ ] `rmsd` — real rotational fit
- [ ] `ext_ti_ana` — real dH/dλ parsing + reference comparison (today falls back to synthetic data)
- [ ] `nhoparam` — **rewrite to the actual gromos++ algorithm** (NMR N-H order parameters S² with
  rotational fit); current code computes Nosé-Hoover thermostat params — wrong algorithm entirely.
  Ref: `.local/gromosPlsPls/gromos++/programs/nhoparam.cc`

**2.4 — Clean up other discovered stubs**
- [ ] `visco` (placeholder formula, not Green-Kubo), `frameout` (hardcoded CA/ALA), `amber2gromos`
  (ignores input), `ener` (hardcoded LJ), `sasa_hasel`, `dssp` (no real H-bonding), `solute_entropy`.

### Priority 3 — py-gromos API & education
Design reference: `.local/polars` (Python API surface, docstrings, method chaining, pyo3 wrapping).
Focus: compositional API for **system/topology construction**.

- Phase 1 — Rust bindings (pyo3-gromos). Done: compositional Simulation API, AlgorithmSequence API,
  Python reference tests (62 passed), `.pyi` stubs.
  - [ ] Expose ForceField evaluation (single-point energy/force)
  - [ ] Expose SHAKE / constraint info
  - [ ] Expose energy decomposition (bonded, LJ, CRF, kinetic, pressure)
  - [ ] Study Polars pyo3 patterns (`PyDataFrame`/`PyExpr`/`PyLazyFrame`): `.local/polars/py-polars/src/`
- Phase 2 — Python API (py-gromos). Done: Topology/Configuration/Simulation wrappers + numpy interop.
  - [ ] Method chaining: `sim.run(steps=1000).energies().plot()`
  - [ ] Energy timeseries as DataFrame (Polars/pandas interop)
  - [ ] `md_runners.py` simplify; `analysis.py` expose gromos-analysis to Python
  - [ ] Rich `__repr__` / `_repr_html_` for Jupyter (Topology, Configuration, Energy)
- Phase 3 — Notebooks & education
  - [ ] Rewrite `py-gromos/notebooks/` (01 inspect+single-point energy; 02 short MD + energy conservation;
    03 NVE/NVT/NPT comparison)
  - [ ] Rewrite `py-gromos/examples/` (17 scripts) on the new API
  - [ ] Fix `test_basic.py`, `test_advanced_features.py`; verify `maturin develop` builds + tests pass

### Priority 4 — Code quality & consistency (last)
- [ ] Clippy (~390: gromos-forces 89, gromos-integrators 77, gromos-io 31, gromos-core 15).
  `cargo clippy --fix --workspace` for auto-fixables, then review manually.
- [ ] Replace bare `unwrap()` in non-test code with `.expect("msg")` or `?`
- [ ] Add missing `#[test]`: constraints (SHAKE — currently 0), improper dihedral
- [ ] Split large files: `nonbonded.rs` (~1500 LOC), `bonded.rs` (~1300 LOC), `gromos-io/topology.rs` (~1200 LOC)
- [ ] Unify CLI error types (`Result<T,String>` → enum); audit `pub` visibility
- [ ] **Benchmarking infra**: save baseline `cargo bench --workspace -- --save-baseline v0.1`; add
  end-to-end MD step / pairlist / SHAKE / bonded benches; consider CI regression tracking; document in
  CONTRIBUTING.md. Existing: `nonbonded_bench`, `math_bench`, `thermostat_bench`, `io_bench`, `scripts/benchmark.sh`.

### Cross-cutting (do continuously) — minimal reference tests
Full tutorial t_01–t_06 end-to-end runs are **deferred to last** (compute-limited). The always-available
substitute: hand-craft **minimal** reference systems inspired by the tutorials and diff gromos-rs output
against the C++ reference, mirroring `crates/gromos-md/tests/gromosXX_references/` (generated via
`run_references.py` against the md++ binary). Stand up an analogous harness for analysis/tools (none exists
today). Every P1 physics feature and every P2 program lands with a minimal reference test.

### Parked / blocked
- `ext_ti_ana`, `nhoparam` — blocked behind Priority 2 (selection + fit + stats + references).
- Tutorials t_01–t_06 end-to-end — compute-limited; replaced near-term by minimal reference tests.

## Key Files

```
crates/gromos-md/src/bin/md.rs               — main MD driver, CLI, simulation setup
crates/gromos-core/src/algorithm.rs          — Algorithm trait, AlgorithmSequence
crates/gromos-core/src/math.rs               — Vec3, BoundaryCondition (Vacuum/Rectangular/Triclinic)
crates/gromos-core/src/topology.rs           — Topology struct
crates/gromos-core/src/selection.rs          — AtomSelection (P2: extend to AtomSpecifier facade)
crates/gromos-forces/src/bonded.rs           — all bonded force calculations
crates/gromos-forces/src/nonbonded.rs        — LJ+CRF, rf_excluded, pairlist loops
crates/gromos-forces/src/electrostatics.rs   — CRF/PME parameters
crates/gromos-integrators/src/algorithms/    — Forcefield, LeapFrog*, Shake, Temperature, Energy
crates/gromos-integrators/src/constraints.rs — SHAKE, SETTLE, LINCS
crates/gromos-integrators/src/thermostats.rs — Berendsen, Nosé-Hoover, Andersen
crates/gromos-integrators/src/barostats.rs   — Berendsen, Parrinello-Rahman
crates/gromos-io/src/topology.rs             — topology file parser
crates/gromos-io/src/coordinate.rs           — coordinate/GENBOX parser
crates/gromos-io/src/imd.rs                  — IMD parameter file parser
crates/gromos-tools/src/bin/topology/        — make_top, com_top, check_top, etc.
crates/gromos-tools/src/bin/box/             — sim_box, build_box, ion, etc.
crates/gromos-analysis/src/bin/structural/rmsd.rs — simple_rotation_fit (P2: replace w/ Kabsch)
crates/gromos-md/tests/test_gromosXX_references.rs — integration tests vs gromosXX
crates/gromos-md/tests/run_references.py           — generate gromosXX reference data
crates/gromos-md/tests/gromosXX_references/        — reference input + expected output
```

## How to Test

```sh
# Build
cargo build --release --bin md

# Run integration tests against gromosXX references
cargo test -p gromos-md --test test_gromosXX_references

# Run including known-failing systems
cargo test -p gromos-md --test test_gromosXX_references -- --include-ignored

# Run a specific system
cargo test -p gromos-md --test test_gromosXX_references -- pair_lj --exact

# Regenerate gromosXX reference data (requires gromosXX md++ binary)
python3 crates/gromos-md/tests/run_references.py --md-binary .local/gromosXX/md++/build/program/md

# Add a new reference system:
# 1. Add to SYSTEMS list in crates/gromos-md/tests/run_references.py
# 2. Run run_references.py to generate expected/ data
# 3. Add ref_test!(name, "dir") in crates/gromos-md/tests/test_gromosXX_references.rs
```
