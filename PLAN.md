# gromos-rs — Status & Plan

Focus: `cargo build --release --bin md`
On commit: update CHANGELOG.md and Cargo.toml version.
DON'T modify `gromosXX_references/*/expected/` — those are ground truth.

References:
- gromosXX source: `.local/gromosXX/md++/src`
- Tutorials: `.local/gromos_tutorial_livecoms/tutorial_files`
- Theory: `.local/doc/gromos_book`
- Force fields: `.local/gromosXX/forcefields`

Doc style: Rust → KaTeX + `[^label]` footnotes; Python → NumPy docstrings + `.. math::`


## Architecture

| Crate | Role |
|-------|------|
| gromos-core | Types, math (Vec3), BoundaryCondition, Topology, Configuration, Algorithm trait |
| gromos-forces | Bonded, nonbonded (LJ+CRF), pairlist, restraints, PME/QM-MM (stubs) |
| gromos-integrators | LeapFrog, SteepestDescent, SHAKE/SETTLE/LINCS, thermostats, barostats, EDS/GaMD/REMD/FEP |
| gromos-io | All file I/O: topology, coordinates, IMD, trajectories, energy, MTB/IFP, PDB, posres |
| gromos-md | Simulation engines: md, md_mpi, remd, eds, gamd, mdf |
| gromos-tools | System construction: make_top, sim_box, pdb2g96, mk_script, ion, prep_posres, ... |
| gromos-analysis | Trajectory analysis: rmsd, rdf, ene_ana, hbond, ... |
| pyo3-gromos | Python bindings (PyO3) |
| py-gromos | Python package (separate maturin build) |

**Design rule:** All file format parsing/writing lives exclusively in `gromos-io`. No duplication across crates.

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

**27 of 27 tests pass.** All levels fully passing.

## What Works

### gromos-core
- Boundary conditions: Vacuum, Rectangular (minimum image), Triclinic (defined but not wired)
- Topology struct with solvent expansion: NSM auto-computed, chargegroups from CGC codes, intra-molecular exclusions

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

## TODO

### PRIORITY — System Construction Programs (tutorial-ready)
Goal: make the molecular system construction pipeline work end-to-end so users can follow
the realistic tutorials at `.local/gromos_tutorial_livecoms/tutorial_files/`.

Programs needed for tutorials t_01 through t_06 (in pipeline order):

#### System building (critical path)
Done: make_top, com_top, check_top, pdb2g96, sim_box, ion, mk_script, make_pt_top, prep_posres

- [ ] **build_box** — PARTIAL, hardcodes water/3-atom; needs real solute file reading
  - Must: read solute coordinates, replicate to fill box at target density
  - Tutorial: used for creating pure solvent boxes

#### Analysis programs (needed by tutorials)
Done: ene_ana, rmsd, rmsf, rdf, hbond, frameout, tser, trs_ana, bar, reweight, filter

- [ ] **ext_ti_ana** — PARTIAL, falls back to placeholder if no input
- [ ] **nhoparam** — verify implementation status (NHO order parameters, t_01)

#### Validation plan
- [ ] Run tutorial t_01 (protein MD) end-to-end with gromos-rs binaries
  - Pipeline: make_top → com_top → check_top → pdb2g96 → sim_box → ion → md → ene_ana/rmsd/rmsf
- [ ] Compare output of each step against gromosXX++ reference output
- [ ] Document which tutorial steps work and which are blocked



### TODO — Triclinic box support (lower priority)
- [ ] Code exists in `math.rs` but md.rs never creates Triclinic periodicity
  - [ ] Wire GENBOX box_type into periodicity selection
  - [ ] Test with truncated octahedron or other non-rectangular boxes

### TODO — Unit conversion audit (topology parsing)
- [ ] HARMBONDANGLETYPE: if/when parsed, must also convert CHT — ×(180/π)²
- Note: IMPDIHEDRALTYPE and BONDANGLEBENDTYPE already converted
- gromosXX converts ALL angle force constants at parse time (in_topology.cc)
- No conversion needed for: CT (cosine angle K), CP (dihedral K), bond K, LJ, charges
- Reference: gromosXX in_topology.cc lines 854, 928, 1055

### TODO — Code quality & consistency
- [ ] Fix clippy warnings (~390 total: gromos-forces 89, gromos-integrators 77, gromos-io 31, gromos-core 15)
  - Mechanical fixes (auto-fixable): borrowed refs, double parens, legacy constants, assign ops
  - Manual: unused variables, dead imports, loop indexing → iterators
- [ ] Run `cargo clippy --fix --workspace` for auto-fixable warnings, review manually after
- [ ] Replace bare `unwrap()` in non-test code with `.expect("msg")` or `?` (2 in CLI arg parsing)
- [ ] Add missing `#[test]` for constraints (SHAKE unit tests — currently 0)
- [ ] Add `#[test]` for improper dihedral (bonded tests missing this)
- [ ] Review large files for possible splitting:
  - `nonbonded.rs` (~1500 LOC) — SIMD / parallel / serial inner loops could be submodules
  - `bonded.rs` (~1300 LOC) — each force type could be its own file
  - `topology.rs` (gromos-io, ~1200 LOC) — parser per block type could be submodules
- [ ] Consistent error types: unify `Result<T, String>` in CLI → proper error enum
- [ ] Audit `pub` visibility — many internal functions are unnecessarily public

### PRIORITY — py-gromos API improvement
Goal: a polished, ergonomic Python API for education and notebooks.
Design reference: `.local/polars` — study its Python API surface, docstrings, method chaining,
DataFrame-style ergonomics, and how it wraps Rust internals via pyo3.

#### Phase 1 — Rust bindings (pyo3-gromos)

Done: Compositional Simulation API (Topology, Configuration, InputParameters, Simulation),
AlgorithmSequence API (preset constructors, sequence manipulation, descriptor resolution),
Python reference tests (62 passed), type stubs (.pyi) for static analysis.

##### TODO — Additional binding types
- [ ] Expose ForceField evaluation (single-point energy/force calculation)
- [ ] Expose SHAKE / constraint info
- [ ] Expose energy decomposition (bonded, LJ, CRF, kinetic, pressure)
- [ ] Study Polars' pyo3 patterns: `PyDataFrame`, `PyExpr`, `PyLazyFrame` wrappers
  - Path: `.local/polars/py-polars/src/` — how they wrap Rust types into Python classes
  - Learn from: `__repr__`, `_repr_html_`, method chaining, `@staticmethod` constructors

#### Phase 2 — Python API (py-gromos)
Done: Topology, Configuration, Simulation wrappers with `__repr__` and numpy interop.

- [ ] Method chaining where natural: `sim.run(steps=1000).energies().plot()`
- [ ] Energy timeseries as DataFrame (Polars or pandas interop)
- [ ] `md_runners.py` — review and simplify (currently wraps CLI `md` binary)
- [ ] `analysis.py` — expose gromos-analysis functions to Python
- [ ] Rich `__repr__` / `_repr_html_` for Jupyter display of Topology, Configuration, Energy

#### Phase 3 — Notebooks & education
- [ ] Rewrite `py-gromos/notebooks/` as clean educational notebooks using the new API
  - [ ] 01: Load topology + coordinates, inspect atoms/bonds, compute single-point energy
  - [ ] 02: Run short MD, plot energy conservation, visualize trajectory
  - [ ] 03: Compare NVE vs NVT vs NPT ensembles, thermostat/barostat effects
- [ ] Rewrite `py-gromos/examples/` to use the new API (currently 17 example scripts)
- [ ] Fix/update Python tests (`test_basic.py`, `test_advanced_features.py`)
- [ ] Verify `maturin develop` builds and tests pass

### TODO — NTIVEL=1 velocity generation (small)
- [ ] Implement Maxwell-Boltzmann velocity generation (NTIVEL=1)
  - [ ] Read NTIVEL, IG (seed), TEMPI (temperature) from INITIALISE block (already parsed in `imd.rs`)
  - [ ] Implement MT19937 RNG matching GSL's `gsl_rng_mt19937`
  - [ ] Implement Gaussian distribution matching GSL's `gsl_ran_gaussian` (polar Box-Muller)
  - [ ] For each atom: σ = sqrt(k_B·T / m_i), v_i = gaussian(σ) for x,y,z
  - [ ] k_Boltzmann = 0.00831441 kJ/(mol·K)
  - [ ] Store in both current().vel and old().vel (gromosXX convention)
  - [ ] gromosXX ref: `util/generate_velocities.cc`, `math/random.h`
  - [ ] Currently worked around with pre-generated velocities in .conf files
  - [ ] Important for reproducibility: users expect NTIVEL=1 + seed to produce deterministic runs

### TODO — Benchmarking infrastructure
Goal: track performance evolution over time, have a reference baseline for regressions.

Existing benchmarks (Criterion):
- `crates/gromos-forces/benches/nonbonded_bench.rs` — nonbonded innerloop
- `crates/gromos-core/benches/math_bench.rs` — Vec3 / math ops
- `crates/gromos-integrators/benches/thermostat_bench.rs` — thermostat scaling
- `crates/gromos-io/benches/io_bench.rs` — topology/coordinate parsing
- `scripts/benchmark.sh` — wrapper for `cargo bench --workspace`

- [ ] Run current benchmarks and save a baseline: `cargo bench --workspace -- --save-baseline v0.1`
- [ ] Add end-to-end MD step benchmark (full forcefield + integrator for water_216_box)
- [ ] Add pairlist construction benchmark (chargegroup-based and atom-based)
- [ ] Add SHAKE constraint benchmark (water_216 system)
- [ ] Add bonded forces benchmark (aladip_solvated — all bonded types)
- [ ] Consider `criterion-compare` or CI integration for automated regression tracking
- [ ] Add Python benchmark: pyo3-gromos single-point energy call overhead
- [ ] Document how to run and compare: `cargo bench -- --baseline v0.1` in CONTRIBUTING.md

### Known Gaps — Future Features

#### TODO — COM rotation removal (medium)
- [ ] Currently only translational COM removal is implemented
- [ ] Need: angular momentum L = Σ mᵢ × rᵢ × (vᵢ - v_com), inertia tensor I (3×3), ω = I⁻¹·L
- [ ] Remove rotation: vᵢ -= ω × rᵢ
- [ ] gromosXX ref: `algorithm/constraints/remove_com_motion.cc`

#### TODO — SETTLE constraints for rigid water (small)
- [ ] Analytical solver for 3-site rigid water (SPC/TIP3P) — no iteration
- [ ] Much faster than SHAKE for water: O(N_water) per step, single-pass
- [ ] Code exists in gromos-rs but not wired
- [ ] gromosXX ref: `algorithm/constraints/settle.cc` (Miyamoto & Kollman 1992)

#### TODO — LINCS constraints (medium)
- [ ] Linear constraint solver using recursion order parameter
- [ ] Alternative to SHAKE — better for long chains, parallelizable
- [ ] Code exists in gromos-rs but not wired
- [ ] gromosXX ref: `algorithm/constraints/lincs.cc`

#### TODO — Nosé-Hoover thermostat (medium)
- [ ] Single bath: ζ̇ = (1/τ²)(T/T₀ - 1), scale = 1 - ζ·dt
- [ ] Chain variant: M coupled thermostats for ergodicity
- [ ] Code exists in gromos-rs but not wired or tested
- [ ] gromosXX ref: `algorithm/temperature/nosehoover_thermostat.cc`

#### TODO — Triclinic box support (medium)
- [ ] Code exists in `math.rs` (Triclinic periodicity) but md.rs never creates it
- [ ] Wire GENBOX box_type into periodicity selection
- [ ] Test with truncated octahedron or other non-rectangular boxes

#### TODO — EDS — Enveloping Distribution Sampling (medium)
- [ ] Multi-state Boltzmann averaging: V_mixed = -1/β · ln(Σ exp(-β(Eᵢ - eir_i)))
- [ ] Per-state force evaluation + blending
- [ ] Adaptive EDS (AEDS) with emax/emin boundaries
- [ ] Code exists in gromos-rs, not tested against references
- [ ] gromosXX ref: `algorithm/integration/eds.cc`

#### TODO — GaMD — Gaussian Accelerated MD (medium)
- [ ] Boost potential: V_boost = k·(V - E_threshold)² when V > E_threshold
- [ ] Welford algorithm for running statistics (mean, variance)
- [ ] Three boost forms: dihedral-only, total-potential, dual
- [ ] Code exists in gromos-rs, not tested against references
- [ ] gromosXX ref: `algorithm/integration/gamd.cc`

#### TODO — FEP / TI — Free Energy Perturbation (medium)
- [ ] Lambda interpolation: K(λ) = (1-λ)K_A + λK_B for all interaction types
- [ ] TI derivative: ∂V/∂λ for thermodynamic integration
- [ ] Soft-core LJ to avoid singularities during atom appearance/disappearance
- [ ] Code exists in gromos-rs (perturbed bond forces), needs testing
- [ ] gromosXX ref: `interaction/bonded/perturbed_*.cc`

#### TODO — REMD — Replica Exchange MD (large)
- [ ] MPI-based parallel tempering across temperature/lambda ladders
- [ ] Exchange acceptance: Δ = (β₁ - β₂)(E₁ - E₂), accept if rand < exp(-Δ)
- [ ] Requires MPI infrastructure (feature-gated)
- [ ] gromosXX ref: `algorithm/integration/replicaExchange/`

## Key Files

```
crates/gromos-md/src/bin/md.rs               — main MD driver, CLI, simulation setup
crates/gromos-core/src/algorithm.rs          — Algorithm trait, AlgorithmSequence
crates/gromos-core/src/math.rs               — Vec3, BoundaryCondition (Vacuum/Rectangular/Triclinic)
crates/gromos-core/src/topology.rs           — Topology struct
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
crates/gromos-md/tests/test_gromosXX_references.rs — integration tests vs gromosXX
crates/gromos-md/tests/run_references.py           — generate gromosXX reference data
crates/gromos-md/tests/gromosXX_references/        — reference input + expected output
```

## How to Test

```sh
# Build
cargo build --release --bin md

# Run integration tests against gromosXX references (14 active, 5 ignored)
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

