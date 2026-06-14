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
- **`FUTURE.md`** — architectural bets where Rust can *overtake* the C++ (SoA core, single
  multi-backend kernel, O(N) cell-list pairlist, solute/solvent representation-vs-role refactor,
  compositional py-gromos topology, GPU/SIMD/determinism, **QM/MM + in-process ML potentials (Dim
  12)**, **Martini/CG ecosystem bridge (Dim 13)**, **the unifying layered architecture (Dim 14) —
  one data+compute core shared by engine and analysis facade**), **plus the differential-audit
  findings** (known gromosXX bugs *not* to port, and one live divergence already in our tree — P1.4).
  PLAN.md = parity & near-term execution; FUTURE.md = where we diverge on purpose.

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

**The unifying architecture (FUTURE.md Dim 14 — the frame for everything below).** The two C++
codebases are *two architectures that don't share a core* (gromosXX = pipeline of mutating algorithms
over `(Topology, Configuration, Simulation)`; gromos++ = ~107 frame-loop `main()`s over a *separate*
`gcore` that re-implements the physics — that re-implementation **is** the 79k LOC). gromos-rs must be
**one machine**: a layered design where the data core (L0: SoA + instancing + roles + spatial-index)
and the *pure compute layer* (L1: `PotentialProvider` for forces, `Observable` for analysis quantities)
are **shared by both the MD engine and the analysis facade** — so an MD step (L0+L1+L2 steppers, loop
over time) and an analysis (L0+L1 only, loop over frames) are the same machinery. **Lock now** (they
shape every refactor): (a) **one owned `(Topology, State)` core, never a second `System`**; (b) **L1
purity — every energy/force/observable is computed in one place, called by both engine and analysis**
(the moment `ener.rs` re-implements LJ, you've recreated gromos++); (c) a **small stable trait
taxonomy** (`PotentialProvider`, `Observable`, `Stepper`, `SpatialIndex`, `Reader`/`Writer`); (d) one
**state-history model** (the `current()`/`old()` question). The QM/MM+ML prerequisites below are the
Dim-12 specialization of this.

**Architectural prerequisites for QM/MM + ML potentials (FUTURE.md Dim 12 — decide on the *next*
force/pairlist refactor, not later).** Dim 12 is the flagship overtake (in-process ML potentials with
no Python boundary; QM/MM as an additive force provider). It costs ~nothing to design for now and is
very expensive to retrofit, so near-term refactors must honor these shapes even before any QM/ML code
exists:
1. **Force evaluation = additive `PotentialProvider`s** (data-in / scattered-forces-out, `&mut self`
   for stateful workers, fallible, extensible return for uncertainty) — *classical* LJ+CRF/bonded
   implement the same trait, so QM/ML is later just another provider, not a fork. (Reinforces the
   single-kernel goal, FUTURE Dim 2.)
2. **The pairlist/cell list (P1.6) is a query-based spatial-index *service*** that can emit the MD
   pairlist, an ML radial graph, and the QM-zone gather from one structure — not a single-purpose
   MD pairlist. (FUTURE Dim 9e.)
3. **A typed units boundary** (kJ/mol·nm ↔ Hartree/Bohr ↔ eV/Å) introduced with the engine, so QM/ML
   unit conversions are checked, not silent.
4. **Region/role membership is a data-model attribute** and "promote a water into the QM zone" is the
   Dim 10 de-instancing op — so the solute/solvent refactor (FUTURE Dim 10) is a Dim 12 prerequisite.
5. **A second, tolerance-based test tier** (validated against reference energies/forces, not
   bit-for-bit) for the non-deterministic QM/ML paths, alongside the classical bit-for-bit suite.

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
| 2   | nacl_1water_settle | 5   | SETTLE (analytical rigid water)      | **PASS** |
| 2   | nacl_1water_lincs | 5    | LINCS (solvent)                      | **PASS** |
| 2   | nacl_3water_box  | 11    | multiple solvent + solute-solvent pairlist | **PASS** |
| 2   | water_3_box_twinrange | 9 | twin-range pairlist (RCUTP<RCUTL, NSNB=5) | **PASS** |
| 2   | water_10_box     | 32    | 2 ions + 10 SPC, positions away from cutoff | **PASS** |
| 2   | nacl_3water_cutoff | 11  | nacl_3water near cutoff boundary     | **PASS** |
| 2   | nacl_water_box   | 62    | ion-water RF in PBC                  | **PASS** |
| 2   | nacl_water_box_shifted | 62 | nacl_water_box with perturbed positions | **PASS** |
| 3   | water_216_box    | 648   | bulk NVE, pairlist, virial           | **PASS** |
| 3   | water_216_box_com| 648   | bulk NVE + COM removal (NTICOM=1, NSCM=10) | **PASS** |
| 3   | water_216_box_com_rot | 648 | COM translation+rotation removal (NTICOM=2, NSCM=-10) | **PASS** |
| 3   | water_216_nvt    | 648   | Berendsen thermostat                 | **PASS** |
| 3   | water_216_nvt_nosehoover | 648 | Nosé-Hoover thermostat (single NHC) | **PASS** |
| 3   | water_216_nvt_nhc_chain | 648 | Nosé-Hoover-Chain (3 chains)        | **PASS** |
| 3   | water_216_npt    | 648   | Berendsen barostat                   | **PASS** |
| 4   | aladip_vacuum_lincs | 12 | LINCS (solute, NTC=2)               | **PASS** |
| 4   | aladip_solvated  | 72    | SHAKE + solute-solvent               | **PASS** |
| 4   | aladip_vacuum_em | 12    | steepest descent EM, vacuum          | **PASS** |
| 4   | aladip_vacuum_em_shake | 12 | SD EM + SHAKE, vacuum             | **PASS** |
| 4   | aladip_solvated_em_noshake | 72 | SD EM, solvated, no SHAKE      | **PASS** |
| 4   | aladip_solvated_em_shake | 72 | SD EM + SHAKE, solvated          | **PASS** |
| 4   | aladip_solvated_em_posres | 72 | SD EM + position restraints     | **PASS** |
| 4   | aladip_solvated_em | 72  | SD EM + SHAKE + posres, solvated    | **PASS** |

**34 of 34 tests pass.** All levels fully passing.

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
  - ⚠️ **O(N²) only** — `StandardPairlistAlgorithm` loops all chargegroup-COG pairs. Correct and
    bit-for-bit, fine for reference-size systems, but no cell/grid list exists. gromosXX ships an
    O(N) `grid_cell`/`extended_grid` pairlist; ours does not. This is the biggest scaling gap —
    see new roadmap item **P1.6** and FUTURE.md Dim 9. (When porting the grid path, do **not** port
    the Martina-flagged solute/solvent misclassification bug at `extended_grid_pairlist_algorithm.cc:1309`.)
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
- [x] **SETTLE** for rigid water — analytical 3-site solver, O(N_water), single-pass; wired and reference-tested (`nacl_1water_settle`).
  - Ref: `algorithm/constraints/settle.cc` (Miyamoto & Kollman 1992)
- [x] **LINCS** — linear constraint solver (recursion order param); better for long chains; wired and reference-tested (`nacl_1water_lincs`, `aladip_vacuum_lincs`).
  - Ref: `algorithm/constraints/lincs.cc`
- [x] **COM rotation removal** — wired and reference-tested (`water_216_box_com_rot`).
  - Fix: use minimum-image (put_into_box) wrapped positions inside `remove_com_rotation`, mirroring gromosXX's `gather_chargegroups` init convention.
  - Fix: suppress periodic COM rotation removal in PBC (gromosXX `configuration.cc:555-560`).
  - Ref: `algorithm/constraints/remove_com_motion.cc`

**1.3 — Thermostat**
- [x] **Nosé-Hoover** — single NHC (algorithm=1) and chain NHC (algorithm=N≥2); reference-tested
  (`water_216_nvt_nosehoover`, `water_216_nvt_nhc_chain`). Both pass against gromosXX.
  - Ref: `algorithm/temperature/nosehoover_thermostat.cc`
  - Note: gromosXX has a *“small flexible constraints hack!”* (`nosehoover_thermostat.cc:129,185`,
    also `berendsen_thermostat.cc:105`) in the DOF/coupling path. Not ported — only applies
    with flex-shake, which we don’t support. **Documented as a deliberate GROMOS quirk.**
  - NTINHT convention documented: 0=read from file, 1=generate from scratch (opposite of intuition).
  - IMD parser bug fixed: algorithm codes now match gromosXX (0=Berendsen, 1=NHC, N=chain length).

**1.4 — Boundary**
- [x] **Triclinic box** — code exists in `math.rs` but md.rs never creates Triclinic periodicity.
  - ⚠️ **KNOWN LIVE DIVERGENCE (FUTURE.md Dim 11, finding #1) — RESOLVED.** Our `Triclinic::nearest_image`
    (`gromos-core/src/math.rs:90`) used the textbook fractional-coordinate `frac - frac.round()`.
    gromosXX (`math/boundary_implementation.cc:285-318`) uses an **iterative `while`-loop reduction
    in z→y→x order** over the lower-triangular box. These are **not equivalent for strongly
    triclinic cells** — fractional `round()` is not always the Cartesian nearest image. Ported the
    `while`-loop z→y→x version (see below); `aladip_trunc_oct` now passes.
  - [x] Add a triclinic gromosXX reference (truncated octahedron) — *before* wiring
    - `aladip_trunc_oct` (aladip.topo/aladip.conf, NTB=-1, RCUTP/RCUTL=0.8/0.9 — the
      box's inscribed-sphere check caps RCUTL around ~1.33 nm for this 3.767 nm cell).
      Reference data generated; Rust test added as `ignore: aladip_trunc_oct` —
      confirmed failing (E_total 133.62 vs gromosXX 132.07), as predicted: md.rs
      ignores NTB and treats the box as Rectangular without the truncoct_triclinic
      rotation gromosXX applies on read.
  - [x] Replace `Triclinic::nearest_image` with the gromosXX while-loop z→y→x reduction
    - Ported `boundary_implementation.cc:285-318` verbatim into `math.rs` (assumes the
      GROMOS lower-triangular box convention). Added unit tests demonstrating that the
      while-loop and the old `frac.round()` reduction diverge by a full lattice vector at
      an exact half-box tie (FUTURE.md Dim 11 #1), and that large displacements still
      reduce correctly. `aladip_trunc_oct` still fails identically (133.62 vs 132.07) —
      expected, since md.rs doesn't construct a triclinic box yet (item below).
  - [x] Wire GENBOX box_type into periodicity selection (incl. truncoct_triclinic_box
    conversion + truncoct_triclinic rotation of positions/velocities on read for NTB=-1)
    - Ported `math::truncoct_triclinic_rotmat`/`truncoct_triclinic_box`/`truncoct_triclinic`
      (`math/transformation.cc`) into `math.rs` as faithful 3-function APIs taking a
      `forward: bool`. Added `Box::truncated_octahedral()`. md.rs now: for NTB=-1, converts
      the cubic BOX block via `truncoct_triclinic_box(.., true)`, rotates positions/velocities
      via `truncoct_triclinic(.., true)`, builds `Periodicity::Triclinic` from the resulting
      lower-triangular box (forcefield.rs's periodicity refresh now matches
      `TruncatedOctahedral` too), and rotates `FREEFORCERED`/`CONSFORCERED` back to the cube
      frame via `truncoct_triclinic_rotmat(false)` on output (matches
      `out_configuration.cc::_print_forcered`).
  - [x] Test with truncated octahedron / non-rectangular boxes — un-ignore `aladip_trunc_oct`
    - **Found and fixed a rotation-matrix bug** while debugging: gromosXX's `product(rot, v)`
      (`gmath.h`) computes `rot^T * v`, not `rot * v` (its `GenericMatrix` constructor stores
      the 3 ctor-arg `Vec`s as *rows*, but `product` contracts over the *first* matrix index).
      Our initial port applied `rot * v` for `forward` and `rot^T * v` for `!forward` —
      exactly backwards. Fixed by swapping the two branches; added a regression test pinned
      to a gromosXX debug-build reference value for `aladip.conf` atom 21
      (`test_truncoct_triclinic_forward_matches_gromosxx_reference`). `aladip_trunc_oct` now
      passes (energies + forces match gromosXX to tolerance); full reference suite
      (35/35) still green.

**1.5 — Scaling: O(N) cell-list pairlist (the biggest scaling gap)** — see FUTURE.md Dim 9 for the
full design; the near-term, bit-for-bit-safe slice:
- [x] **Charge-group-aware cell (linked) list** as `CellListPairlistAlgorithm`, drop-in behind the
  existing `update<BC>()` interface, selectable by system size. Bin **chargegroup COGs** (not atoms)
  to preserve GROMOS's neutral-group cutoff + RF correctness.
  - Implemented in `gromos-core/src/pairlist.rs`: bins each chargegroup's reference position
    (solute CG → COG, solvent CG → first-atom position, matching `StandardPairlistAlgorithm`'s own
    distance conventions) into a grid that exactly tiles the box, cells ≥ `long_range_cutoff +
    skin`, with periodic-wrapped 27-cell neighbor search (deduplicated, so it stays correct even
    when `grid_dim` is 1 or 2 along an axis). For `BoxType::Rectangular` this is the real O(N) path;
    for vacuum/triclinic/truncated-octahedron boxes (where an axis-aligned grid can't be made
    periodicity-safe without extra machinery) it transparently falls back to
    `StandardPairlistAlgorithm`'s O(N²) path — still correct, just not yet accelerated. Exported
    from `gromos-core` alongside `StandardPairlistAlgorithm`. Replaces the old, unexported,
    non-chargegroup-aware `GridCellPairlistAlgorithm` stub. Not yet wired into the `md`/
    `pyo3-gromos` binaries' algorithm selection — the system-size heuristic is the remaining piece.
- [x] **Validate by set-equality** against `StandardPairlistAlgorithm` (the O(N²) oracle) on every
  reference system — identical pair set, not just identical energy — before trusting it.
  - Added `pairlist::tests::test_cell_list_matches_standard_{solvent_only,single_cell,with_solute}`:
    sort+min/max-normalize all four pair lists (`solute_short/long`, `solvent_short/long`) from
    both algorithms and assert equality. Covers pure-solvent with periodic wrap (grid_dim=4),
    the single-cell fallback (grid_dim=1), and mixed solute+solvent with exclusions (grid_dim=2,
    the trickiest neighbor-cell dedup case). Full `gromosXX_references` suite (35/35) still
    passes unchanged (still runs on `StandardPairlistAlgorithm`).
- [x] Keep O(N²) as the always-correct reference path; cell list is an accelerator, not a replacement.
  - `StandardPairlistAlgorithm` remains the only algorithm wired into `md.rs` and `pyo3-gromos`;
    `CellListPairlistAlgorithm` is an additional, independently-tested option.
- [x] Do **not** reproduce the Martina solute/solvent misclassification bug
  (`extended_grid_pairlist_algorithm.cc:1309`); classify pairs by **both** atoms' roles. (Dim 10
  removes the root cause entirely.)
  - `process_cg_pair` classifies each CG-CG pair by `cg1=min(cg_a,cg_b)` / `cg2=max(cg_a,cg_b)`
    against `n_solute_cg` — i.e. by both chargegroups' roles, not just one.
- Deferred to FUTURE.md (post-parity): spatial reorder / Z-order of atom arrays (Dim 9b),
  displacement-triggered rebuild (Dim 9c), charge groups as a first-class primitive (Dim 9d),
  triclinic/truncated-octahedron cell-list acceleration, and system-size-based algorithm
  selection in the MD binaries.

**1.6 — Restraints & special interactions** — today **only position restraints** exist
(`What Works`). gromosXX has ~24 special-interaction types in `interaction/special/`; the
scientifically essential ones each need porting + a minimal reference test, and most have a
*perturbed* variant (couples to FEP, P1.5). Per-term ground-truth available from gromosXX's own
`check/aladip*.t.cc` (see cross-cutting test note). Priority within this group:
- [ ] **Distance restraints** (`distance_restraint_interaction.cc`) — NMR NOE, the most common. +perturbed.
- [ ] **Dihedral restraints** (`dihedral_restraint_interaction.cc`) — +perturbed. (Note gromosXX caveat:
  `perturbed_dihedral_restraint_interaction.cc:152` — "may go wrong if phi0_A/phi0_B > 2π apart"; Dim 11 — second-source.)
- [ ] **Angle restraints** (`angle_restraint_interaction.cc`) — +perturbed.
- [ ] **J-value restraints** (`jvalue_restraint_interaction.cc`) — NMR ³J couplings.
- [ ] **Order-parameter restraints** (`order_parameter_restraint_interaction.cc`) — NMR S². (Two
  already-fixed bugs documented in-source at `:91,:213` — read them before porting; Dim 11.)
- [ ] **Distance-field** (`distance_field_interaction.cc`) + **local elevation**
  (`local_elevation_interaction.cc`) — enhanced sampling. +perturbed distance-field.
- Lower priority within group: RDC, X-ray, symmetry, colvar, electric-field, NEMD, adde_reweighting.
- Each lands with a minimal gromosXX reference test, mirroring the MD harness.

**1.7 — Virtual atoms** (`algorithm/virtualatoms/`, tool `addvirt_top`) — massless interaction
sites placed by geometry (e.g. aromatic centroids, lone pairs). **Primary motivation here: TI/free
energy** (P1.5) and certain force fields. Touches the data model — coordinate with the Dim 10
instancing refactor (FUTURE.md) so virtual sites aren't hard-coded around the old solute/solvent
split. Reference: `algorithm/virtualatoms/prepare_virtualatoms.cc`, `create_virtualatoms.cc`.

**1.8 — Advanced sampling (bigger; code exists, untested, OPTIONAL)**
- [ ] **EDS** — V_mixed = −1/β·ln(Σ exp(−β(Eᵢ−eir_i))); per-state force eval + blending; AEDS emax/emin.
  Ref: `algorithm/integration/eds.cc`
- [ ] **GaMD** — V_boost = k·(V−E_threshold)² when V>E_threshold; Welford running stats; dihedral/total/dual.
  Ref: `algorithm/integration/gamd.cc`
- [ ] **FEP / TI** — K(λ) = (1−λ)K_A + λK_B; ∂V/∂λ; soft-core LJ. Perturbed bond forces exist, need testing.
  Ref: `interaction/bonded/perturbed_*.cc`
  - ⚠️ **Second-source the perturbed RF self-term — do NOT transcribe the C++.** The gromosXX authors
    themselves flagged it as untrusted: `perturbed_nonbonded_term.cc:596,749` (*"Chris: CHECK! I'm
    not sure if the self-term correction is not wrong…"*) and `:1444` (*"there is a bug here from the
    previous version"*). Derive from the GROMOS book, implement to that, then diff against the C++;
    where they disagree, investigate. (FUTURE.md Dim 11 findings #3/#4.)
- [ ] **REMD** (large) — MPI parallel tempering; Δ = (β₁−β₂)(E₁−E₂), accept if rand < exp(−Δ); feature-gated MPI.
  Ref: `algorithm/integration/replicaExchange/`

### Priority 2 — Analysis foundations (correlated with P1; the no-duplication layer)
The gromosPlsPls facade built on gromosXX primitives. Order matters: foundations before consumers.

> **Architectural note (FUTURE.md Dim 10): the solute/solvent split.** gromosXX uses a rigid index
> threshold (`i >= num_solute_atoms()`); gromos++ is worse — *separate containers* (`mol()` vs
> `sol()`) forcing dual code paths in every analysis program and an `m:`/`s:` namespace split. Our
> `solvate()` is currently the worst hybrid (keeps a solvent template *but also* expands flat
> `mass/charge/iac` arrays). The target model decouples **representation** (a `MoleculeType` registry
> + instances — dedup any repeated molecule, not just water) from **role** (solvent-ness as a
> queryable attribute, not a partition), enabling cheap `promote()` of one water to solute. **Do the
> selection work below in a way that treats solvent-ness as an attribute, so it doesn't calcify
> around the old split before the Dim 10 refactor lands.**

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
- [ ] **PBC gathering / molecule unwrapping — a unify-the-logics target, not a one-liner.** This is
  the canonical place where gromosXX and gromos++ have *separate, divergent* implementations of the
  same idea (the engine's `Periodicity`/min-image vs gromos++'s `bound/` gather methods: gbond, gref,
  gcog, gltime, …). It's also notoriously bug-prone — a molecule broken across the box silently
  corrupts every downstream analysis. **Build one gathering primitive in `gromos-core` (on the
  engine's boundary code) and have both the engine and the gromos++-style tools (`gathtraj`, `inbox`,
  `unify_box`, `frameout`) consume it** — the no-duplication principle applied to the single most
  duplicated, most error-prone shared algorithm. Ref: `.local/gromosPlsPls/gromos++/src/bound/`,
  `utils/Gather.cc`.
- [ ] Single-point energy entry point so `ener.rs` calls `gromos-forces` instead of hardcoding LJ σ/ε
  (requires `gromos-analysis` to depend on the engine crates — the no-duplication change).

**2.3 — Make the blocked/real programs tractable & verifiable** (once 2.1/2.2 land)
- [ ] `rmsd` — real rotational fit
- [ ] `ext_ti_ana` — real dH/dλ parsing + reference comparison (today falls back to synthetic data)
- [ ] `nhoparam` — **rewrite to the actual gromos++ algorithm** (NMR N-H order parameters S² with
  rotational fit); current code computes Nosé-Hoover thermostat params — wrong algorithm entirely.
  Ref: `.local/gromosPlsPls/gromos++/programs/nhoparam.cc`

**2.3b — Free-energy analysis estimators** (the *analysis* side of free energy; the *engine* side —
FEP/TI, EDS, GaMD — is P1.5). A coherent family worth doing together once dH/dλ parsing (`ext_ti_ana`)
and statistics/`ee()` land:
- [ ] `bar` — Bennett Acceptance Ratio. Ref: `.local/gromosPlsPls/gromos++/programs/bar.cc`
- [ ] `ext_ti_merge` — merge/extrapolate TI curves.
- [ ] `reweight` — configurational reweighting.
- [ ] `m_widom` — Widom test-particle insertion (excess chemical potential).
- [ ] `dg_ener` — free-energy from energy differences.

**2.4 — Clean up other discovered stubs**
- [ ] `visco` (placeholder formula, not Green-Kubo), `frameout` (hardcoded CA/ALA), `amber2gromos`
  (ignores input), `ener` (hardcoded LJ), `sasa_hasel`, `dssp` (no real H-bonding), `solute_entropy`.

### Priority 3 — py-gromos API & education
Design reference: `.local/polars` (Python API surface, docstrings, method chaining, pyo3 wrapping).
Focus: compositional API for **system/topology construction**.

> **Design target (FUTURE.md "Compositional topology"):** beat the gromos++ file-pipeline workflow
> (make_top → com_top → check_top → sim_box → ion, stringly-typed, file round-trips, late stdout
> errors). Model a topology as an **algebra**: `ForceField` (resolution context) → `BuildingBlock`
> (MTB template) → `Topology` (composable value with `+`/`*`/`solvate()`/`check()`/`from_sequence`).
> **Python expresses verbs; the Rust core owns every invariant** (exclusions, 1-4, LJ matrix, charge
> groups — the *same* functions the engine uses, never a second copy). Payoff: build → minimize →
> simulate in one address space, zero files; eager typed validation; reproducible recipes. The `*`
> operator (`water * 216`, `na * 3`) rides on the Dim 10 instancing model. A `LazyTopology`
> (Polars-style deferred build plan) is the second iteration. Refactor `gromos-tools` make_top/com_top
> into reusable library functions (not binary-only `main`s) so pyo3 can call them.

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

**Free win — mine gromosXX's own `check/*.t.cc` regression suite.** The gromosXX devs hard-code
per-term reference energies in `md++/src/check/`: `aladip.t.cc` carries `QuarticBond=18.053811`,
`NonBonded_newRF=-84.092443`, `DistanceRestraint=257.189539`, **and the perturbed/soft-core terms**
(`PerturbedQuarticBond`, `PerturbedSoftBond`, `PerturbedDihedral`, …); also `c16_cg.t.cc`
(coarse-grained), `lambdas.t.cc`, `scaling.t.cc`, `aladip_ls.t.cc` (lattice-sum). Porting these as
unit tests gives **per-term validation independent of running the md binary** — and is a genuine
*second source of truth* for exactly the perturbed terms Dim 11 says not to trust the C++ on (P1.5,
P1.7). High value, low cost; do it alongside the FEP/restraints work.

### Cross-cutting (do continuously) — differential audit (don't port the bugs) — FUTURE.md Dim 11
The reference suite is a bug oracle **only for wired paths.** Every *unwired* path (triclinic, FEP,
grid pairlist, Nosé-Hoover, EDS/GaMD, lattice-sum/PME) is un-audited surface where both inherited
gromosXX bugs *and* silent gromos-rs divergences hide. Rules, applied to every port:
1. **Reference test BEFORE wiring**, not after. (Would have caught the triclinic divergence, P1.4.)
2. **Grep the C++ for self-flagged defects before porting a subsystem:**
   `grep -rniE 'bug|fixme|wrong|hack' interaction/ algorithm/ math/` — the authors flag their own.
3. **Second-source uncertain/subtle physics** (RF self-terms, virial, Ewald): derive from the GROMOS
   book independently, implement to that, then diff against the C++; on disagreement, investigate —
   don't assume the C++ is right.
4. **Reproduce genuine GROMOS quirks as named, documented decisions** (like the "BONDANGLETYPE
   intentionally unsupported" entries), so a deliberate GROMOS-ism is never mistaken for a port bug.
5. Eventually a `--gromos-compat` vs `--corrected` split so "faithful" and "correct" can coexist.

Known items already logged: triclinic nearest-image divergence (P1.4, **live in our tree**), Martina
grid-pairlist bug (P1.6), perturbed RF self-term uncertainty (P1.5), thermostat flexible-constraints
hack (P1.3). Full table in FUTURE.md Dim 11.

### Parked / blocked
- `ext_ti_ana`, `nhoparam` — blocked behind Priority 2 (selection + fit + stats + references).
- Tutorials t_01–t_06 end-to-end — compute-limited; replaced near-term by minimal reference tests.

### Deferred breadth (tracked, not scheduled — low priority by design)
Subsystems present in gromosXX/gromos++ but intentionally low on our list, recorded so they're not
forgotten:
- **PME / lattice-sum electrostatics** — **RF stays the focus by design** (it's why GROMOS chose RF;
  PME was never well-optimized in gromosXX). NAMD/GROMACS optimize PME with *different logics*
  (FFT-based mesh, decomposition) — **needs a focused investigation before committing** whether to
  port gromosXX's lattice-sum or adopt a modern approach. Lower priority than RF-based physics.
  Ref: `interaction/nonbonded/interaction/latticesum.{h,cc}` (note the `// wrong!!!` traps — Dim 11).
- **Stochastic / Langevin dynamics (SD integrator)** — common but ~same priority tier as PME for our
  goals. `random_force` scaffolding exists; the SD leap-frog + friction/stochastic terms are unported.
  Ref: `algorithm/integration/` (stochastic), `math/random.*`.
- **Coarse-grained → Martini bridge** — *wanted* (Martini interop = ecosystem reach), but it's a
  research bet gated on a nonbonded-conventions investigation. Tracked as **FUTURE.md Dim 13**, not a
  near-term PLAN item. CG data model should ride the Dim 10 instancing refactor (beads = heavy repeats).
- **Polarisable / charge-on-spring force fields** — **explicitly out of scope** (low relevance for
  target use cases). Listed only to mark the decision.

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
