# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).

## [0.0.22] (2026-06-27)

### Features

- **gromos-core:** `PairlistAlgorithm` enum (`Standard`, `CellList`) with `from_imd()` dispatch
  and `update<BC>()` delegation via match — zero heap allocation on the hot path (9a-0)
- **gromos-core:** `from_imd()` slots match gromosXX exactly: 0=standard (forced), 1=grid
  (ExtendedGrid fallback→Standard), 2=grid_cell (CellList); IMD parser now recognises
  `"grid_cell"` keyword; no auto-heuristic — every value is an explicit instruction
- **gromos-core:** 9 unit tests covering full dispatch matrix (algorithm × box_type × n_atoms × chargegroups)
- **gromos-integrators, gromos-md, pyo3-gromos:** `StandardPairlistAlgorithm` field swapped to
  `PairlistAlgorithm` everywhere; `md.rs` calls `from_imd(imd.algorithm, …)` (9a-0)

### Testing

- **gromos-md:** `test_pairlist_margin.rs` — empirical CellList vs Standard margin test on
  `water_216_box` (100 steps); observed margin = 0.000e0 (bit-identical) (9a-1)
- **gromos-md:** `water_1000_spc_gridcell` reference — 1000 equilibrated SPC molecules
  (3.1057 nm box, vol7 tutorial), gromosXX Grid_Cell_Pairlist (7×7×7 grid); directly
  validates gromos-rs CellList against gromosXX; completes the proof chain

### Reference tests

- 38/38 pass (37 original + `water_1000_spc_gridcell`)

## [0.0.21] (2026-06-27)

### Features

- **gromos-core:** `gather.rs` — PBC molecule gathering: `gather_chain`, `gather_bond`,
  `gather_molecules`, `centre_of_geometry`, `centre_of_mass`; shared primitive for MD engine
  and analysis tools; 5 unit tests
- **gromos-core:** `fit.rs` moved to `gromos-analysis` (analysis-only); `gather.rs` stays in
  `gromos-core` (shared with MD engine)
- **gromos-analysis:** `fit.rs` — Kabsch rotational fit (Horn 1987 quaternion method);
  `superimpose()`, `rmsd()`, `kabsch_rotation()`, `weighted_cog()`; 7 unit tests
- **gromos-analysis:** `rmsd` binary — real Kabsch fit; @atomspec, @ref, @nofit, @pbc
- **gromos-analysis:** `nhoparam` — N-H bond order parameters S²; rotational fit, window
  averaging, `ee()`; port of GROMOS nhoparam algorithm
- **gromos-forces:** `energy.rs` — `single_point_energy()` + `EnergyParams`; computes
  bonded + nonbonded potential without running MD; used by `ener` binary
- **gromos-analysis:** `ener` — real energy recalculation from trajectory; per-frame
  E_bond, E_lj, E_crf, E_pot; optional PBC gathering
- **gromos-analysis:** `bar` — BAR iteration (numerically stable log-sum-exp, GROMOS
  reference, Shirts 2003); bootstrap error estimation
- **gromos-analysis:** `ext_ti_merge` — linear interpolation between λ windows;
  trapezoidal ΔG from merged curve
- **gromos-analysis:** `frameout` — full GROMOS feature parity: PBC gathering (@pbc),
  atom filtering (@include SOLUTE/SOLVENT/ALL), rotational fit (@ref @atomsfit),
  @spec ALL/EVERY/SPEC, @time range, cnf/pdb/trc output, @single; 8 integration tests
- **gromos-io:** `TrajectoryWriter::write_trc_frame()` — write POSITIONRED frames directly
  from positions (no Configuration needed); standard 3-column GROMOS format
- **gromos-io:** `write_pdb_positions()` in `pdb.rs` — write PDB from raw positions +
  optional topology for atom names; re-exported from lib.rs
- **gromos-io:** `TrajectoryReader` now handles both 3-column (standard GROMOS) and legacy
  7-column POSITIONRED formats; GENBOX made optional (vacuum trajectories)
- **gromos-core:** `selection.rs` — `m:` prefix: selects solute molecules by role
  (mirrors `s:` for solvent); improved error messages with full syntax cheatsheet and
  specific hints for common mistakes (`r:`, `mol:`, bare atom names as prefix)
- **gromos-analysis:** tests — 8 `frameout` integration tests + 7 trajectory round-trip
  unit tests in `gromos-io`

### Refactor

- **gromos-analysis:** `fit.rs` relocated from `gromos-core` → `gromos-analysis`;
  `gromos-forces` added as dependency of `gromos-analysis`
- All analysis binaries use only `gromos-io` for I/O — zero hand-rolled BufWriter/File
  in binary entry points
- PLAN.md restructured: Dim 10 section dissolved into Priority 2; priority order updated
  (CellList wiring → code quality → dihedral restraints → Dim 10 Phase 3)

### Reference tests

- 37/37 pass (unchanged); `frameout` integration tests added to `gromos-analysis`

## [0.0.20] (2026-06-26)

### Refactor

- **gromos-core (Dim 10 Phase 2a–3):** Dissolve the solute/solvent split — full instancing model.
  Direct field access everywhere: `topo.moltypes[0].bonds`, no wrapper methods.
  `MoleculeType` carries atoms + all bonded terms; `Atom = MolTypeAtom` (type alias);
  `Solute` reduced to ZST shell; `Solvent` struct and `Vec<Solvent>` field removed.
  `rebuild_flat_arrays()`: iac/mass/charge derived from instance registry (Phase 2a).
  `all_bonds_global()` / `all_angles_global()` / etc.: instance-iterating force-loop
  iterators with local→global index translation — one loop handles flexible solute,
  flexible solvent (future), and any repeated molecule type (Phase 3).
  `Topology::new()` pre-initializes `moltypes[0]` (SOLUTE) + `instances[0]` so no
  guard methods are needed at write sites.
  Follows GROMACS `gmx_moltype_t` / `gmx_molblock_t` convention exactly.
  37/37 reference tests pass byte-identical.

## [0.0.19] (2026-06-24)

### Features

- **gromos-forces:** perturbed bonded forces (P1.7 Step 2) — quartic bond, cos-harmonic
  angle, improper dihedral, proper dihedral; gromosXX-faithful λ-interpolation + dE/dλ;
  `calculate_perturbed_bonded_forces` wired into `Forcefield::apply` step 2b
- **gromos-forces:** perturbed nonbonded corrections (P1.7 Step 3) — soft-core LJ+CRF
  dual-topology via `perturbed_pairlist_correction`; RF self-energy, excluded-pair, 1-4,
  and PERTATOMPAIR corrections; `ch4_water_fep` (CH4→dummy in 999 SPC water, λ=0.5)
  passes to <1e-6 kJ/mol vs gromosXX
- **gromos-io:** fix PERTURBATION block multi-line parser (ALPHLJ/ALPHC silently dropped);
  fix effective alpha = per-atom × global (matches `in_perturbation.cc:1308`)
- **gromos-io:** `FreeEnergyWriter` + `read_free_energy_trajectory` (FREEENERGY03 blocks);
  md binary writes `@trg` when NTG≠0 at NTWE frequency (P1.7 Step 4)
- **gromos-core:** `Stat` — `ave()`, `rmsd()`, `ee()` block averaging, port of gromos++
  `gmath/Stat` (Allen–Tildesley); 4 unit tests
- **gromos-analysis:** `ext_ti_ana` — real TI analysis; reads N `.trg` files, ⟨dH/dλ⟩ ±
  ee() per window, trapezoidal ΔG integration (P1.7 Step 5 + P2.2 stat + P2.3 program)
- **gromos-core:** full gromos++ AtomSpecifier grammar — `a:NAME`, `1:name,name`,
  `1:res(nr:atom)`, `1:res(name:atom)`, `not(spec)`, `minus(spec)`, `;`-union, `all`/`no`;
  routes via `topology.molecules` (no `num_solute_atoms()` threshold, Dim 10 ready);
  30 tests, every index confirmed against gromos++ `atominfo` on aladip topology
- **gromos-core:** per-atom metadata accessors — `atom_name(i)`, `residue_nr(i)`,
  `residue_name(i)`, `molecule_nr(i)` covering solute + solvent uniformly
- **gromos-analysis:** `atominfo` binary — gromos++-compatible TITLE+ATOMS output
- **gromos-core:** Dimension 10 Phase 1 — instancing model alongside legacy structs;
  `Role`, `MolTypeAtom`, `MoleculeType`, `MoleculeInstance`; `Topology::moltypes` +
  `instances` populated by `init_solute_moltype()` and `solvate()`; per-atom accessors
  prefer moltype path; `s:` in AtomSelection uses `role == Role::Solvent`; `promote()`;
  37/37 reference tests byte-identical

### Test infrastructure

- Shared topology/coord/ptp files moved to `tests/gromosXX_references/shared/`;
  all `input.toml` paths updated from `../../` to `../shared/`
- `.trg` free-energy trajectory tracked in `ch4_water_fep` (dH/dλ at 1e-6 rel tol)
- `run_references.py` removed (both copies); Rust runner is the only harness
- 37 tests pass (up from 36); 1 ignored

## [0.0.18] (2026-06-19)

### Features

- **gromos-core:** add `PerturbedAtom` and `PerturbedAtomPair` structs; extend
  `PerturbedSolute` with `atoms` + `atom_pairs` fields; add `is_perturbed: Vec<bool>`
  to `Topology` (P1.7 Step 1 scaffolding)
- **gromos-io:** replace ptp.rs stub with a full gromosXX-faithful `.pttopo` reader
  (`read_pttopo`); parses PERTATOMPARAM, PERTATOMPAIR, PERTBONDSTRETCH(H),
  PERTBONDANGLE(H), PERTIMPROPERDIH(H), PERTPROPERDIH(H); validated against
  `aladip.pttopo` (3 atoms, 1 pair, 2 bonds, 2 angles, 2 impropers, 2 dihedrals)
- **gromos-md:** wire `@pttopo` into md binary — when NTG≠0, read perturbation topology
  and populate `topo.perturbed_solute` + `topo.is_perturbed`
- **gromos-core:** fix `Energies::total()` to include `distanceres_total` (matches
  gromosXX: `total = potential_total + kinetic_total + special_total` where
  `special_total` includes distance restraints)
- **gromos-md:** add `nacl_1water_distres` reference test — instantaneous distance
  restraint (NTDIR=2) on Na-Cl pair; validates total energy accounting and step-0
  forces; 36 of 36 tests pass

### Reference test matrix

- 36 tests pass (up from 34); new: `nacl_1water_distres` (distance restraint)

## [0.0.17] (2026-06-17)

### Features

- **gromos-forces:** port gromosXX distance restraints (PLAN.md §1.6, first priority)
  - `DistanceRestraint` / `DistanceRestraints` — faithful translation of
    `distance_restraint_interaction.cc`: RAH dimensionality encoding
    (`dim_base ∈ {0,10,20,30,40,50,60}` × form `∈ {-1,0,+1}`), harmonic/linear
    switching (gromosXX `r_linear = DIR0 = 0.3 nm`), mode 1 / mode 2 (NTDIR=2 ×w0)
  - `PerturbedDistanceRestraint` / `PerturbedDistanceRestraints` — translation of
    `perturbed_distance_restraint_interaction.cc`: λ-interpolated r0/w0,
    hidden-restraint prefactor `2^(n+m)·λⁿ·(1-λ)ᵐ`
  - Reference test `test_distance_restraint_gromosxx_reference` validates both at
    gromosXX's own hard-coded values: DistanceRestraint = 257.189539 kJ/mol,
    PerturbedDistanceRestraint = 195.899012 kJ/mol (aladip.distrest, tol 1e-3)
- **gromos-core:** add `DistanceRestraintSpec` / `PerturbedDistanceRestraintSpec` to
  `topology.rs`; extend `Topology` with `distance_restraints` /
  `perturbed_distance_restraints` Vec fields
- **gromos-core:** add `distanceres_total: f64` to `Energy` struct
- **gromos-io:** new `distanceres.rs` parser for DISTANCERESSPEC / PERTDISRESSPEC
  blocks (virtual atom type 0 only; unsupported types logged + skipped)
- **gromos-io:** parse DISTANCERES block (`NTDIR NTDIRA CDIR DIR0 TAUDIR FORCESCALE
  VDIR NTWDIR`) and PERTURBATION block (`NTG NRDGL RLAM DLAMT ALPHLJ ALPHC NLAM
  NSCALE`) in `imd.rs`; add `ImdParameters::lambda_and_derivative()` helper
- **gromos-integrators:** wire `DistanceRestraints` / `PerturbedDistanceRestraints`
  into `Forcefield` (step 4c after position restraints); result accumulated into
  `conf.current_mut().energies.distanceres_total`

### Reference test matrix

- 13 restraint tests pass (13 restraint-specific, up from 0 gromosXX-faithful)
- All other workspace tests unchanged

## [0.0.16] (2026-06-14)

### Features

- **gromos-core:** add `CellListPairlistAlgorithm` (FUTURE.md Dim 9a) — a charge-group-aware
  linked-cell pairlist, drop-in via the existing `update<BC>()` interface
  - Bins each chargegroup's reference position (solute CG → COG, solvent CG → first-atom
    position, matching `StandardPairlistAlgorithm`'s own distance conventions) into a grid that
    exactly tiles the box, with cells ≥ `long_range_cutoff + skin` and a periodic-wrapped,
    deduplicated 27-cell neighbor search (correct even when `grid_dim` is 1 or 2 along an axis)
  - O(N) for `BoxType::Rectangular`; vacuum/triclinic/truncated-octahedron boxes fall back to
    `StandardPairlistAlgorithm`'s O(N²) path (correct, not yet accelerated)
  - Pairs are classified by `min/max(cg_a, cg_b)` against `num_solute_chargegroups` — both
    chargegroups' roles — avoiding the Martina solute/solvent misclassification bug
    (`extended_grid_pairlist_algorithm.cc:1309`)
  - Replaces the old, unexported, non-chargegroup-aware `GridCellPairlistAlgorithm` stub
- **gromos-core:** add set-equality tests (`test_cell_list_matches_standard_*`) validating that
  `CellListPairlistAlgorithm` produces identical `solute_short/long`/`solvent_short/long` pair
  sets to `StandardPairlistAlgorithm`, covering periodic wrap, the single-cell fallback, and
  mixed solute+solvent with exclusions

### Reference test matrix

- 35 of 35 tests pass (unchanged — `StandardPairlistAlgorithm` remains the only algorithm wired
  into the `md`/`pyo3-gromos` binaries)

## [0.0.15] (2026-06-14)

### Features

- **gromos-core:** port `truncoct_triclinic_rotmat` / `truncoct_triclinic_box` / `truncoct_triclinic`
  (gromosXX `math/transformation.cc`) for NTB=-1 truncated-octahedron boxes
  - `forward`/`!forward` APIs convert the legacy cube-edge BOX block to/from the lower-triangular
    triclinic box vectors and rotate positions/velocities between the cube and triclinic frames
  - Add `Box::truncated_octahedral()`
- **gromos-core:** replace `Triclinic::nearest_image`'s textbook `frac - frac.round()` reduction
  with gromosXX's iterative while-loop z→y→x reduction (`boundary_implementation.cc:285-318`) —
  not equivalent for strongly triclinic cells (FUTURE.md Dim 11, finding #1)
- **gromos-integrators:** `forcefield.rs` periodicity refresh now also matches
  `BoxType::TruncatedOctahedral`
- **gromos-md:** wire NTB=-1 in the `md` binary — convert the BOX block, rotate positions/velocities
  on read, build `Periodicity::Triclinic` from the resulting box, and rotate
  `FREEFORCERED`/`CONSFORCERED` back to the cube frame on output
  (`out_configuration.cc::_print_forcered`)
- **gromos-md:** add `aladip_trunc_oct` gromosXX reference test (truncated octahedron, NTB=-1)

### Bug fixes

- **gromos-core:** fix `truncoct_triclinic_rotmat`'s forward/backward rotation — gromosXX's
  `product(rot, v)` (`gmath.h`) computes `rot^T * v`, not `rot * v`, because `GenericMatrix`'s
  3-`Vec` constructor stores its arguments as rows while `product()` contracts over the matrix's
  first index. The initial port had the two branches exactly swapped. Added a regression test
  pinned to a gromosXX debug-build reference value for `aladip.conf` atom 21.

### Reference test matrix

- 35 of 35 tests pass (up from 34); new: `aladip_trunc_oct` (truncated octahedron / NTB=-1)

## [0.0.14] (2026-06-12)

### Features

- **gromos-integrators:** implement Nosé-Hoover thermostat — single NHC and NHC chain
  - `NoseHooverThermostat` struct with `new_single_bath` (nhc=1) and `new_chain_bath` (nhc≥2)
  - Single NHC: `ζ[0] += dt/τ² · (T_free/T₀ - 1)`, `scale = 1 - ζ[0]·dt`
  - Chain NHC: sequential tail-to-head ζ update (gromosXX `calc_chain_scaling`)
  - Wired in `md.rs` via `imd.temp_bath[0].algorithm` dispatch (0=Berendsen, 1=NHC, N=chain)
- **gromos-io:** fix MULTIBATH algorithm code mapping to match gromosXX exactly
  (0=weak-coupling, 1=nose-hoover, N≥2=nose-hoover-chains); add `nhc_chain` field
- **gromos-md:** add `water_216_nvt_nosehoover` and `water_216_nvt_nhc_chain` reference tests

### Reference test matrix

- 34 of 34 tests pass (up from 32); new: NHC single, NHC chain

## [0.0.13] (2026-06-12)

### Features

- **gromos-integrators:** wire SETTLE, LINCS, and COM rotation removal into the MD algorithm
  sequence; all three are now reference-tested and passing
  - `SettleAlgorithm` / `LincsAlgorithm` wired for solvent (NTCS=settle/lincs) and solute
    (NTCP=lincs); reference tests: `nacl_1water_settle`, `nacl_1water_lincs`, `aladip_vacuum_lincs`
  - COM rotation removal (`RemoveCOMMotion`) fully wired; reference test: `water_216_box_com_rot`

### Bug fixes

- **gromos-integrators:** fix COM rotation removal for PBC systems
  - Use minimum-image (put_into_box) wrapped positions when computing angular momentum and
    inertia tensor, mirroring gromosXX's `gather_chargegroups` init convention
    (`periodicity.cc:175`); without this the position-COM and angular momentum diverged
  - Suppress *periodic* COM rotation removal when the box is periodic, matching gromosXX
    `configuration.cc:555-560` (`param.centreofmass.remove_rot = false` for non-vacuum PBC)
  - Remove leftover debug `eprintln!` that printed all 648 atom positions at step 0
- **gromos-io:** fix `imd.rs` parsing for NTCP/NTCS lincs/settle keywords and LINCS order params

### Reference test matrix

- 32 of 32 tests pass (up from 28); new: SETTLE, LINCS (solvent+solute), COM rotation removal

## [0.0.12] (2026-06-07)

### Features

- **gromos-core:** add `random` module — GSL-compatible `mt19937` RNG (`GslMt19937`) and
  `gsl_ran_gaussian` (Marsaglia polar method), plus `generate_velocities` for NTIVEL=1
  Maxwell-Boltzmann initial velocity generation, matching gromosXX `util::generate_velocities` /
  `math::RandomGeneratorGSL` bit-for-bit
- **gromos-md:** wire NTIVEL=1 into `md` binary's velocity setup (reads NTIVEL/IG/TEMPI from
  INITIALISE block); add `water_single_genvel` reference test verifying generated velocities
  reproduce gromosXX trajectories exactly

## [0.0.11] (2026-06-07)

### Features

- **gromos-integrators:** add steepest descent energy minimization algorithm
- **gromos-md:** add energy minimization reference tests (vacuum + solvated, with/without SHAKE)
- **gromos-io:** add IFP, MTB, posres, mk_script, jobs, and script_template parsers (~2,600 lines)
- **gromos-tools:** add `prep_posres` utility; rewrite `ion`, `sim_box`, `pdb2g96`, `com_top`, `make_top`, and `mk_script` against the new parsers
- **gromos-tools:** implement `build_box` (replaces hardcoded placeholder) — generates a condensed-phase system on a regular grid sized to a target density, matching gromosXX `build_box.cc`
- **core/forces:** update configuration and position restraints handling

### Performance

- **forces/shake:** optimize nonbonded hot path and SHAKE with cache-friendly data structures
  - Replace `HashSet` exclusions with sorted `Vec` + binary search in topology
  - Flatten LJ parameter matrix into a contiguous `LJParamMatrix` (removes double indirection)
  - Add `ShakeBuffers` for precomputed constraint lists and reusable skip arrays
  - Add CG-grouped nonbonded kernel (`CGPairGroup`): compute PBC nearest-image once per charge-group pair
  - Add parallel (`_parallel`) and no-virial (`_novirial`) innerloop variants and `ForceStorage::merge` for thread-local reduction
  - Add Criterion benchmark suite (`gromos-md/benches/md_bench.rs`)
  - Fix duplicate energy/force accumulation in `calculate_dihedral_new_forces`

### Refactor

- **workspace:** remove the `gromos-cli` crate (superseded by the per-crate binary layout from the 0.0.10 restructure)
- **workspace:** condense and reorganize `PLAN.md` into a priority-ordered roadmap (MD engine physics → analysis foundations → py-gromos/education → code quality)

### Bug Fixes

- **ci:** fix Rust workflow steps using an incorrect `working-directory: gromos-rs` (the workspace lives at the repo root); fix coverage lcov output path and the crates.io publish job directory (`crates/gromos`)

### Chores

- **workspace:** clean up workspace `Cargo.toml`, `Makefile`, and update `py-gromos` reference test paths

## [0.0.10] (2026-05-17)

### Refactor

- **workspace:** restructure crates into focused responsibilities
  - `gromos-md`: 8 simulation engine binaries (md, md_mpi, md_mpi_cuda, mdf, remd, repex_mpi, eds, gamd) + integration tests
  - `gromos-tools`: 30 system construction binaries organized in subdirectories (topology/, box/, conversion/, utilities/)
  - `gromos-analysis`: 66 analysis binaries organized in subdirectories (structural/, energy/, distribution/, dynamics/, free_energy/, trajectory/, noe/, clustering/, xray/, special/) + existing library code
  - `gromos-cli`: slimmed to thin unified `gromos` multicall binary (clap only)
  - `gromos` facade: removed `gromos-analysis` re-export to avoid circular dependency; analysis bins use `gromos_core`/`gromos_io` directly
  - All 21 reference tests pass, workspace compiles clean

## [0.0.9](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.8...v0.0.9) (2026-03-29)

### Bug Fixes

- **forces:** fix CRF energy mismatch (~3.82 kJ/mol) for systems with solvent near cutoff boundary
  - Expand solvent long-range pairlist to all atom pairs (was: first-atom only with shared PBC shift)
  - Remove HEAVISIDE truncation from `lj_crf_innerloop` to match gromosXX `#undef XXHEAVISIDE` default
  - Long-range solvent now processed with `lj_crf_innerloop` (per-atom nearest_image) instead of `solvent_innerloop`

### Features

- **tests:** add water_10_box and nacl_water_box_shifted reference systems
  - water_10_box: 2 ions + 10 SPC waters (32 atoms), positions away from cutoff boundaries
  - nacl_water_box_shifted: nacl_water_box with perturbed positions near cutoff
- **tests:** promote nacl_water_box, nacl_water_box_shifted, nacl_3water_cutoff from ignored to active
  - 14 of 19 reference tests now pass (was 11)
- **io:** add verbose logging and md_output.log capture to run_references.py
- **io:** improve force trajectory writer robustness

## [0.0.8](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.7...v0.0.8) (2026-03-08)

### Features

- **core:** add solute/solvent boundary fields to Topology
  - Add `num_solute_chargegroups` field and `num_solute_atoms()` method
  - Set during `build_topology()` for solute/solvent dispatch

- **forces:** gromosXX-compatible nonbonded force architecture
  - Add HEAVISIDE truncation to `lj_crf_innerloop` (skip pairs beyond cutoff²)
  - Add `solvent_innerloop` with shared PBC shift from O-O nearest_image
  - Split `rf_excluded_interactions` into solute/solvent paths:
    - Solute: self-term + excluded pair forces + full RF energy
    - Solvent: no self-term, no forces, only distance-dependent energy
  - Add `cutoff_sq` to CRFParameters

### Refactor

- **pairlist:** separate solute/solvent CG pairs with correct distance metrics
  - Solute CGs: center-of-geometry distance, exclusion checks
  - Solvent CGs: first-atom position distance (gromosXX convention)
  - Intra-CG non-excluded pairs for solute CGs
  - Solvent-solvent stores first-atom pairs only (expanded in innerloop)
  - Add debug logging for CG positions, distances, pair classification

## [0.0.7](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.6...v0.0.7) (2026-03-02)

### Refactor

- **io:** gromosXX-compatible topology init order with solvate() ([60207d9](https://github.com/yerkoescalona/gromos-rs/commit/60207d9584e67c0ae622ed6267f19fd38965467e))
  - Follow gromosXX initialization: read_parameter → read_topology → solvate → read_configuration
  - Add Topology::solvate(nsm) method matching gromosXX topo.solvate(0, nsm)
  - Add SolventAtomTemplate/SolventConstraintTemplate types on Topology
  - build_topology() now stores solvent template without expanding
  - Remove build_topology_with_solvent() — NSM comes from IMD SYSTEM block
  - Reorder md.rs: IMD → topology+solvate → coordinates (was: IMD → coords → topo)
  - Store chargegroup_codes on Topology for later use
  - LJ matrix now includes solvent IAC types at build time
  - All 6 reference tests pass (pair_lj, pair_lj_mixed, nacl_pair, water_single, benzene_vacuum, water_3_box)

## [0.0.6](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.5...v0.0.6) (2026-02-27)

### Features

- **forces:** gromosXX-compatible bonded forces, topology parsing, and SHAKE ([5d409e7](https://github.com/yerkoescalona/gromos-rs/commit/5d409e7da1ffb9528db09d122e68515bda4df6bc))
  - Fix bonded force vector conventions to match gromosXX (v = pos(i) - pos(j))
  - Simplify quartic bond force: avoid unnecessary r division, use r^2 directly
  - Add NTF flag support to calculate_bonded_forces_ntf() for selective force terms
  - Add debug logging for individual bonded energy components
  - RF excluded corrections: add forces (not just energy), rename to rf_excluded_interactions
  - Parse dihedral/improper types and terms from topology (TORSDIHEDRALTYPE, DIHEDRAL, IMPDIHEDRAL)
  - Fix SOLUTEATOM parser: exclusions on same line (col 8+), skip INE14 line
  - Rewrite SHAKE to match gromosXX convergence (tolerance on |constraint_length - r|)
  - Add shake_algorithm.rs module
  - Add ForceStorage to nonbonded for combined force+energy tracking

## [0.0.5](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.4...v0.0.5) (2026-02-27)

### Refactor

- **io:** gromosXX-compatible CLI and IO cleanup ([c0e9ea2](https://github.com/yerkoescalona/gromos-rs/commit/c0e9ea2603942bed73944dab073cd58045966f5d))
  - Move coordinate reader from md.rs to gromos-io (read_coordinates + CoordinateData)
  - Add GENBOX, POSITIONRED, VELOCITYRED block support in coordinate reader
  - Rewrite IMD parser for gromosXX positional format (comment headers + data lines)
  - Add NTF force flags, COMTRANSROT, shake_tol to ImdParameters
  - Refactor md.rs CLI to match gromosXX @-parameters: @topo @conf @input (required), @fin @trc @tre @trf @trv (output), @pttopo @posresspec @refpos @distrest @angrest @dihrest @gamd (input), @verb @print @version @develop (control)
  - All simulation parameters now come from @input .imd/.in file
  - Fix pre-existing f64/usize type errors in energy_binary and trajectory_binary tests
  - 9 new/updated tests for coordinate and IMD parsing
  - pair_lj reference still passes at ~1e-11 precision

## [0.0.4](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.3...v0.0.4) (2026-02-24)

### Refactor

- **core:** implement gromosXX Algorithm/AlgorithmSequence pattern ([c61cd57](https://github.com/yerkoescalona/gromos-rs/commit/c61cd57c5b7b1eb29c0f067fd02038560b923ab9))
  - Add Algorithm trait and AlgorithmSequence in gromos-core with 7 unit tests
  - Implement algorithm wrappers: Forcefield, LeapFrogVelocity, LeapFrogPosition, TemperatureCalculation, EnergyCalculation
  - TemperatureCalculation uses gromosXX averaged kinetic energy formula
  - Refactor md.rs to use AlgorithmSequence instead of manual force+integrate loop
  - Energy convention follows gromosXX: write to current(), exchange_state moves to old(), read from old() for output
  - Energy output uses full f64 scientific notation for exact comparison
  - All 10 steps of pair_lj match gromosXX reference to ~1e-11 precision
  - Fix minor clippy warnings in gamd, remd, integrator, pairlist

## [0.0.3](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.2...v0.0.3) (2026-02-22)

### Refactor

- **forces:** enhance nonbonded interactions with RF corrections and update math structures ([51b8681](https://github.com/yerkoescalona/gromos-rs/commit/51b86819e74ffcbe918b30b00f180b4c5be2c9ae))
- **core:** support f64 instead of f32 ([ade5423](https://github.com/yerkoescalona/gromos-rs/commit/ade5423ca49fc8974dfdc624f71770a4b7c10a92))

### Chores

- remove dead logging module and replace black with ruff ([de003df](https://github.com/yerkoescalona/gromos-rs/commit/de003dfef3066db053255dbc58087739becf87f7))

## [0.0.2](https://github.com/yerkoescalona/gromos-rs/compare/v0.0.1...v0.0.2) (2025-12-21)

### Refactor

- **workspace:** reorganize into crate-based workspace structure ([7d97a54](https://github.com/yerkoescalona/gromos-rs/commit/7d97a54096d902529e878374395b837717e9343f))

## [0.0.1](https://github.com/yerkoescalona/gromos-rs/commits/v0.0.1) (2025-11-22)

### Refactor

- update force writer integration tests for new API ([ec11525](https://github.com/yerkoescalona/gromos-rs/commit/ec11525f4e08aabc8d621a187cd8875afd3b90f6))

### Chores

- clean and update README files ([27f80bc](https://github.com/yerkoescalona/gromos-rs/commit/27f80bc3abcad0df3d19bd20d5a5291ce7f5fb8e))
- add GROMOS and claude entries to .gitignore ([926ff5c](https://github.com/yerkoescalona/gromos-rs/commit/926ff5c41d8658b964a09aa9d93682f7627ba264))

### Initial

- initial commit ([c0a551c](https://github.com/yerkoescalona/gromos-rs/commit/c0a551cfcac12b0a60daaad908d151f689395732))
