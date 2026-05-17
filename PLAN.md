# gromos-rs — Status & Plan

we will focus on `md`, so therefore focus on `cargo build --release --bin md`

letsgo reviewing things with compare references
modify `PLAN.md` in case you advance
in case we commit update CHANGELOG.md and update the Cargo.toml version

DONT MODIFY THE REFERENCES in gromosXX_references/*/expected/ — those are gromosXX ground truth.
To regenerate references: `python3 gromosXX_references/run_references.py --md-binary /path/to/gromosXX/md`
New reference systems should be added via run_references.py (add to SYSTEMS list).

Validation is done via Rust integration tests (replaces the old compare_reference.py):
  `cargo test -p gromos-md --test test_gromosXX_references`
New test systems: add a `ref_test!(name, "system_dir")` line in:
  `crates/gromos-md/tests/test_gromosXX_references.rs`

Here you can find the gromosXX code: .local/gromosXX/md++/src
Here you can find realistic tutorials: .local/gromos_tutorial_livecoms/tutorial_files
Here you can find theory behind everything: .local/doc/gromos_book
Here you can find forcefields: .local/gromosXX/forcefields


## Architecture

```
gromos-core        → types, math, Vec3, BoundaryCondition trait, Topology, Configuration, Algorithm trait
gromos-forces      → bonded (quartic/harmonic bonds, angles, dihedrals, impropers, cross-dihedrals)
                     nonbonded (LJ + CRF, rf_excluded_interactions, pairlist)
                     PME, QM/MM, restraints, polarization (stubs/partial)
gromos-integrators → LeapFrog, VelocityVerlet, StochasticDynamics, SteepestDescent
                     SHAKE, SETTLE, LINCS constraints
                     Berendsen thermostat/barostat, Nosé-Hoover, Andersen, Parrinello-Rahman
                     EDS, GaMD, REMD, FEP
gromos-io          → topology parser, coordinate reader (GENBOX, POSITION, VELOCITY),
                     IMD parser, trajectory/energy/force writers
gromos-md          → simulation engines (md, md_mpi, remd, eds, gamd, mdf)
gromos-tools       → system construction & preparation (make_top, sim_box, pdb2g96, mk_script, ...)
gromos-analysis    → trajectory analysis & post-processing (rmsd, rdf, ene_ana, hbond, ...)
pyo3-gromos        → Python bindings
py-gromos          → Python package (separate build)
```

**Crate responsibilities:**
- **gromos-md** — simulation engines only. Depends on the full physics stack (forces, integrators, constraints, io).
- **gromos-tools** — system construction programs (topology building, solvation, ion placement, coordinate conversion, script generation). Depends mainly on gromos-io + gromos-core.
- **gromos-analysis** — trajectory analysis and post-processing. Depends on gromos-io + gromos-core + analysis library.

**CRITICAL DESIGN RULE — no I/O duplication:**

All file format parsing and writing lives ONLY in `gromos-io`. No binary in gromos-md, gromos-tools, or
gromos-analysis re-implements reading/writing of topology, coordinates, trajectories, energy files, or IMD
parameters. This is the #1 lesson from gromosXX/gromos++: every program had its own bespoke parser for
the same file formats, leading to subtle incompatibilities, duplicated bugs, and unmaintainable code.

Concretely:
- **Topology (.top)** — `gromos_io::topology::parse_topology()` is the ONE parser. All crates use it.
- **Coordinates (.cnf/.g96)** — `gromos_io::coordinate::read_configuration()` is the ONE reader.
- **Trajectories (.trc/.tre/.trf/.trv)** — `gromos_io::trajectory` readers/writers used by everyone.
- **IMD parameters (.imd/.in)** — `gromos_io::imd::parse_imd()` is the ONE parser.
- **PDB (.pdb)** — `gromos_io::pdb` is the ONE parser/writer.
- **Building blocks (.mtb/.ifp)** — will live in `gromos_io::building_blocks` (not in gromos-tools bins).
- **Energy trajectory** — `gromos_io::energy` reader used by both md (writing) and analysis (reading).

If a binary needs to read/write a format, it calls `gromos-io`. If the format isn't supported yet,
the parser gets added to `gromos-io`, NOT locally in the binary. Zero exceptions.

## MD Loop (AlgorithmSequence)

```
RemoveCOMMotion (if NTICOM/NSCM) → Forcefield → LeapFrogVelocity → BerendsenThermostat (if TAU>0) → LeapFrogPosition → SHAKE (if ntc>1) → TemperatureCalculation → EnergyCalculation
```

- Forcefield: pairlist update → bonded (NTF-controlled) → nonbonded (lj_crf_innerloop) → rf_excluded_interactions
- Energy convention: write to current(), exchange_state moves to old(), read old() for output
- Kinetic energy: averaged between old/new velocities (gromosXX formula)

## Decisions Taken

- f64 everywhere (not f32)
- gromosXX `@` CLI convention: `@topo @conf @input @fin @trc @tre @trf @trv @verb ...`
- All simulation parameters from `@input` .imd/.in file (not CLI flags)
- Bonded force vectors: `v = pos(i) - pos(j)` (gromosXX convention)
- Boundary condition from box_dims: vacuum if (0,0,0), rectangular otherwise
- GENBOX block parsed (box_type + dims) but box_type ignored — only dims used
- Energy output: full f64 scientific notation for exact comparison
- Tolerances: force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9

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

**21 of 21 tests pass.** All levels fully passing.

## What Works

- LJ + CRF nonbonded forces (vacuum and PBC code paths)\n- 1-4 interactions: cs6/cs12 LJ + scaled CRF, parsed from LJPARAMETERS (CS12/CS6) and SOLUTEATOM (INE14)
- All bonded force types: quartic/harmonic bonds, cos-harmonic/harmonic angles, dihedrals, impropers, cross-dihedrals
- NTF flag control for selective bonded force evaluation
- RF excluded interactions (forces + energy + self-terms)
- Pairlist: StandardPairlistAlgorithm with configurable update frequency
  - Chargegroup-based mode: CG center-of-geometry distance check, atom pairs between CGs
  - Atom-based mode: direct atom-atom distance check
  - Twin-range: short (RCUTP) and long (RCUTL) pairlists with force caching
- Boundary conditions: Vacuum, Rectangular (minimum image), Triclinic (defined but not wired)
- SHAKE constraints: solute bonds (NTC>1) + solvent constraints (NTCS>0)
  - NtcMode enum: SolventOnly (NTC=1), HydrogenBonds (NTC=2), AllBonds (NTC=3)
  - Constraint force accumulation into conf.old().constraint_force
  - Virial tensor contribution: ref_r ⊗ ref_r · λ/dt²
  - skip_now/skip_next optimization for converged constraints
  - shake_positions() and shake_velocities() for initial configuration (NTISHK)
  - Iterates SOLVENTCONSTR template over all solvent molecules
  - Velocity correction applied to all constrained atoms
- COM motion removal: NTICOM (initial) + NSCM (periodic every N steps)
  - Placed first in algorithm sequence (before Forcefield), matching gromosXX convention
  - Translational COM velocity subtracted from all atom velocities
- SETTLE, LINCS constraints (implemented but not wired)
- Berendsen thermostat & barostat (parsed from MULTIBATH/PRESSURESCALE IMD blocks)
- Pressure / virial architecture matching gromosXX:
  - prepare_virial (KE tensor) + atomic_to_molecular_virial in Forcefield
  - PressureCalculation algorithm: P = (KE + 0.5*virial) * 2/V
  - BerendsenBarostat algorithm: isotropic box+position scaling
  - Pressure groups parsed from PRESSUREGROUPS topology block
- Topology parser: all block types including SOLUTEATOM exclusions, SOLVENTCONSTR, LJPARAMETERS, PRESSUREGROUPS
- Solvent expansion: SOLVENTATOM/SOLVENTCONSTR parsed, NSM auto-computed from coord count
  - build_topology(parsed, n_coord_atoms) — no manual NSM parameter needed
  - Solvent atoms expanded into flat arrays (mass, charge, iac)
  - Chargegroups built from CGC codes + one CG per solvent molecule
  - Intra-molecular solvent exclusions generated
- Coordinate reader: POSITION, POSITIONRED, VELOCITY, VELOCITYRED, GENBOX
- Energy/trajectory/force writers (ENERTRJ, POSITIONRED, FREEFORCERED/CONSFORCERED blocks)
- Force trajectory output (@trf): writes FREEFORCERED + CONSFORCERED matching gromosXX format
  - Wired in md binary: `@trf <file>` produces per-step force output
  - Tests compare forces atom-by-atom with tolerance 1e-6 kJ/(mol·nm)

## TODO

### PRIORITY — System Construction Programs (tutorial-ready)
Goal: make the molecular system construction pipeline work end-to-end so users can follow
the realistic tutorials at `.local/gromos_tutorial_livecoms/tutorial_files/`.

Programs needed for tutorials t_01 through t_06 (in pipeline order):

#### System building (critical path)
- [ ] **make_top** — PARTIAL, needs real MTB/IFP building block parsing
  - Currently generates a hardcoded single-type topology (N identical atoms)
  - Must: read MTB files (MTBUILDBLSOLUTE/MTBUILDBLSOLVENT), read IFP (54A7/54A8 forcefield)
  - Must: assemble topology from building blocks + linking rules
  - Tutorials: t_01 (GB3 protein), t_02 (ligands), t_06 (methanol)
  - gromosXX++ ref: `.local/gromosXX/md++/` has no make_top (it's a gromos++ tool)
  - Forcefield files: `.local/gromosXX/forcefields/` (54A7, 54A8)
- [ ] **build_box** — PARTIAL, hardcodes water/3-atom; needs real solute file reading
  - Must: read solute coordinates, replicate to fill box at target density
  - Tutorial: used for creating pure solvent boxes
- [ ] **com_top** — REAL ✓ (merges multiple topologies)
- [ ] **check_top** — REAL ✓ (validates topology)
- [ ] **pdb2g96** — REAL ✓ (PDB → GROMOS coordinate format)
- [ ] **sim_box** — REAL ✓ (solvation: place solute in solvent box, remove clashes)
- [ ] **ion** — REAL ✓ (place counter-ions — verify it does placement not just analysis)
- [ ] **mk_script** — REAL ✓ (generate IMD input / run scripts)
- [ ] **make_pt_top** — REAL ✓ (perturbation topology for TI/FEP)

#### Analysis programs (needed by tutorials, already real)
- [ ] **ene_ana** — REAL ✓ (energy trajectory statistics)
- [ ] **rmsd** — REAL ✓ (RMSD calculation)
- [ ] **rmsf** — REAL ✓ (per-atom fluctuations)
- [ ] **rdf** — REAL ✓ (radial distribution function)
- [ ] **hbond** — REAL ✓ (hydrogen bond analysis)
- [ ] **frameout** — REAL ✓ (extract frames from trajectory)
- [ ] **tser** — REAL ✓ (time series extraction)
- [ ] **trs_ana** — REAL ✓ (trajectory statistics)
- [ ] **bar** — REAL ✓ (Bennett Acceptance Ratio for FEP)
- [ ] **ext_ti_ana** — PARTIAL, falls back to placeholder if no input
- [ ] **reweight** — REAL ✓ (GaMD reweighting)
- [ ] **filter** — REAL ✓ (trajectory frame filtering)
- [ ] **nhoparam** — verify implementation status (NHO order parameters, t_01)

#### Validation plan
- [ ] Run tutorial t_01 (protein MD) end-to-end with gromos-rs binaries
  - Pipeline: make_top → com_top → check_top → pdb2g96 → sim_box → ion → md → ene_ana/rmsd/rmsf
- [ ] Compare output of each step against gromosXX++ reference output
- [ ] Document which tutorial steps work and which are blocked

#### Notes
- Many analysis bins in gromos-cli duplicate what should live in gromos-analysis;
  keep them here for now to ensure methodologies match, refactor later.
- Focus order: (1) make_top (blocks everything), (2) build_box, (3) validate pipeline

### Binary Placement — gromos-md / gromos-tools / gromos-analysis

**gromos-md** (8 bins) — simulation engines:

| Programs |
|----------|
| `md`, `md_mpi`, `md_mpi_cuda`, `mdf`, `remd`, `repex_mpi`, `eds`, `gamd` |

**gromos-tools** (27 bins) — system construction & preparation:

| Category | Programs |
|----------|----------|
| Topology construction | `make_top`, `com_top`, `con_top`, `link_top`, `red_top`, `check_top`, `addvirt_top`, `pert_top`, `make_pt_top`, `top2build` |
| Box/coord preparation | `build_box`, `sim_box`, `ran_box`, `ran_solvation`, `ion`, `copy_box`, `bin_box`, `inbox`, `unify_box`, `check_box` |
| Coordinate conversion | `pdb2g96`, `traj2pdb`, `amber2gromos` |
| Script generation | `mk_script`, `mpi_scaling` |
| Simulation prep | `prep_eds`, `prep_noe`, `prep_xray` |
| Utilities | `atominfo`, `gch` |

**gromos-analysis** (72 bins) — trajectory analysis & post-processing:

| Category | Programs |
|----------|----------|
| Structural | `rmsd`, `rmsf`, `rmsdmat`, `rmsdhel`, `rgyr`, `dssp`, `fitmol`, `cog`, `distmat` |
| Energy | `ene_ana`, `ener`, `ener_box`, `int_ener`, `temperature`, `pot_aver`, `propertyaver` |
| Pair/distribution | `rdf`, `rdf_matrix`, `hbond`, `contactnum`, `close_pair`, `iondens`, `sasa`, `sasa_hasel` |
| Dynamics/transport | `diffus`, `visco`, `tcf`, `angaver`, `rot_rel`, `dipole`, `epsilon`, `eps_field` |
| Free energy | `bar`, `ext_ti_ana`, `trs_ana`, `reweight`, `rep_ana`, `eds_ana`, `gamd_ana` |
| Trajectory manipulation | `frameout`, `filter`, `tstrip`, `gathtraj`, `follow`, `tser` |
| NOE/NMR | `noe`, `post_noe`, `predict_noe`, `nhoparam` |
| Clustering | `cluster`, `postcluster` |
| X-ray/crystallography | `r_factor`, `r_real_factor`, `structure_factor`, `xray_map`, `xrayts`, `cry`, `cry_rms` |
| Special | `disicl`, `ditrans`, `solute_entropy`, `m_widom`, `shake_analysis`, `swd`, `epath`, `pairlist` |

**Migration plan:** ✅ DONE
- [x] Created `crates/gromos-md/` with simulation engines (8 bins) + integration tests
- [x] Created `crates/gromos-tools/` with system construction bins (30 bins, subdirs: topology/box/conversion/utilities)
- [x] Moved analysis bins into `crates/gromos-analysis/` as `[[bin]]` targets (66 bins, subdirs: structural/energy/distribution/dynamics/free_energy/trajectory/noe/clustering/xray/special)
- [x] Removed `gromos-analysis` from `gromos` facade to avoid circular dependency; analysis bins use `gromos_core`/`gromos_io` directly
- [x] `crates/gromos-cli/` kept as thin unified `gromos` multicall binary (clap only)
- [x] Shared library code stays in gromos-io / gromos-core
- [x] gromos-analysis library implementations (`rdf`, `rmsd`, `hbond`, `gyration`, `diffusion`) available for bins to call into
- [x] All 21 reference tests pass, workspace compiles clean

**Target directory structure:**

```
crates/gromos-md/
├── Cargo.toml
├── src/
│   ├── bin/
│   │   ├── md.rs              (main MD simulation driver)
│   │   ├── md_mpi.rs
│   │   ├── md_mpi_cuda.rs
│   │   ├── mdf.rs
│   │   ├── remd.rs
│   │   ├── repex_mpi.rs
│   │   ├── eds.rs
│   │   └── gamd.rs
│   └── lib.rs                 (shared: algorithm sequence builder, CLI @-convention)
└── tests/
    ├── test_gromosXX_references.rs
    └── gromosXX_references/

crates/gromos-tools/
├── Cargo.toml
├── src/
│   ├── bin/
│   │   ├── topology/
│   │   │   ├── make_top.rs
│   │   │   ├── com_top.rs
│   │   │   ├── con_top.rs
│   │   │   ├── link_top.rs
│   │   │   ├── red_top.rs
│   │   │   ├── check_top.rs
│   │   │   ├── addvirt_top.rs
│   │   │   ├── pert_top.rs
│   │   │   ├── make_pt_top.rs
│   │   │   └── top2build.rs
│   │   ├── box/
│   │   │   ├── build_box.rs
│   │   │   ├── sim_box.rs
│   │   │   ├── ran_box.rs
│   │   │   ├── ran_solvation.rs
│   │   │   ├── ion.rs
│   │   │   ├── copy_box.rs
│   │   │   ├── bin_box.rs
│   │   │   ├── inbox.rs
│   │   │   ├── unify_box.rs
│   │   │   └── check_box.rs
│   │   ├── conversion/
│   │   │   ├── pdb2g96.rs
│   │   │   ├── traj2pdb.rs
│   │   │   └── amber2gromos.rs
│   │   └── utilities/
│   │       ├── mk_script.rs
│   │       ├── mpi_scaling.rs
│   │       ├── atominfo.rs
│   │       ├── gch.rs
│   │       ├── prep_eds.rs
│   │       ├── prep_noe.rs
│   │       └── prep_xray.rs
│   └── lib.rs                 (shared: building block parser, MTB/IFP readers, box algorithms)
└── tests/

crates/gromos-analysis/
├── Cargo.toml
├── src/
│   ├── lib.rs                 (library: reusable analysis algorithms)
│   ├── rdf.rs
│   ├── rmsd.rs
│   ├── hbond.rs
│   ├── gyration.rs
│   ├── diffusion.rs
│   ├── clustering.rs
│   ├── fitting.rs             (rotational fit, Kabsch)
│   ├── selection.rs           (atom selection grammar)
│   ├── timeseries.rs          (statistics, autocorrelation, block averaging)
│   └── bin/
│       ├── structural/
│       │   ├── rmsd.rs
│       │   ├── rmsf.rs
│       │   ├── rmsdmat.rs
│       │   ├── rmsdhel.rs
│       │   ├── rgyr.rs
│       │   ├── dssp.rs
│       │   ├── fitmol.rs
│       │   ├── cog.rs
│       │   └── distmat.rs
│       ├── energy/
│       │   ├── ene_ana.rs
│       │   ├── ener.rs
│       │   ├── ener_box.rs
│       │   ├── int_ener.rs
│       │   ├── temperature.rs
│       │   ├── pot_aver.rs
│       │   └── propertyaver.rs
│       ├── distribution/
│       │   ├── rdf.rs
│       │   ├── rdf_matrix.rs
│       │   ├── hbond.rs
│       │   ├── contactnum.rs
│       │   ├── close_pair.rs
│       │   ├── iondens.rs
│       │   ├── sasa.rs
│       │   └── sasa_hasel.rs
│       ├── dynamics/
│       │   ├── diffus.rs
│       │   ├── visco.rs
│       │   ├── tcf.rs
│       │   ├── angaver.rs
│       │   ├── rot_rel.rs
│       │   ├── dipole.rs
│       │   ├── epsilon.rs
│       │   └── eps_field.rs
│       ├── free_energy/
│       │   ├── bar.rs
│       │   ├── ext_ti_ana.rs
│       │   ├── trs_ana.rs
│       │   ├── reweight.rs
│       │   ├── rep_ana.rs
│       │   ├── eds_ana.rs
│       │   └── gamd_ana.rs
│       ├── trajectory/
│       │   ├── frameout.rs
│       │   ├── filter.rs
│       │   ├── tstrip.rs
│       │   ├── gathtraj.rs
│       │   ├── follow.rs
│       │   └── tser.rs
│       ├── noe/
│       │   ├── noe.rs
│       │   ├── post_noe.rs
│       │   ├── predict_noe.rs
│       │   └── nhoparam.rs
│       ├── clustering/
│       │   ├── cluster.rs
│       │   └── postcluster.rs
│       ├── xray/
│       │   ├── r_factor.rs
│       │   ├── r_real_factor.rs
│       │   ├── structure_factor.rs
│       │   ├── xray_map.rs
│       │   ├── xrayts.rs
│       │   ├── cry.rs
│       │   └── cry_rms.rs
│       └── special/
│           ├── disicl.rs
│           ├── ditrans.rs
│           ├── solute_entropy.rs
│           ├── m_widom.rs
│           ├── shake_analysis.rs
│           ├── swd.rs
│           ├── epath.rs
│           └── pairlist.rs
└── tests/
```

### DONE — Fix nacl_3water_box (SHAKE multi-solvent) ✓
- [x] Added solvent constraint loop using `solvent_constraint_template` over all molecules

### DONE — Twin-range pairlist (NSNB>1) ✓
- [x] Wire RCUTP as short-range cutoff, long-range forces cached on NSNB update steps

### DONE — Fix nacl_water_box CRF mismatch ✓
- [x] Root cause: two issues in nonbonded innerloop
  1. Solvent long-range pairlist stored only first-atom pairs and used `solvent_innerloop` with shared
     PBC shift. At CG distances 0.8–0.9 nm, H atoms can land in different periodic images than the
     O-O shift gives. Fix: expand long-range solvent CG pairs to all 9 atom pairs, use `lj_crf_innerloop`
     with per-atom `nearest_image`.
  2. `lj_crf_innerloop` had HEAVISIDE truncation (`if r2 > cutoff2 { continue; }`) but gromosXX
     reference build uses `#undef XXHEAVISIDE` (default). Fix: removed HEAVISIDE check.
- [x] nacl_water_box CRF: -103.572 → -107.392 (ref=-107.394, delta=0.002 kJ/mol)
- [x] nacl_water_box_shifted and nacl_3water_cutoff also now pass
- [x] All 11 previously passing tests still pass (no regressions)

### DONE — Fix water_216_box kinetic energy mismatch ✓
- [x] Root cause: NTIVEL=1 + COM removal (see Resolved Investigations below)
- [x] water_216_box now passes all 10 frames within 5e-7 kJ/mol

### DONE — Fix aladip_vacuum / aladip_solvated (improper dihedral energy) ✓
- [x] Root cause: improper dihedral force constant K not converted from kJ/(mol·deg²) to kJ/(mol·rad²)
  - gromosXX multiplies K by (180/π)² ≈ 3283 during topology parsing (in_topology.cc:1055)
  - gromos-rs was storing K raw, giving energies ~3283× too small
- [x] Fix: multiply CQ by 180²/π² in parse_improper_dihedral_types()
- [x] Both aladip_vacuum and aladip_solvated now pass

### DONE — COM motion removal → unblocks `water_216_box_com` ✓
- [x] Implement `RemoveCOMMotion` algorithm
  - [x] Parse NTICOM from INITIALISE block (line 1, field 3)
  - [x] Parse NSCM from COMTRANSROT block (already done)
  - [x] `remove_com_translation()`: COM_v = Σ(m_i·v_i)/Σ(m_i), then v_i -= COM_v
  - [x] Wire into AlgorithmSequence as FIRST algorithm (before Forcefield, gromosXX convention)
  - [x] At step 0: apply if NTICOM>=1 (initial removal)
  - [x] At step > 0: apply every NSCM steps when NSCM>0
  - [x] Only modifies current().vel (gromosXX convention); placement before exchange_state ensures both states are COM-free
  - [x] gromosXX ref: `algorithm/constraints/remove_com_motion.cc`, `algorithm/create_md_sequence.cc`
- [x] Reference test: `water_216_box_com` (NTICOM=1, NSCM=10) passes

### DONE — SHAKE constraint force, virial, NTC, skip optimization ✓
- [x] Added `NtcMode` enum (SolventOnly, HydrogenBonds, AllBonds) for NTC 1–3 control
- [x] Added `ntc` field to `ShakeParameters`
- [x] Constraint force accumulation: f(c) = (r_new - r_uc) · m / dt², stored in `conf.old().constraint_force`
- [x] Virial tensor contribution: ref_r ⊗ ref_r · λ/dt² added to virial
- [x] skip_now/skip_next optimization: converged constraints skipped on next iteration
- [x] `shake_positions()` and `shake_velocities()` for initial configuration shaking
- [x] `shake_algorithm.rs`: `shake_initial_positions`/`shake_initial_velocities` flags, `init()` method
- [x] `md.rs`: NTC mode from `imd.ntc`, NTISHK controls initial shaking (0=none, 1=pos, 2=vel, 3=both)
- [x] Constraint forces passed to force trajectory writer (`@trf`)
- [x] All 16 tests pass, zero regressions

### DONE — Berendsen thermostat → unblocks `water_216_nvt` ✓
- [x] Implemented `BerendsenThermostat` as Algorithm trait in `berendsen_thermostat.rs`
  - [x] Single-bath support with configurable T0, τ, DOF
  - [x] Scale velocities: λ = sqrt(1 + dt/τ · (T₀/T_free - 1)), v *= λ
  - [x] T_free = 2·E_kin_new / (DOF·k_B) using "new" kinetic energy (gromosXX: multibath.bath.ekin)
  - [x] τ < 0 → no coupling (NVE), τ = 0 → instantaneous scaling
  - [x] Placed between LeapFrogVelocity and LeapFrogPosition (gromosXX convention)
- [x] Added `kinetic_energy_new` field to Energy struct for thermostat scaling
- [x] `TemperatureCalculation::init()` calls `apply()` to pre-compute initial E_kin_new
  - This matches gromosXX `Temperature_Calculation::init()` behavior — NOT a bug:
    the thermostat needs real kinetic energy data at step 0 to compute the first
    scaling factor. Without it, scale=1.0 at step 0 and the error cascades.
- [x] DOF = 3·N_atoms - N_solvent_constraints - NDFMIN
- [x] Reference test: `water_216_nvt` (TAU=0.1, weak-coupling) passes all 4 frames

### DONE — Virial / pressure architecture (matching gromosXX) ✓
- [x] Constraint virial contribution now computed in SHAKE (ref_r ⊗ ref_r · λ/dt²)
- [x] Refactored to match gromosXX architecture:
  - [x] `prepare_virial` in Forcefield: computes KE tensor on `conf.current()` at start of force calc
    - Molecular: COM velocities per pressure group (KE_ij = 0.5 * v_COM_i * v_COM_j / M_total)
    - Atomic: per-atom KE_ij = 0.5 * Σ m_k * v_ki * v_kj
  - [x] `atomic_to_molecular_virial` in Forcefield: corrects virial at end of force calc
    - For each pressure group: COM position via chain gather, corrP(b,a) += f(a) * r(b), virial -= corrP
  - [x] `PressureCalculation` simplified: only computes P_tensor = (KE + 0.5*virial) * 2/V from conf.old()
- [x] `nearest_image` already exists in math.rs (Periodicity enum with Rectangular/Triclinic/Vacuum)
- [x] Pressure groups infrastructure:
  - [x] Changed `topo.pressure_groups` from `Vec<Vec<usize>>` to `Vec<Range<usize>>` (like molecules)
  - [x] Added PRESSUREGROUPS parser in topology reader (cumulative boundary vector format)
  - [x] Solvent molecules added as individual pressure groups during `solvate()`
  - [x] Fallback: all solute atoms as one pressure group if PRESSUREGROUPS block absent
- [x] `Forcefield.virial_type` field wired from IMD PRESSURESCALE.VIRIAL setting

### DONE — Berendsen barostat → unblocks `water_216_npt` ✓
- [x] `BerendsenBarostat` algorithm: isotropic scaling μ = (1 - κ·dt/τ·(P₀-P))^(1/3)
- [x] Wired in md.rs: PressureCalculation → BerendsenBarostat in algorithm sequence
- [x] PRESSURESCALE parser: keyword format (COUPLE/SCALE/COMP/TAUP/VIRIAL/SEMI/pres0)
- [x] **Fixed: missing long-range virial in twin-range.**  When twin-range is active
      (RCUTP < RCUTL), the long-range ForceStorage.virial was computed but never added
      to the main nonbonded_storage. This caused the virial trace to be ~1000 kJ/mol
      too small, producing wildly wrong pressures and barostat box-collapse.
      Fix: cache `longrange_virial` alongside forces/energies, add it each step.
- [x] **Fixed: stale periodicity under NPT.**  `Forcefield.periodicity` (Rectangular struct
      caching box_size, half_box, inv_box) was set once at construction and never updated.
      After the barostat scales the box, all subsequent PBC nearest_image calls used the
      old box dimensions, causing subtly wrong forces and accumulating energy drift.
      Fix: refresh `self.periodicity` from `conf.current().box_config` at the start of
      every `Forcefield.apply()` call (guarded for vacuum/zero-dim boxes).
- [x] Reference test: `water_216_npt` passes all 4 frames

### TODO — Triclinic box support (lower priority)
- [ ] Code exists in `math.rs` but md.rs never creates Triclinic periodicity
  - [ ] Wire GENBOX box_type into periodicity selection
  - [ ] Test with truncated octahedron or other non-rectangular boxes

### DONE — 1-4 interactions (cs6/cs12 + INE14 parsing) ✓
- [x] Parse CS12/CS6 (columns 5-6) from LJPARAMETERS block in topology
- [x] Parse INE14 line from SOLUTEATOM block → stored in `Topology.one_four_pairs`
- [x] Added `LJParameters::new_with_14(c6, c12, cs6, cs12)` constructor
- [x] Separated exclusion lists (gromosXX convention):
  - `exclusions` = 1-2 + 1-3 (used for RF excluded interactions)
  - `is_excluded_or_14()` = 1-2 + 1-3 + 1-4 (used for pairlist exclusion)
- [x] Implemented `one_four_interaction_loop` in nonbonded.rs:
  - Uses cs6/cs12 for LJ, coulomb_scaling for CRF, no cutoff, full virial
- [x] Wired into forcefield.rs after RF excluded interactions
- [x] Added cs6/cs12 fields to nonbonded `LJParameters` struct
- [x] Created `butane_vacuum` test system (4 UA atoms, gauche conformation)
- [x] Reference tests: butane_vacuum and benzene_vacuum pass

### TODO — Unit conversion audit (topology parsing)
- [x] IMPDIHEDRALTYPE: CQ converted from kJ/(mol·deg²) to kJ/(mol·rad²) — ×(180/π)²
- [x] BONDANGLEBENDTYPE: CHT (k_harmonic) converted — ×(180/π)²
- [ ] HARMBONDANGLETYPE: if/when parsed, must also convert CHT — ×(180/π)²
- gromosXX converts ALL angle force constants at parse time (in_topology.cc)
- No conversion needed for: CT (cosine angle K), CP (dihedral K), bond K, LJ, charges
- Reference: gromosXX in_topology.cc lines 854, 928, 1055

### TODO — Code quality & consistency
- [ ] Fix clippy warnings (~390 total: gromos-forces 89, gromos-integrators 77, gromos-io 31, gromos-core 15)
  - Mechanical fixes (auto-fixable): borrowed refs, double parens, legacy constants, assign ops
  - Manual: unused variables, dead imports, loop indexing → iterators
- [ ] Run `cargo clippy --fix --workspace` for auto-fixable warnings, review manually after
- [x] Fix test compilation: integer-to-float casts (`i as f64`), f32→f64 literal suffixes,
      cs6/cs12 fields in LJParameters test constructors, IMD block format (keyword→positional),
      nearest_image test assertions (ri-rj direction), bspline loop bound cast
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

##### DONE — Compositional Simulation API (OpenMM-inspired) ✓
- [x] `Topology("system.topo")` — wraps Rust Topology, exposes n_atoms, masses, charges, solvate()
- [x] `Configuration("initial.cnf")` — wraps coordinate data, exposes positions/velocities/box as numpy
- [x] `InputParameters("run.imd")` — wraps ImdParameters, exposes dt, nstlim, nsm, cutoff, temperature, ntc, ntb, ntwx, ntwe
- [x] `Simulation(topo, conf, params)` — compositional constructor from Python objects
- [x] `Simulation("file.topo", "file.cnf", "file.imd")` — string-path constructor (backward-compatible)
- [x] `Simulation.from_files(...)` — explicit static method alternative
- [x] `sim.algorithm_names` — read-only list of algorithms in the MD sequence
- [x] Shared `build_simulation()` function — single source of truth for both constructor paths
- [x] Solvation handled automatically: topology solvated only if not already solvated
- [x] All 62 Python reference tests pass (21 systems × 3 checks, 1 skip)

##### DONE — Python reference tests via Simulation API ✓
- [x] `test_gromosXX_references.py` — uses compositional API (Topology, Configuration, InputParameters)
- [x] Energy, position, and force comparison against gromosXX reference data
- [x] Write frequency parsing via `InputParameters.nstlim`/`.ntwe`/`.ntwx` (replaces regex)
- [x] Minimum image convention for position comparison (handles periodic wrapping)
- [x] 62 passed, 1 skipped (water_216_box_com positions — no reference trajectory)

##### DONE — Compositional AlgorithmSequence API ✓
- [x] `AlgorithmSequence` — Python-constructable, ordered list of algorithm descriptors
- [x] Individual algorithm descriptor classes: `Forcefield`, `LeapFrogIntegrator`,
      `LeapFrogVelocity`, `LeapFrogPosition`, `BerendsenThermostat`, `BerendsenBarostat`,
      `ShakeConstraints`, `TemperatureCalculation`, `PressureCalculation`, `EnergyCalculation`,
      `RemoveCOMMotion`
- [x] Each algorithm has typed keyword parameters with defaults from `InputParameters`
- [x] Sequence manipulation: `add()`, `insert_after()`, `insert_before()`, `remove()`, `replace()`
- [x] Preset constructors: `AlgorithmSequence.nve(topo, params)`, `.nvt()`, `.npt()`, `.from_parameters()`
- [x] `Simulation.from_sequence(topo, conf, params, sequence)` — creates simulation with custom sequence
- [x] Advanced Leap-Frog split: `LeapFrogVelocity` + `LeapFrogPosition` for thermostat placement
- [x] Descriptor resolution: lightweight Python objects → real Rust algorithms at Simulation creation
- [x] `__repr__`, `__len__`, `__contains__`, `.names` for inspection
- [x] Validated: identical results between standard `Simulation()` and `Simulation.from_sequence()`

##### DONE — Type stubs for static analysis ✓
- [x] `gromos.pyi` stub file: type signatures for all PyO3 bindings (Vec3, Energy, Frame,
      Topology, Configuration, InputParameters, Simulation, AlgorithmSequence, all algorithm
      descriptor classes, rmsd, rdf)
- [x] Resolves `ty`/`mypy`/`pyright` "Cannot resolve imported module `.gromos`" error
- [x] Pattern follows Polars: `.pyi` stub alongside native `.abi3.so` extension module

##### TODO — Additional binding types
- [ ] Expose ForceField evaluation (single-point energy/force calculation)
- [ ] Expose SHAKE / constraint info
- [ ] Expose energy decomposition (bonded, LJ, CRF, kinetic, pressure)
- [ ] Study Polars' pyo3 patterns: `PyDataFrame`, `PyExpr`, `PyLazyFrame` wrappers
  - Path: `.local/polars/py-polars/src/` — how they wrap Rust types into Python classes
  - Learn from: `__repr__`, `_repr_html_`, method chaining, `@staticmethod` constructors

#### Phase 2 — Python API (py-gromos)
- [x] `gromos.Topology("file.topo")` → Topology object with `__repr__`, masses, charges
- [x] `gromos.Configuration("file.cnf")` → Configuration with positions as numpy arrays
- [x] `gromos.Simulation(topo, conf, params)` — high-level compositional runner
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

## Resolved Investigations

### water_216_box E_kin mismatch — FIXED ✓

**Symptom:** E_kin=139.3 vs reference E_kin=2605.1 at step 0. E_pot matched exactly.

**Root cause:** Two interacting issues:
1. **NTIVEL=1 (runtime velocity generation).** The original .in files had NTIVEL=1, telling gromosXX
   to generate Maxwell-Boltzmann velocities at runtime using GSL MT19937 RNG (seed=210185, T=300K).
   gromos-rs doesn't implement NTIVEL=1, so it read the all-zero velocities from the .conf file.
2. **COM motion removal.** gromosXX removes initial COM translation (NTICOM=1) and periodic COM
   motion (NSCM=10) during init and the MD loop. gromos-rs doesn't implement COM removal yet.
   Even with correct velocities, this caused a ~1.2 kJ/mol E_kin discrepancy.

**Fix:**
- Pre-generated exact Maxwell-Boltzmann velocities using a C program linked against GSL,
  matching gromosXX's RNG output bit-for-bit (MT19937, `gsl_ran_gaussian` polar Box-Muller).
- Embedded velocities in `water_216_box.conf` using scientific notation (`%.15e`) with proper
  24-character label alignment (gromosXX skips first 24 chars of VELOCITY lines).
- Set NTIVEL=0 so both codes read velocities from file.
- Set NTICOM=0 and NSCM=0 to disable COM removal in both codes.
- Regenerated gromosXX references with these settings.

**Files changed:**
- `gromosXX_references/water_216_box/water_216_box.conf` — full-precision velocities embedded
- `gromosXX_references/water_216_{box,nvt,npt}/*.in` — NTIVEL=0, NTICOM=0, NSCM=0
- `gromosXX_references/water_216_{box,nvt,npt}/expected/` — regenerated from gromosXX

**Result:** water_216_box passes all 10 frames, max abs diff = 5e-7 kJ/mol.
water_216_nvt/npt still fail because they need thermostat/barostat (separate TODO).

### nacl_water_box CRF mismatch — FIXED ✓

**Symptom:** E_pot off by ~3.82 kJ/mol at step 0. LJ matched exactly. CRF was the sole discrepancy.

**Root cause:** Two issues in the nonbonded innerloop:
1. **Solvent long-range pairlist used shared PBC shift.** The pairlist stored only first-atom (O-O) pairs
   for solvent CGs and `solvent_innerloop` expanded them using the O-O periodic shift for all 9 atom pairs.
   At CG distances 0.8–0.9 nm, H atoms can be in different periodic images than what the O-O shift gives.
   gromosXX standard mode (no fast SPC loops) uses per-atom `nearest_image` for all pairs.
   **Fix:** Expand long-range solvent CG pairs to all 9 atom pairs in the pairlist; process with
   `lj_crf_innerloop` (per-atom nearest_image) instead of `solvent_innerloop`.
2. **HEAVISIDE truncation was enabled but gromosXX reference has it disabled.** `lj_crf_innerloop` had
   `if r2 > cutoff2 { continue; }` (the `XXHEAVISIDE` check), but gromosXX was built with
   `#undef XXHEAVISIDE` (default in config.h). This caused expanded long-range atom pairs beyond
   the cutoff to be skipped when they shouldn't be.
   **Fix:** Removed HEAVISIDE check from `lj_crf_innerloop`.

**Files changed:**
- `crates/gromos-core/src/pairlist.rs` — solvent long-range now expands CG pairs to all atom pairs
- `crates/gromos-forces/src/nonbonded.rs` — removed HEAVISIDE from `lj_crf_innerloop`
- `crates/gromos-integrators/src/algorithms/forcefield.rs` — long-range solvent uses `lj_crf_innerloop`

**Result:** nacl_water_box CRF delta dropped from 3.82 kJ/mol to 0.002 kJ/mol.
Also fixed nacl_water_box_shifted and nacl_3water_cutoff. No regressions on other tests.

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