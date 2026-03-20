# gromos-rs — Status & Plan

we will focus on `md`, so therefore focus on `cargo build --release --bin md`

letsgo reviewing things with compare references
modify `PLAN.md` in case you advance
in case we commit update CHANGELOG.md and update the Cargo.toml version

DONT MODIFY THE REFERENCES in gromosXX_references/*/expected/ — those are gromosXX ground truth.
To regenerate references: `python3 gromosXX_references/run_references.py --md-binary /path/to/gromosXX/md`
New reference systems should be added via run_references.py (add to SYSTEMS list).

Validation is done via Rust integration tests (replaces the old compare_reference.py):
  `cargo test -p gromos-cli --test test_gromosXX_references`
New test systems: add a `ref_test!(name, "system_dir")` line in:
  `crates/gromos-cli/tests/test_gromosXX_references.rs`

Here you can find the gromosXX code: .local/gromosXX/md++/src
Here you can find realistic tutorials: .local/gromos_tutorial_livecoms/tutorial_files
Here you can find theory behind everything: .local/doc/gromos_book


## Architecture

```
gromos-core      → types, math, Vec3, BoundaryCondition trait, Topology, Configuration, Algorithm trait
gromos-forces    → bonded (quartic/harmonic bonds, angles, dihedrals, impropers, cross-dihedrals)
                   nonbonded (LJ + CRF, rf_excluded_interactions, pairlist)
                   PME, QM/MM, restraints, polarization (stubs/partial)
gromos-integrators → LeapFrog, VelocityVerlet, StochasticDynamics, SteepestDescent
                     SHAKE, SETTLE, LINCS constraints
                     Berendsen thermostat/barostat, Nosé-Hoover, Andersen, Parrinello-Rahman
                     EDS, GaMD, REMD, FEP
gromos-io        → topology parser, coordinate reader (GENBOX, POSITION, VELOCITY),
                   IMD parser, trajectory/energy/force writers
gromos-cli       → md binary (main simulation driver)
gromos-analysis  → RDF, RMSD, gyration, MSD, H-bonds
pyo3-gromos      → Python bindings
py-gromos        → Python package (separate build)
```

## MD Loop (AlgorithmSequence)

```
Forcefield → LeapFrogVelocity → LeapFrogPosition → SHAKE (if ntc>1) → TemperatureCalculation → EnergyCalculation
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

Validation: `cargo test -p gromos-cli --test test_gromosXX_references`
Generate refs: `python3 crates/gromos-cli/tests/run_references.py --md-binary .local/gromosXX/md++/build/program/md`
Test file: `crates/gromos-cli/tests/test_gromosXX_references.rs`
Ref data: `crates/gromos-cli/tests/gromosXX_references/`

| Lvl | System           | Atoms | Isolates                              | Status   |
|-----|------------------|-------|---------------------------------------|----------|
| 0   | pair_lj          | 2     | Pure LJ, no PBC                      | **PASS** |
| 0   | pair_lj_mixed    | 2     | LJ combination rules                 | **PASS** |
| 0   | nacl_pair        | 2     | Coulomb + LJ ions                    | **PASS** |
| 1   | water_single     | 3     | bond + angle + intramolecular CRF    | **PASS** |
| 1   | benzene_vacuum   | 12    | aromatic ring + improper + torsion   | **PASS** |
| 1   | nacl_pair_box    | 2     | Coulomb + LJ in PBC with RF (no solvent) | **PASS** |
| 1   | aladip_vacuum    | 12    | all bonded + exclusions + 1-4        | FAIL — missing topo file |
| 2   | water_3_box      | 9     | PBC + min image + pairlist + CRF     | **PASS** |
| 2   | nacl_1water_box  | 5     | minimal solute-solvent + SHAKE       | **PASS** |
| 2   | nacl_3water_box  | 11    | multiple solvent + solute-solvent pairlist | **PASS** |
| 2   | water_3_box_twinrange | 9 | twin-range pairlist (RCUTP<RCUTL, NSNB=5) | **PASS** |
| 2   | water_10_box     | 32    | 2 ions + 10 SPC, positions away from cutoff | **PASS** |
| 2   | nacl_3water_cutoff | 11  | nacl_3water near cutoff boundary     | **PASS** |
| 2   | nacl_water_box   | 62    | ion-water RF in PBC                  | **PASS** |
| 2   | nacl_water_box_shifted | 62 | nacl_water_box with perturbed positions | **PASS** |
| 3   | water_216_box    | 648   | bulk NVE, pairlist, virial           | FAIL — E_kin wrong |
| 3   | water_216_nvt    | 648   | Berendsen thermostat                 | FAIL — E_kin wrong |
| 3   | water_216_npt    | 648   | Berendsen barostat                   | FAIL — E_kin wrong |
| 4   | aladip_solvated  | 72    | SHAKE + solute-solvent               | FAIL — missing topo file |

**14 of 19 tests pass.** Levels 0-2 fully passing.

## What Works

- LJ + CRF nonbonded forces (vacuum and PBC code paths)
- All bonded force types: quartic/harmonic bonds, cos-harmonic/harmonic angles, dihedrals, impropers, cross-dihedrals
- NTF flag control for selective bonded force evaluation
- RF excluded interactions (forces + energy + self-terms)
- Pairlist: StandardPairlistAlgorithm with configurable update frequency
  - Chargegroup-based mode: CG center-of-geometry distance check, atom pairs between CGs
  - Atom-based mode: direct atom-atom distance check
  - Twin-range: short (RCUTP) and long (RCUTL) pairlists with force caching
- Boundary conditions: Vacuum, Rectangular (minimum image), Triclinic (defined but not wired)
- SHAKE constraints: solute bonds (NTC>1) + solvent constraints (NTCS>0)
  - Iterates SOLVENTCONSTR template over all solvent molecules
  - Velocity correction applied to all constrained atoms
- SETTLE, LINCS constraints (implemented but not wired)
- Berendsen thermostat & barostat (parsed from MULTIBATH/PRESSURESCALE IMD blocks)
- Topology parser: all block types including SOLUTEATOM exclusions, SOLVENTCONSTR, LJPARAMETERS
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

### Immediate — Fix aladip_vacuum / aladip_solvated (missing topology)
- [ ] `aladip.topo` referenced as `../../aladip.topo` but doesn't exist
  - Either the topo file was never committed, or the relative path is wrong
  - Need to locate or regenerate the aladip topology file
  - Both aladip_vacuum and aladip_solvated depend on it

### Immediate — Fix water_216 kinetic energy mismatch
- [ ] water_216_box, water_216_nvt, water_216_npt all fail with E_kin wrong
  - E_pot matches perfectly at step 0: ours=-6876.5, ref=-6876.5 ✓
  - E_kin way off: ours=139.3, ref=2605.1 (should be ~300K worth)
  - Our TemperatureCalculation uses `0.5 * m * (|v_new|² + |v_old|²) / 2`
  - At step 0, `v_old` = initial velocities (from conf at t=0), `v_new` = v after first leap-frog step
  - gromosXX reports E_kin=2605 at step 0, which corresponds to the input velocities at full T=300K
  - Need to check how gromosXX handles step 0: likely uses input velocities for both old/new,
    or computes E_kin before the first velocity update
  - Check `molecular_translational_ekin()` in `.local/gromosXX/md++/src/configuration/state_properties.cc`
  - water_216 has NSM=0 (all 648 atoms are SOLUTEATOM), no solvent expansion
- [ ] water_216_nvt additionally needs Berendsen thermostat wired and working
- [ ] water_216_npt additionally needs Berendsen barostat + pressure coupling

### Known Gaps
- Triclinic box: code exists in math.rs but md.rs never creates Triclinic periodicity
- 1-4 interactions: parsed (INE14) but need verification in force loop
- Virial / pressure calculation: needs verification for NPT
- COM motion removal (COMTRANSROT): parsed but not wired

## Resolved Investigations

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
crates/gromos-cli/src/bin/md.rs              — main MD driver, CLI, simulation setup
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
crates/gromos-cli/tests/test_gromosXX_references.rs — integration tests vs gromosXX
crates/gromos-cli/tests/run_references.py           — generate gromosXX reference data
crates/gromos-cli/tests/gromosXX_references/        — reference input + expected output
```

## How to Test

```sh
# Build
cargo build --release --bin md

# Run integration tests against gromosXX references (14 active, 5 ignored)
cargo test -p gromos-cli --test test_gromosXX_references

# Run including known-failing systems
cargo test -p gromos-cli --test test_gromosXX_references -- --include-ignored

# Run a specific system
cargo test -p gromos-cli --test test_gromosXX_references -- pair_lj --exact

# Regenerate gromosXX reference data (requires gromosXX md++ binary)
python3 crates/gromos-cli/tests/run_references.py --md-binary .local/gromosXX/md++/build/program/md

# Add a new reference system:
# 1. Add to SYSTEMS list in crates/gromos-cli/tests/run_references.py
# 2. Run run_references.py to generate expected/ data
# 3. Add ref_test!(name, "dir") in crates/gromos-cli/tests/test_gromosXX_references.rs
```