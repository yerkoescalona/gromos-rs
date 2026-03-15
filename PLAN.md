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
| 1   | aladip_vacuum    | 12    | all bonded + exclusions + 1-4        | TODO     |
| 2   | water_3_box      | 9     | PBC + min image + pairlist + CRF     | **PASS** |
| 2   | nacl_1water_box  | 5     | minimal solute-solvent + SHAKE       | **PASS** |
| 2   | nacl_3water_box  | 11    | multiple solvent + solute-solvent pairlist | **PASS** |
| 2   | water_3_box_twinrange | 9 | twin-range pairlist (RCUTP<RCUTL, NSNB=5) | **PASS** |
| 2   | nacl_water_box   | 62    | ion-water RF in PBC                  | FAIL     |
| 3   | water_216_box    | 648   | bulk NVE, pairlist, virial           | FAIL     |
| 3   | water_216_nvt    | 648   | Berendsen thermostat                 | TODO     |
| 3   | water_216_npt    | 648   | Berendsen barostat                   | TODO     |
| 4   | aladip_solvated  | 72    | SHAKE + solute-solvent               | TODO     |

Levels 0-1 pass — pair forces, bonded, and PBC match gromosXX to ~1e-9.
Level 2 water_3_box passes. nacl_1water_box passes (1 solvent molecule + SHAKE).
nacl_3water_box: **FIXED** — SHAKE was only iterating solute bonds, not solvent constraints.
  Now applies SOLVENTCONSTR to all solvent molecules. All 10 steps match to ~1e-10.
water_3_box_twinrange: **FIXED** — twin-range force caching implemented.
  RCUTP wired as short-range cutoff, long-range forces cached on NSNB update steps.
  All 10 steps match gromosXX to ~5e-11.
nacl_water_box runs but E_pot off by ~3.82 kJ/mol (ours=-112.57, ref=-116.39).
  - LJ energy matches: our=-8.998, ref=-8.998 ✓
  - CRF mismatch: ~3.82 kJ/mol — likely `put_chargegroups_into_box` issue
  - With RCUTP=0.8 now wired: solute_short=12, solute_long=24, solvent_short=70, solvent_long=19
  - RF excluded already split: solute gets self+full, solvent gets only DD term ✓
  - See **Investigation: nacl_water_box CRF mismatch** below
Level 3 water_216_box fails (large energy discrepancy — all atoms are SOLUTEATOM, NSM=0).

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
- Energy/trajectory/force writers (ENERTRJ, POSITIONRED blocks)

## TODO

### DONE — Fix nacl_3water_box (SHAKE multi-solvent) ✓
- [x] **Root cause:** SHAKE only iterated `topo.solute.bonds`, never applied solvent constraints
- [x] **Fix:** Added solvent constraint loop using `solvent_constraint_template` over all molecules
- [x] **Fix:** `shake_enabled` now checks `imd.ntcs > 0 && imd.nsm > 0` (not just `ntc > 1`)
- [x] Refactored into `shake_one_constraint()` helper for code reuse

### DONE — Twin-range pairlist (NSNB>1) ✓
- [x] Wire RCUTP as short-range cutoff (`PairlistContainer::new(imd.rcutp, imd.rcutl, 0.0)`)
- [x] Forcefield detects twin-range: `twin_range_active = rcutp < rcutl - 1e-10`
- [x] On pairlist update: evaluate long-range pairs, cache forces + energies
- [x] On non-update steps: add cached long-range forces/energies to short-range
- [x] Non-twin-range systems (RCUTP==RCUTL): evaluate long+short together as before
- [x] All 10 steps of water_3_box_twinrange match gromosXX to ~5e-11

### Immediate — Fix nacl_water_box CRF mismatch
- [ ] Fix nacl_water_box E_pot mismatch (~3.82 kJ/mol off)
  - [x] Implement gromosXX pairlist architecture: solute vs solvent separation
  - [x] Add intra-CG non-excluded pairs to pairlist (solute CGs only)
  - [x] Solvent CG distance: use first-atom position, not center-of-geometry
  - [x] Solvent innerloop: compute nearest_image once per O-O pair, reuse shift for all 9 pairs
  - [x] Implement XXHEAVISIDE: atom-level cutoff truncation in innerloop
  - [x] Cutoff comparison: gromosXX uses `d > cutoff²` (strictly greater, pairs AT cutoff included)
  - [x] Wire RCUTP as short-range cutoff (twin-range)
  - [x] Split RF excluded: solvent gets only DD term (no self, no DI, no forces)
  - [ ] **Implement `put_chargegroups_into_box`** — gromosXX gathers CGs into primary image before pairlist
    - Pairlist pair counts still off: solute_short=12+long=24, solvent_short=70+long=19
    - Box is [0, 2) nm; some CGs may need wrapping for correct distances
    - This is the most likely remaining cause of the ~3.82 kJ/mol discrepancy

### Later — Level 3 (bulk + ensembles)
- [ ] water_216_box NVE: bulk nonbonded scaling, pairlist updates, virial
- [ ] water_216_nvt: Berendsen thermostat integration
- [ ] water_216_npt: Berendsen barostat + pressure coupling + box rescaling

### Later — Level 4 (complex)
- [ ] aladip_vacuum: all bonded + exclusions + 1-4 interactions
- [ ] aladip_solvated: SHAKE + mixed solute-solvent system
- [ ] Triclinic boundary condition wiring in md.rs

### Known Gaps
- Triclinic box: code exists in math.rs but md.rs never creates Triclinic periodicity
- 1-4 interactions: parsed (INE14) but need verification in force loop
- Virial / pressure calculation: needs verification for NPT
- COM motion removal (COMTRANSROT): parsed but not wired
- `put_chargegroups_into_box`: not implemented — may affect large systems with atoms outside [0, L)

## Investigation: nacl_water_box CRF mismatch

**Symptom:** E_pot off by ~3.82 kJ/mol at step 0. E_lj matches exactly.
Our step 0: E_pot=-112.570, E_kin=0.0138, E_crf=-103.572
Reference:  E_pot=-116.390, E_kin=0.0151, E_crf=-107.392
(E_kin now correct after SHAKE fix; E_crf difference is ~3.82 kJ/mol)

### gromosXX Source Code Analysis (from Seafile/GROMOS/gromos_git/gromosXX/gromosXX/src/)

**Source files read:**
- `interaction/nonbonded/pairlist/standard_pairlist_algorithm.cc` (640 lines)
- `interaction/nonbonded/interaction/nonbonded_term.cc` (719 lines)
- `interaction/nonbonded/interaction/nonbonded_innerloop.cc` (1987 lines)
- `interaction/nonbonded/interaction/solvent_innerloop.cc` (1721 lines)
- `interaction/nonbonded/interaction/nonbonded_set.cc`
- `interaction/nonbonded/interaction/nonbonded_outerloop.cc`
+
#### 1. Twin-range pairlist: TWO separate cutoffs
gromosXX uses `m_cutoff_short_2` (RCUTP²) and `m_cutoff_long_2` (RCUTL²):
- `d > m_cutoff_long_2` → SKIP (outside)
- `d > m_cutoff_short_2` → LONGRANGE list (interactions recalculated every NSNB steps)
- `d <= m_cutoff_short_2` → SHORTRANGE list (interactions recalculated every step)

**Our bug:** We set `short_range_cutoff = long_range_cutoff = RCUTL (0.9)`.
The nacl_water_box input has `RCUTP=0.8, RCUTL=0.9`. We ignore RCUTP entirely!
This means we're NOT using twin-range — all pairs within 0.9 go to short-range,
none to long-range. This shouldn't cause a step-0 energy error (same pairs are
evaluated), but it will cause divergence on subsequent steps.

#### 2. `put_chargegroups_into_box` — called BEFORE pairlist
In `prepare()` (line ~40-42): `periodicity.put_chargegroups_into_box(conf, topo);`
This wraps all CG positions into the primary box image [0, L) BEFORE computing
distances. Then solute CG COGs are calculated from the wrapped positions.

#### 3. Solvent-solvent pairlist stores ALL atom pairs
gromosXX `_solvent_solvent()` (line ~320-365) stores individual atom pairs:
```cpp
for(a1 = cg(cg1); a1 < cg(cg1+1); ++a1)
  for(a2 = cg(cg2); a2 < cg(cg2+1); ++a2)
    pairlist.solvent_short[a1].push_back(a2);
```
**Our code** stores only first-atom pairs and expands in `solvent_innerloop()`.
This is functionally equivalent as long as the innerloop expansion matches.

#### 4. HEAVISIDE truncation uses `cutoff_long²`
In `nonbonded_term.cc` (line ~183):
```cpp
m_cut2 = sim.param().pairlist.cutoff_long * sim.param().pairlist.cutoff_long;
#ifdef XXHEAVISIDE
#define HEAVISIDETRUNC(rlen, vars) if (rlen > m_cut2) { vars = 0.0; return; }
```
HEAVISIDE uses `cutoff_long²` (RCUTL²), NOT `cutoff_short²`.
**Our code** uses `crf.cutoff_sq` which is set from the same RCUTL. This matches.

#### 5. RF excluded interactions: solute vs solvent DIFFER
**Solute RF excluded** (`RF_excluded_interaction_innerloop`, line 1291):
- Self-term: `rf_interaction(r=0, qi*qi, f, e_crf)` → `e_self = 0.5 * qi² * FPEPSI * (-crf_cut)`
- Excluded pairs: `rf_interaction(r, qi*qj, f, e_crf)` → full RF correction with force + energy
- `rf_interaction` formula: `force = q * FPEPSI * crf_cut3i * r`, `energy = q * FPEPSI * (-crf_2cut3i * r² - crf_cut)`

**Solvent RF excluded** (`RF_solvent_interaction_innerloop`, line 1547):
- **NO self-term** for solvent atoms (comment: "distance independent parts should add up to zero")
- **NO forces** for excluded solvent pairs (rigid molecules → forces are zero)
- **Only energy:** `e_crf = -qi*qj * FPEPSI * crf_2cut3i * r²` (ONLY the distance-dependent part)
- The constant `-crf_cut` term is deliberately OMITTED for solvent

**THIS IS A KEY DIFFERENCE FROM OUR CODE!**
Our `rf_excluded_interactions()` applies the same formula to ALL atoms (solute + solvent):
- Same self-term for all atoms: `-0.5 * qi² * FPEPSI * crf_cut`
- Same excluded-pair energy: `qi*qj * FPEPSI * (-crf_2cut3i * r² - crf_cut)`
- Same forces for all excluded pairs

In gromosXX, solvent excluded pairs get only `-qi*qj * FPEPSI * crf_2cut3i * r²`
(no `-crf_cut` constant, no force, no self-term).

The reasoning: for neutral solvent molecules (charge group sums to zero),
the self-terms and `-crf_cut` constants cancel exactly across all atoms in the
molecule. gromosXX skips them for efficiency. But the distance-dependent
`-crf_2cut3i * r²` part still matters for the energy.

#### 6. `lj_crf_interaction` formula comparison
gromosXX (nonbonded_term.cc line ~193-210):
```cpp
e_lj = (c12_dist6i - c6) * dist6i;
e_crf = q_eps * (disti - crf_2cut3i * dist2 - crf_cut);
force = (c12_dist6i + c12_dist6i - c6) * 6.0 * dist6i * dist2i
      + q_eps * (disti * dist2i + crf_cut3i);
```
**Our code matches this exactly.** ✓

#### 7. Solvent innerloop: shared PBC shift
gromosXX `solvent_innerloop()` (solvent_innerloop.cc line ~34-38):
```cpp
periodicity.nearest_image(*pos_i, *pos_j, r);
tx = r(0) - (*pos_i)(0) + (*pos_j)(0);
ty = r(1) - (*pos_i)(1) + (*pos_j)(1);
tz = r(2) - (*pos_i)(2) + (*pos_j)(2);
```
Then for each atom pair: `x = (pos_i+atom_i)(0) + tx - (pos_j+atom_j)(0)`
**Our code matches this.** ✓

### Current Debug Output (step 0)
```
Pairlist: 22 CGs total, 2 solute CGs, cutoff_short²=0.640000, cutoff_long²=0.810000
solute_short=12, solute_long=24, solvent_short=70, solvent_long=19
RF solute: self=-2.2968e2, excl=2.0643e2 (1 pairs)
RF solvent: excl=-4.1676e0 (60 pairs)
Bond: 0  LJ: -8.998  CRF: -103.572
```

Note: twin-range now active (RCUTP=0.8, RCUTL=0.9). RF excluded correctly split:
solvent gets only the DD term. The ~3.82 kJ/mol gap is in the pairlist CRF (fewer
pairs evaluated due to missing `put_chargegroups_into_box`).

### Diagnostic Systems for Bug Isolation

Created 4 minimal systems to bisect the nacl_water_box CRF mismatch:

| System | Atoms | What it tests | Notes |
|--------|-------|---------------|-------|
| nacl_pair_box | 2 | Ions in PBC box with RF, no solvent (NSM=0) | Isolates: does PBC + RF work for solute-only? |
| nacl_1water_box | 5 | Na+Cl- + 1 SPC water (NSM=1) | Isolates: minimal solvent expansion + SHAKE |
| nacl_3water_box | 11 | Na+Cl- + 3 SPC waters (NSM=3) | Isolates: solvent-solvent + solute-solvent pairlist |
| water_3_box_twinrange | 9 | 3 waters with RCUTP=0.4, RCUTL=0.8, NSNB=5 | Isolates: twin-range pairlist mechanics |

All use the nacl_water_box topology (SOLVENTATOM block) or water_3_box topology.
Reference data generated with gromosXX md++ (Release, no OMP/MPI).


## GROMOS Book Vol2 — Key Conventions and Their Relation to the Bugs

### 1. Reaction-Field Electrostatics (Lines 6306–6420)

The RF electrostatic energy is decomposed into three terms (vol2.tex):

$$V^{\text{ele,pair}}_{\text{RF}} = V^{\text{CB}} + V^{\text{DD}} + V^{\text{DI}}$$

- **Coulombic (CB):** $\frac{q_i q_j}{4\pi\epsilon_0 r_{ij}}$
- **Distance-Dependent (DD):** $-\frac{q_i q_j C_{rf}}{4\pi\epsilon_0 \cdot 2R_{rf}^3} r_{ij}^2$
- **Distance-Independent (DI):** $-\frac{q_i q_j (1 - C_{rf}/2)}{4\pi\epsilon_0 R_{rf}}$

The critical passage is at **line 6391–6394**:

> In GROMOS, [the Coulombic term] is **not evaluated for excluded atoms**, while [the DD and DI terms] **are evaluated for these atoms as well**, unless the simulation is performed in the GROMOS96 compatibility mode.

And the **self-term** (line 6403–6406):

$$\phi^{\text{self}} = -\frac{1 - C_{rf}/2}{R_{rf}}$$

This contributes $\frac{1}{8\pi\epsilon_0} \sum_i q_i^2 \cdot \phi^{\text{self}}$ to the total energy.

**What this means for your bug:** The book says excluded pairs get DD + DI, and all atoms get the self-term. However, gromosXX *optimizes* the solvent case by recognizing that for neutral charge groups, the DI constants and self-terms cancel exactly within each molecule. So gromosXX skips them for solvent, computing only the DD piece (`-qi*qj*FPEPSI*crf_2cut3i*r²`). Your code applies the full formula (self + DD + DI) to ALL atoms including solvent. The extra constant terms don't cancel numerically because you're summing them separately—they produce a non-zero residual that accounts for the ~4.35 kJ/mol CRF mismatch.

---

### 2. Charge Group Positions: Solute vs Solvent (Lines 1319–1335)

The book explicitly states two different conventions:

- **Solute CG position** (line 1324): Centre of geometry: $R_{cg} = \sum_{i=1}^{N_{cg}} \mathbf{r}_i / N_{cg}$
- **Solvent CG position** (line 1333): Position of the **first atom** of the solvent molecule. A solvent molecule may only contain one charge group.

**What this means:** Your pairlist code must use different CG position calculations for solute vs solvent. This is already noted in your PLAN as implemented (`[x] Solvent CG distance: use first-atom position`).

---

### 3. Twin-Range Pairlist (Lines 1353–1415)

The book defines (line 1368–1377):

- **Short-range cutoff (RCUTP):** Used for the pairlist. Interactions within [0, RCUTP] recalculated every step.
- **Long-range cutoff (RCUTL):** Interactions in [RCUTP, RCUTL] evaluated every NSNB steps and held constant between updates.

**What this means:** ~~Your code currently sets both cutoffs to RCUTL, ignoring RCUTP.~~ **FIXED:** RCUTP is now wired as the short-range cutoff, and long-range forces are cached and reused between pairlist updates.

---

### 4. Heaviside Cutoff Truncation (Lines 6322–6341)

Two modes defined at lines 6322 and 6332:

- **AT (atom-based):** $H(R_{cut} - r_{ij})$ — truncation based on atom-atom distance
- **CHG (charge-group-based):** $H(R_{cut} - r_{cg,ij})$ — truncation based on CG-CG distance, but the interaction function uses the actual atom-atom distance $r_{ij}$

The Heaviside convention (line 6341): $H(\xi) = 1$ when $\xi \geq 0$ (i.e., pairs **at** exactly the cutoff ARE included). This matches gromosXX's `d > cutoff²` check (strictly greater means equal is included).

---

### 5. Gathering / Put-into-Box (Lines 2866–2884, Algorithm Step 1B at Line 18741)

The book states (line 2870):
> Solute charge group atoms and solvent molecules are translated, applying periodic boundary conditions such that the **first atom** of a solute charge group or of a solvent molecule lies within the central periodic box.

The algorithm (step 1B, line 18741):
> If required, apply the periodic boundary conditions to put the solute charge groups and solvent molecules in the central computational box.

This happens **before** force calculation (step 2). This is `put_chargegroups_into_box` in gromosXX—it gathers all CGs into the primary image before the pairlist is built, ensuring consistent distance calculations.

**What this means:** If your initial coordinates have atoms outside [0, L), the pairlist distances will be wrong without this gathering step. For small systems with atoms already in the box, it doesn't matter. For `nacl_water_box` (20 waters), some atoms might be outside.

---

### 6. SHAKE and Velocity Correction (Lines 15298–15350)

The book (lines 15311–15340) explains constrained velocities in the leap-frog scheme:

1. Compute unconstrained position: $\mathbf{r}_i^{uc}(t+\Delta t) = \mathbf{r}_i(t) + \mathbf{v}_i(t+\Delta t/2)\Delta t$
2. SHAKE: $\text{SHAKE}(\mathbf{r}(t); \mathbf{r}^{uc}(t+\Delta t); \mathbf{r}(t+\Delta t))$
3. Constrained velocity: $\mathbf{v}_i(t+\Delta t/2) = [\mathbf{r}_i(t+\Delta t) - \mathbf{r}_i(t)] / \Delta t$

**What this means for nacl_3water_box:** The velocity is simply the position displacement divided by dt. If SHAKE corrects positions for all solvent molecules but the velocity recalculation (step 8C in the algorithm, line 18874) is only applied to the first molecule (or DOF counting is wrong), you get the 30x E_kin discrepancy. The book (line 15532–15537) shows that for a 3-atom rigid water:

$$N_{dof}(\text{solvent}) = 3N_{atoms} - N_{constraints} = 9 - 3 = 6$$

per water molecule. If you're dividing by wrong DOF or not applying the velocity correction to all molecules, E_kin blows up.

---

### 7. Energy Reporting — Kinetic Energy Averaging (Line 18940)

The algorithm step 13 (line 18940):
> The kinetic energy at $t_n$ is calculated as the **average** of the kinetic energies at $t_n - \Delta t/2$ and $t_n + \Delta t/2$.

This is a key GROMOS convention — the reported E_kin is the average of old and new half-step kinetic energies.

---

### Summary of Root Causes

| Bug | Book Reference | Status | Explanation |
|-----|---------------|--------|-------------|
| **nacl_3water_box E_kin 30x too large** | Lines 15311–15340, 15532–15537 | **FIXED** | SHAKE only looped over `topo.solute.bonds`, never applied constraints to solvent molecules. Added solvent constraint loop over all molecules. |
| **Twin-range drift after step 0** | Lines 1368–1377 | **FIXED** | RCUTP wired as short-range cutoff, long-range forces cached on NSNB steps and reused on intermediate steps. |
| **nacl_water_box CRF ~3.82 kJ/mol off** | Lines 6391–6394, 2866–2884 | OPEN | RF excluded split done (solvent: DD only). RCUTP wired. Remaining: `put_chargegroups_into_box` not implemented — pair counts still off. |
| **Pairlist pair count mismatch** | Lines 2866–2884, step 1B | OPEN | `put_chargegroups_into_box` gathers CGs into primary image before pairlist. Without it, distance calculations may differ. |


### Action Items
1. ~~**Fix RF excluded for solvent:**~~ ✓ Split into solute/solvent paths. Solvent: no self, no DI, no forces, only DD.
2. ~~**Wire RCUTP:**~~ ✓ RCUTP used as short-range cutoff. Twin-range caching implemented.
3. **Implement `put_chargegroups_into_box`:** Wrap CG positions into primary box before pairlist update.
   This is the remaining fix needed for nacl_water_box.
   May affect pair counts when atoms are outside [0, L).
4. **Run diagnostic systems:** Use nacl_pair_box → nacl_1water_box → nacl_3water_box
   progression to pinpoint where CRF diverges. If nacl_pair_box passes, the bug is
   in solvent RF handling. If it fails, the bug is in PBC+RF for ions.

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

# Run integration tests against gromosXX references (10 active, 7 ignored)
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