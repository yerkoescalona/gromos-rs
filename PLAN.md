# gromos-rs — Roadmap & Reference Tests

Focus: `cargo build --release --bin md`
On commit: update CHANGELOG.md and Cargo.toml version.
DON'T modify `gromosXX_references/*/expected/` — those are ground truth.

References:
- gromosXX source (MD engine, "md++"): `.local/gromosXX/md++/src`
- gromosPlsPls source (analysis/tools, "gromos++"): `.local/gromosPlsPls/gromos++/src`
- Tutorials: `.local/gromos_tutorial_livecoms/tutorial_files`
- Theory: `.local/doc/gromos_book`
- Force fields: `.local/gromosXX/forcefields`
- **`FUTURE.md`** — architectural bets (SoA core, O(N) pairlist, QM/MM + ML potentials, Martini
  bridge, unifying layered architecture) + differential-audit findings (known gromosXX bugs not to
  port, live divergences). PLAN.md = near-term execution; FUTURE.md = where we diverge on purpose.

**Per-crate status, key files, and crate-specific rules live in each crate's stage contract:**
`crates/<crate>/.claude/CONTEXT.md` — read it before touching that crate.

---

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
| 2   | nacl_1water_distres | 5  | distance restraint on Na-Cl pair (NTDIR=2, CDIR*w0) | **PASS** |
| 4   | ch4_water_fep | 2998 | CH4→dummy in 999 SPC water, λ=0.5, twin-range NB FEP | **PASS** |

**36 of 36 tests pass.** All levels fully passing (1 ignored: aladip_vacuum_fep, pending Step 4).

(No reference tests yet for `gromos-analysis` / `gromos-tools` — see P2 + cross-cutting below.)

---

## Roadmap (priority order)

Overarching principle: **gromosXX as the reference, no duplication between the engine (gromosXX) and
the analysis/tools facade (gromosPlsPls).** Every feature lands with a minimal gromosXX reference test.

### Priority 1 — MD engine physics (gromosXX-faithful, reference-tested)
Wire the already-coded-but-unwired physics; keep implementations in `gromos-forces`/`gromos-integrators`
(reusable), never duplicated into binaries. Each item gets a minimal reference test.

**1.1 — Reproducibility & correctness** ✓ complete
- [x] NTIVEL=1 velocity generation (Maxwell-Boltzmann) — `water_single_genvel` passes
- [x] Unit-conversion audit (topology parsing) — all conversions verified, no bugs found

**1.2 — Constraints** ✓ complete
- [x] SETTLE — `nacl_1water_settle` passes
- [x] LINCS — `nacl_1water_lincs`, `aladip_vacuum_lincs` pass
- [x] COM rotation removal — `water_216_box_com_rot` passes

**1.3 — Thermostat** ✓ complete
- [x] Nosé-Hoover single NHC + chain NHC — `water_216_nvt_nosehoover`, `water_216_nvt_nhc_chain` pass

**1.4 — Boundary** ✓ complete
- [x] Triclinic nearest-image: replaced fractional `round()` with gromosXX while-loop z→y→x reduction
- [x] NTB=-1 truncated-octahedron: `truncoct_triclinic_box` + position/velocity rotation + output rotation — `aladip_trunc_oct` passes

**1.5 — O(N) cell-list pairlist** ✓ complete (not yet wired into md binary)
- [x] `CellListPairlistAlgorithm` — bins chargegroup COGs; validated by set-equality vs `StandardPairlistAlgorithm`; O(N) for rectangular boxes, falls back to O(N²) for triclinic/vacuum
- [x] Martina bug NOT reproduced; solute/solvent classification by both atoms' roles
- [ ] Deferred: system-size heuristic to wire CellList into md.rs; spatial reorder; displacement-triggered rebuild; triclinic acceleration

**1.6 — Restraints & special interactions** ✓ distance restraints wired; rest parked
- [x] Distance restraints — `nacl_1water_distres` passes (NTDIR=2, CDIR*w0, instantaneous).
  Physics in `gromos-forces/restraints.rs`, wired into `Forcefield`. `Energies::total()` fixed
  to include `distanceres_total` (mirrors gromosXX `special_total`). Perturbed variant
  unit-tested vs aladip reference values (257.19 / 195.90 kJ/mol).
- [ ] Dihedral restraints — parked; return after P1.7. Note gromosXX caveat at `:152`
  re phi0_A/phi0_B > 2π apart — second-source before porting.
- [ ] Angle restraints — parked
- [ ] J-value restraints — parked
- [ ] Order-parameter restraints — parked; two already-fixed bugs at `:91,:213` in gromosXX
- [ ] Distance-field + local elevation — parked
- [ ] RDC, X-ray, symmetry, colvar, electric-field, NEMD — parked

**1.7 — FEP / TI** ← **NEXT PRIORITY**

Goal: run a full TI simulation of a simple alchemical transformation and produce a `.trg`
free-energy trajectory that `ext_ti_ana` can integrate into ΔG.

Virtual atoms assessment: NOT required for core FEP/TI alchemical transformations. Standard
GROMOS FEP uses explicit atoms with λ-scaled LJ/CRF. Port virtual atoms only when a concrete
use case requires them (e.g., united-atom CH groups in a perturbation, TIP4P water).

Steps (each with a reference test before wiring):

- [x] **Step 1 — `.pttopo` reader** (`gromos-io/ptp.rs`) ✓
  - `read_pttopo()` parses PERTATOMPARAM, PERTATOMPAIR, PERTBONDSTRETCH(H),
    PERTBONDANGLE(H), PERTIMPROPERDIH(H), PERTPROPERDIH(H); all 1→0 index conversions fixed
  - Perturbed terms removed from `topo.solute` in `md.rs` to avoid double-counting
  - `lambda_and_derivative` wired into Forcefield; `dhdl_total` added to `Energies`
  - Soft-variant blocks (PERTBONDSOFT, PERTANGLESOFT, PERTIMPROPERDIHSOFT) log a warning;
    the atoms they refer to are NOT yet removed from `topo.solute` — deferred to Step 3

- [x] **Step 2 — Perturbed bonded forces** (`gromos-forces/bonded/perturbed.rs`) ✓
  - All four terms rewritten faithfully to gromosXX sources:
    quartic bond `¼K(r²−r0²)²`, cos-harmonic angle `½K(cosθ−cos0)²`,
    improper `½K(ζ−ζ0)²`, proper dihedral (states A+B computed separately, then combined)
  - dH/dλ computed per term; accumulated in `conf.energies.dhdl_total`
  - `calculate_perturbed_bonded_forces` wired into Forcefield step 2b
  - `aladip_vacuum_fep` reference exists (ignored pending Step 3): bonded matches
    gromosXX to <0.1 kJ/mol; nonbonded off by ~13 kJ/mol — perturbed atoms still
    use state-A charges in the pairlist until Step 3 implements λ-scaled charges

- [x] **Step 3 — Perturbed nonbonded** (`gromos-forces/nonbonded/perturbed.rs`) ✓
  - Soft-core LJ+CRF dual-topology interaction implemented (`perturbed_lj_crf_interaction`)
  - Correction approach: state-A innerloop runs on all pairs; `perturbed_pairlist_correction`
    subtracts state-A and adds full soft-core perturbed for each pair where `pert[i].is_some()`
  - Self-energy, excluded-pair, 1-4, and PERTATOMPAIR corrections all wired
  - `ch4_water_fep` reference test (CH4→dummy in 999 SPC water, λ=0.5, twin-range,
    RCUTP=0.8/RCUTL=1.4) passes to <1e-6 kJ/mol vs gromosXX
  - **Bugs fixed during Step 3 — read before touching this area:**
    1. **PERTURBATION block multi-line** (`gromos-io/src/imd.rs`): the block splits NTG/NRDGL/RLAM/DLAMT on line 1 and ALPHLJ/ALPHC/NLAM/NSCALE on line 2; the parser was calling `data_lines.first()` and silently dropped ALPHLJ=0. Fixed by combining all data lines.
    2. **Effective alpha = per-atom × global** (`forcefield.rs:build_pert_info`): the `.ptp` ALJ/ACRF column is a scaling factor; gromosXX multiplies it by the global ALPHLJ/ALPHC (`in_perturbation.cc:1308`). `build_pert_info` now does `alpha_lj = pa.lj_soft * self.global_alphlj`. Without this, soft-core distances were wrong and close-range repulsive pairs gave ~2× too large a correction.
    3. Bug 1 masked Bug 2: ALPHLJ=0 from the broken parser zeroed the whole effective alpha, producing a different wrong answer. Both fixes are needed together.
  - **Debugging methodology that worked:** (a) verify NTG=0 agrees first to confirm topology/innerloop are correct; (b) add per-pair `eprintln!` in both codes for 5 pairs (one close, one mid, one long-range); (c) if per-pair match but totals differ → pairlist structure; (d) if per-pair differ → parameter mismatch → back-calculate what alpha gromosXX must be using from its output energy for one pair, then trace that value in gromosXX C++ source.
  - ⚠️ **Perturbed RF self-term** still needs second-sourcing from GROMOS book (gromosXX authors flagged at `perturbed_nonbonded_term.cc:596,749,1444`). CRF in `ch4_water_fep` is zero-charge so this wasn't exercised.

- [x] **Step 4 — ∂V/∂λ accumulation + output** ✓
  - `dhdl_total` accumulated in `Energies` from bonded (`calculate_perturbed_bonded_forces`) + NB (`pert_nb_dhdl`) each step
  - `FreeEnergyWriter` in `gromos-io/src/free_energy.rs`: writes FREEENERGY03 blocks to `.trg`
  - `md` binary wires `@trg` output when NTG != 0; writes at `NTWE` frequency
  - Reference test `ch4_water_fep` verified: dH/dλ at λ=0.5 (CH4→dummy) matches self-consistent baseline to 1e-6 rel; `expected/free_energy.trg` committed
  - Shared test data moved to `tests/gromosXX_references/shared/`; test runner tracks `.trg` via `free_energy` key in `input.toml`

- [x] **Step 5 — `ext_ti_ana` tool** ✓
  - `gromos-core/src/stat.rs`: `Stat` with `ave()`, `rmsd()`, `msd()`, `ee()` (Allen–Tildesley block averaging, faithful port of gromos++ `gmath/Stat`)
  - `gromos-io/src/free_energy.rs`: `read_free_energy_trajectory()` reader for FREEENERGY03 blocks
  - `ext_ti_ana`: reads N `.trg` files (@trg), optional @lambda override, @skip equilibration, trapezoidal ΔG, per-window ⟨dH/dλ⟩ ± ee()

**1.8 — Virtual atoms** (deferred; not blocking FEP/TI for standard use cases)
- [ ] Port `algorithm/virtualatoms/` (aromatic centroids, lone pairs, TIP4P site)
- [ ] Needed for: united-atom NMR restraints, TIP4P water model, some perturbed topologies
- [ ] Coordinate with Dim 10 instancing refactor (FUTURE.md) so virtual sites aren't
  hard-coded around the old solute/solvent split

**1.9 — Advanced sampling** (stubs exist; delegatable after FEP/TI is solid)
- [ ] **EDS** — V_mixed = −1/β·ln(Σ exp(−β(Eᵢ−eir_i))); per-state force eval + blending; AEDS.
- [ ] **GaMD** — V_boost = k·(V−E_threshold)²; Welford running stats; dihedral/total/dual.
- [ ] **REMD** — MPI parallel tempering; Δ = (β₁−β₂)(E₁−E₂); feature-gated MPI.

---

## Dimension 10 — Instancing model (FUTURE.md Dim 10) — IN PROGRESS

> Dissolve the solute/solvent split: separate *representation* (store once, instance N times)
> from *role* (per-instance attribute). See FUTURE.md §Dim10 for full design rationale.

### Why `MoleculeType.bonds` is correct (not scope creep)

`role` and `MoleculeType.bonds` are **independent knobs** — this is what makes the design
powerful:

- **`role`** controls: thermostat/pressure group, `s:`/`m:` selection, analysis membership,
  fast-path dispatch (e.g. SETTLE for rigid water vs SHAKE for flexible)
- **`MoleculeType.bonds`** controls: whether bonded *forces* are computed for a molecule type
  (`empty` = geometry kept by constraints only; `non-empty` = explicit bond/angle forces)

This handles every combination without special-casing:

| Simulation | `mt.bonds` | `role` |
|------------|-----------|--------|
| Standard GROMOS: flexible solute + rigid SPC | solute: full / water: empty | Solute / Solvent |
| Flexible water (SPC/Fw): both have bond forces | both: full | Solute / Solvent |
| DMSO solvent: full bonding, grouped as solvent | full | Solvent |
| Promoted catalytic water: solute for restraints, rigid geometry | empty | **Solute** after `promote()` |
| Mixed: some waters flexible, most rigid | two water MoleculeTypes | Solvent |

**The current gromosXX limitation** (DMSO must be faked as "solute atoms") disappears: any
molecule can have full bonded forces regardless of its role.

**The final bonded force loop (Phase 3 target):**
```rust
for inst in &topo.instances {
    let mt = &topo.moltypes[inst.moltype_id];
    for bond in &mt.bonds {
        let i = inst.atom_offset + bond.i;  // local → global
        let j = inst.atom_offset + bond.j;
        // compute bond force — one loop, no solute/solvent special-casing
    }
}
```

### Phase status

**Phase 1 — COMPLETE ✓ (37/37 reference tests pass, byte-identical output):**
- [x] `Role`, `MolTypeAtom`, `MoleculeType`, `MoleculeInstance` types
- [x] `moltypes` + `instances` fields on `Topology`
- [x] `init_solute_moltype()`, `solvate()` populates instances
- [x] Per-atom accessors (`atom_name`, `residue_nr`, `residue_name`) prefer moltype path
- [x] `role_of_atom(i)`, `is_solvent_atom(i)`, `promote(mol_idx)` API
- [x] `s:` in `AtomSelection` uses `role == Role::Solvent`

**Phase 2a ✓** — `rebuild_flat_arrays()`: flat `iac`/`mass`/`charge` derived from instances.

**Phase 2b ✓** — `MoleculeType` carries all bonded terms. Read accessors `solute_bonds()` etc.
route through `moltypes[0]`. ~100 read-access call sites migrated. (Uncommitted diff — pending review.)

**Phase 2c ✓** — `Solute` now contains only `atoms: Vec<Atom>`. All bonded write sites use
`mut_solute_bonds()` / `mut_solute_angles()` / etc. (Uncommitted diff — pending review.)

**Phase 2d — DEFERRED:** Remove `solute.atoms: Vec<Atom>`.
  - `Atom` has boolean flags (`is_perturbed`, `is_coarse_grained`, `is_polarisable`) absent from
    `MolTypeAtom` — merge the two types first
  - ~36 read + ~15 write call sites; topology writer iterates with full `Atom` struct
  - Unblock by: add flags to `MolTypeAtom`, add `mut_solute_atoms()`, migrate iteratively

**Phase 2e — DEFERRED:** Remove `Vec<Solvent>` (12 accesses in constraints/settle/md.rs);
replace with moltype lookup via instances + constraint template.

**Phase 2f — DEFERRED:** `promote(mol_idx)` with CG/exclusion renumbering (needs Dim 9d).

**Phase 3 — THE GOAL:** Replace per-term bonded force loops in `bonds.rs`, `angles.rs`,
`dihedrals.rs` with the instance-iterating loop above. This is when the DMSO/flexible-solvent
use case becomes live and the `solute_bonds()` bridge accessors can be removed entirely.

**Current state:** `Solute` has only `atoms`; all bonded topology is in `MoleculeType`; flat
arrays are derived caches; `s:` uses role. Remaining duality: `solute.atoms` + `Vec<Solvent>`.
Phase 2d+2e remove those; Phase 3 wires the instance loop into the force calculations.

---

### Priority 2 — Analysis foundations (the no-duplication layer)

> **Architectural note (FUTURE.md Dim 10):** do selection/gathering work in a way that treats
> solvent-ness as an attribute, not a partition — don't calcify around the old split.

**2.1 — Atom selection** (gromosXX primitives first, gromos++ facade on top)
- [x] Solidify queryable per-atom metadata — `Topology::atom_name(i)`, `residue_nr(i)`, `residue_name(i)`, `molecule_nr(i)` cover all atoms (solute + solvent) uniformly without `num_solute_atoms()` threshold ✓
- [x] Fix known AtomSelection gaps + full gromos++ grammar — `a:name`, `1:name,name`, `1:res(nr:atom)`, `1:res(name:atom)`, `not(spec)`, `minus(spec)`, `;`-union, `all`/`no`; all routes use `topology.molecules` (no `num_solute_atoms()` threshold, Dim 10 ready) ✓
- [x] `atominfo` binary (`gromos-analysis`) — reads topology + AtomSpecifier, prints TITLE+ATOMS block; output verified against gromos++ `atominfo` on aladip (12 atoms, all selection forms) ✓
- [x] 30 reference tests in `selection.rs`, every expected index confirmed by `gromos++ atominfo` ✓
- [ ] AtomSpecifier grammar facade (deferred — `s:`/`m:` unification pending Dim 10 role-attribute model)

**2.2 — Shared analysis infra into lower crates** (kills duplication)
- [ ] Rotational fit (Kabsch/SVD) — ref: `.local/gromosPlsPls/gromos++/src/fit/RotationalFit.cc`
- [x] Statistics + error-estimate — `gromos-core/src/stat.rs`: `Stat` with `ave()`, `rmsd()`, `ee()` (block averaging, gromos++ faithful) ✓
- [ ] PBC gathering / molecule unwrapping — one primitive in gromos-core consumed by both engine and gromos++-style tools
- [ ] Single-point energy entry point so `ener.rs` calls gromos-forces

**2.3 — Real program implementations** (once 2.1/2.2 land)
- [ ] `rmsd` — real rotational fit (Kabsch)
- [x] `ext_ti_ana` — reads N `.trg` files, integrates ⟨dH/dλ⟩ over λ (trapezoidal), reports ΔG ± ee() ✓
- [ ] `nhoparam` — rewrite to actual gromos++ algorithm (NMR N-H order parameters S²)

**2.3b — Free-energy analysis estimators** (once ext_ti_ana + ee() land)
- [ ] `bar` (Bennett Acceptance Ratio), `ext_ti_merge`, `reweight`, `m_widom`, `dg_ener`

**2.4 — Stub cleanup**
- [ ] `visco`, `frameout`, `amber2gromos`, `ener`, `sasa_hasel`, `dssp`, `solute_entropy`

### Priority 3 — py-gromos API & education
Design reference: `.local/polars` (Python API, method chaining, pyo3 patterns).

> **Design target (FUTURE.md "Compositional topology"):** `ForceField` → `BuildingBlock` → `Topology`
> as an algebra (`+`/`*`/`solvate()`/`check()`); Python expresses verbs; Rust core owns every
> invariant. Build → minimize → simulate in one address space, zero files.

- [ ] Phase 1 — Rust bindings (pyo3-gromos): see pyo3-gromos CONTEXT.md for remaining items
- [ ] Phase 2 — Python API (py-gromos): method chaining, energy DataFrame, rich reprs
- [ ] Phase 3 — Notebooks & education: rewrite `py-gromos/notebooks/` + `examples/`

### Priority 4 — Code quality (last)
- [ ] Clippy (~390 warnings: gromos-forces 89, gromos-integrators 77, gromos-io 31, gromos-core 15)
- [ ] Replace bare `unwrap()` in non-test code with `.expect("msg")` or `?`
- [ ] Add missing `#[test]`: constraints (SHAKE — 0 today), improper dihedral
- [ ] Split large files: `nonbonded.rs` (~1500 LOC), `bonded.rs` (~1300 LOC), `gromos-io/topology.rs` (~1200 LOC)
- [ ] Unify CLI error types; audit `pub` visibility
- [ ] Benchmarking infra: baseline `cargo bench --workspace -- --save-baseline v0.1`; add end-to-end MD step / pairlist / SHAKE / bonded benches; document in CONTRIBUTING.md

### Cross-cutting — minimal reference tests (do continuously)
Full tutorial t_01–t_06 end-to-end runs are **deferred to last** (compute-limited). Substitute:
hand-craft **minimal** reference systems and diff gromos-rs output against the C++ reference.
Stand up an analogous harness for analysis/tools (none exists today). Every P1 physics feature
and every P2 program lands with a minimal reference test.

**Free win — mine gromosXX's own `check/*.t.cc` regression suite.** The gromosXX devs hard-code
per-term reference energies in `md++/src/check/`: `aladip.t.cc` carries `QuarticBond=18.053811`,
`NonBonded_newRF=-84.092443`, `DistanceRestraint=257.189539`, and the perturbed/soft-core terms.
Also `c16_cg.t.cc`, `lambdas.t.cc`, `scaling.t.cc`, `aladip_ls.t.cc` (lattice-sum). Porting
these as unit tests gives per-term validation independent of running the md binary — and is a
genuine second source of truth for exactly the perturbed terms porting.md says not to trust the
C++ on. High value, low cost; do it alongside the FEP/restraints work.

### Cross-cutting — differential audit (do continuously) — FUTURE.md Dim 11
The reference suite is a bug oracle **only for wired paths.** Rules applied to every port:
1. **Reference test BEFORE wiring**, not after.
2. **Grep the C++ for self-flagged defects:** `grep -rniE 'bug|fixme|wrong|hack' interaction/ algorithm/ math/`
3. **Second-source uncertain physics** (RF self-terms, virial, Ewald): derive from the GROMOS book; diff against C++; investigate disagreements.
4. **Reproduce genuine GROMOS quirks as named, documented decisions.**
5. Eventually a `--gromos-compat` vs `--corrected` split so both modes can coexist.

Known items: triclinic nearest-image divergence (P1.4, resolved), Martina grid-pairlist bug
(P1.5, not ported), perturbed RF self-term uncertainty (P1.8), thermostat flexible-constraints
hack (P1.3, documented). Full table in FUTURE.md Dim 11.

### Parked / blocked
- [ ] `ext_ti_ana` — blocked behind P1.7 step 4 (`.trg` output) + P2.2 (`ee()` statistics)
- [ ] `nhoparam` — blocked behind P2 (selection + fit + stats + references)
- [ ] Remaining restraints (P1.6: dihedral, angle, J-value, order-parameter, distance-field,
  RDC, X-ray, colvar) — parked; resume after P1.7 FEP/TI is solid
- [ ] EDS, GaMD, REMD (P1.9) — stubs exist; delegatable to others after FEP/TI lands
- [ ] Tutorials t_01–t_06 end-to-end — compute-limited; replaced near-term by minimal reference tests

### Deferred breadth (tracked, not scheduled)
- [ ] **PME / lattice-sum electrostatics** — RF stays the focus by design. PME needs a focused investigation before committing (port gromosXX's lattice-sum or adopt a modern approach). Note the `// wrong!!!` traps in `interaction/nonbonded/interaction/latticesum.{h,cc}` (Dim 11).
- [ ] **Stochastic / Langevin dynamics** — `random_force` scaffolding exists; SD leap-frog + friction terms unported.
- [ ] **Coarse-grained → Martini bridge** — tracked as FUTURE.md Dim 13; gated on nonbonded-conventions investigation.
- [ ] **Polarisable / charge-on-spring force fields** — explicitly out of scope.

---

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
