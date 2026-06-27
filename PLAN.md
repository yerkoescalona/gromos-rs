# gromos-rs — Roadmap & Reference Tests

Focus: `cargo build --release --bin md`
On commit: update CHANGELOG.md and Cargo.toml version.
           Update the .claude files (in the root dir: .claude and within crates: crates/gromos-*/.claude).
DON'T modify `gromosXX_references/*/expected/` — those are ground truth.

References:
- gromosXX source (MD engine, "md++"): `.local/gromosXX/md++/src`
- gromosPlsPls source (analysis/tools, "gromos++"): `.local/gromosPlsPls/gromos++/src`
- Tutorials: `.local/gromos_tutorial_livecoms/tutorial_files`
- Theory: `.local/doc/gromos_book`
- Force fields: `.local/gromosXX/forcefields`
- **`FUTURE.md`** — architectural bets (SoA core, O(N) pairlist, QM/MM + ML potentials, Martini
  bridge, unifying layered architecture) + differential-audit findings (known GROMOS bugs not to
  port, live divergences). PLAN.md = near-term execution; FUTURE.md = where we diverge on purpose.

**Per-crate status, key files, and crate-specific rules live in each crate's stage contract:**
`crates/<crate>/.claude/CONTEXT.md` — read it before touching that crate.

---

## Reference Test Status

Validation: `cargo test -p gromos-md --test test_gromosXX_references`
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

**37 of 37 tests pass.** (1 ignored: aladip_vacuum_fep)

(No reference tests yet for `gromos-analysis` / `gromos-tools` — see P2 + cross-cutting below.)

---

## Roadmap (priority order)

Overarching principle: **GROMOS as the reference, no duplication between the MD engine and
the analysis/tools facade.** Every feature lands with a minimal GROMOS reference test.

### Priority 1 — MD engine physics (GROMOS-faithful, reference-tested)
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
- [x] Triclinic nearest-image: replaced fractional `round()` with GROMOS while-loop z→y→x reduction
- [x] NTB=-1 truncated-octahedron: `truncoct_triclinic_box` + position/velocity rotation — `aladip_trunc_oct` passes

**1.5 — O(N) cell-list pairlist** — algorithm complete, **wiring deferred** ← **NEXT (highest impact)**
- [x] `CellListPairlistAlgorithm` — bins chargegroup COGs; validated by set-equality vs `StandardPairlistAlgorithm`; O(N) for rectangular boxes, falls back to O(N²) for triclinic/vacuum
- [x] Martina bug NOT reproduced; solute/solvent classification by both atoms' roles
- [ ] **Wire CellList into md.rs** — add system-size heuristic (e.g. `n_atoms > 500` + rectangular box → CellList). This is the single biggest performance win: gromos-rs is currently O(N²) for all simulations, making 50k-atom production runs unusable. The algorithm is done, this is just the dispatch logic.
- [ ] Validate: 37/37 tests pass byte-identical after wiring
- [ ] Benchmark: confirm O(N) scaling on `water_216_box`

**1.6 — Restraints & special interactions** — distance done; dihedral next
- [x] Distance restraints — `nacl_1water_distres` passes (NTDIR=2, CDIR*w0, instantaneous).
  Physics in `gromos-forces/restraints.rs`, wired into `Forcefield`. Perturbed variant
  unit-tested vs aladip reference values (257.19 / 195.90 kJ/mol).
- [ ] **Dihedral restraints** — needed for practically every protein folding/conformational study. Second-source the `phi0_A/phi0_B > 2π` edge case (GROMOS bug at `:152`) before porting. Add reference test.
- [ ] Angle restraints — parked
- [ ] J-value, order-parameter, distance-field, local elevation, RDC, X-ray, colvar — parked

**1.7 — FEP / TI** ✓ complete
- [x] `.pttopo` reader, perturbed bonded forces, perturbed nonbonded (soft-core LJ+CRF)
- [x] dH/dλ accumulation, `.trg` output, `ext_ti_ana` integration
- [x] `ch4_water_fep` passes to <1e-6 kJ/mol vs GROMOS; `ch4_water_fep` tracks dH/dλ in reference test
- Note: perturbed RF self-term still needs second-sourcing from GROMOS book (flagged at `perturbed_nonbonded_term.cc:596,749,1444`). Zero-charge `ch4_water_fep` doesn't exercise it.

**1.8 — Virtual atoms** — skip for now; not blocking any common use case
- [ ] Port `algorithm/virtualatoms/` (aromatic centroids, lone pairs, TIP4P site)
- [ ] Only needed for: united-atom NMR restraints, TIP4P water, some perturbed topologies

**1.9 — Advanced sampling** — stubs exist; delegatable
- [ ] EDS — `V_mixed = −1/β·ln(Σ exp(−β(Eᵢ−eir_i)))`
- [ ] GaMD — `V_boost = k·(V−E_threshold)²`
- [ ] REMD — MPI parallel tempering; feature-gated MPI

---

### Priority 2 — Architecture + Analysis

> **Dim 10** (FUTURE.md): dissolve the solute/solvent split — representation (store once, instance N times)
> separate from role (per-instance attribute). Phases 1–2e complete ✓ (0.0.20):
> `moltypes[0]` = solute (bonds/atoms), `moltypes[1..]` = solvent types, `instances[k].role` = Solute/Solvent,
> `s:`/`m:` route through role, flat arrays derived from instances.

**2.0 — Dim 10 remaining**
- [x] Phases 1–2e: instancing model, role attribute, `moltypes[0].bonds` direct access, `Solvent` struct removed ✓
- [ ] **Phase 3 — bonded force loop replacement** — replace per-term loops in `bonds.rs`/`angles.rs`/`dihedrals.rs` with instance-iterating form (enables flexible DMSO solvent, closes architecture); validate 37/37 byte-identical
- [ ] `promote()` with CG/exclusion renumbering — needs Dim 9d charge-group primitives

**2.1 — Atom selection** ✓ complete
- [x] Full gromos-rs grammar: `a:name`, `1:name,name`, `1:res(nr:atom)`, `1:res(name:atom)`, `not()`, `minus()`, `;`-union
- [x] `m:`/`s:` route through `role` attribute (Dim 10). 32 reference tests, all confirmed by `atominfo`.
- [x] `atominfo` binary, better error messages with syntax hints

**2.2 — Shared analysis infra** ✓ complete
- [x] Kabsch rotational fit — `gromos-analysis/src/fit.rs` (Horn 1987 quaternion), 7 unit tests
- [x] Statistics + block-averaging `ee()` — `gromos-core/src/stat.rs`
- [x] PBC gathering — `gromos-core/src/gather.rs`: chain + bond-connectivity + molecule gathering
- [x] Single-point energy — `gromos-forces/src/energy.rs`: `single_point_energy()`, used by `ener` binary

**2.3 — Real program implementations** ✓ complete
- [x] `rmsd` — Kabsch fit, @atomspec, @ref, @pbc; all I/O through gromos-io
- [x] `nhoparam` — N-H order parameters S², rotational fit, window averaging, `ee()`
- [x] `ext_ti_ana` — trapezoidal ΔG ± ee(), reads `.trg` files
- [x] `frameout` — full GROMOS feature parity: PBC gather, @include SOLUTE/SOLVENT/ALL, rotational fit, @spec EVERY/SPEC/ALL, cnf/pdb/trc output; 8 integration tests

**2.3b — Free-energy estimators** — core done
- [x] `bar` — BAR iteration (numerically stable log-sum-exp), bootstrap error
- [x] `ext_ti_merge` — linear interpolation between λ windows, trapezoidal ΔG
- [ ] `reweight`, `m_widom`, `dg_ener` — stubs exist; skip for now (rarely used)

**2.4 — Code quality** ← **NEXT after CellList** (390 warnings hide real bugs)
- [ ] Clippy pass: `gromos-forces` (89), `gromos-integrators` (77), `gromos-io` (31), `gromos-core` (15)
- [ ] Replace bare `unwrap()` in non-test code with `.expect("msg")` or `?`
- [ ] Add missing `#[test]`: SHAKE constraints (0 today), improper dihedral unit test
- [ ] Split large files: `nonbonded.rs` (~1500 LOC), `gromos-io/topology.rs` (~1200 LOC)

**2.5 — Stub cleanup** — parked
- [x] `frameout`, `ener`, `rmsd`, `nhoparam`, `ext_ti_ana`, `bar`, `ext_ti_merge` — real implementations
- [ ] `visco`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` — stubs; parked

### Priority 3 — py-gromos API & education

> **Design target:** `ForceField` → `BuildingBlock` → `Topology` as an algebra; Python expresses verbs; Rust core owns every invariant.

- [ ] Phase 1 — Rust bindings (pyo3-gromos): see `pyo3-gromos/.claude/CONTEXT.md`
- [ ] Phase 2 — Python API: method chaining, energy DataFrame, rich reprs
- [ ] Phase 3 — Notebooks & education: rewrite `py-gromos/notebooks/` + `examples/`

### Priority 4 — Benchmarking
- [ ] Baseline `cargo bench --workspace -- --save-baseline v0.1`
- [ ] End-to-end MD step / pairlist / SHAKE / bonded benches
- [ ] Confirm O(N) scaling after CellList wiring

---

### Cross-cutting — reference tests (do continuously)

Every P1 physics feature and every P2 program lands with a minimal reference test.

**Free win — mine GROMOS's own `check/*.t.cc` regression suite.** The GROMOS devs hard-code
per-term reference energies in `md++/src/check/`: `aladip.t.cc` carries `QuarticBond=18.053811`,
`NonBonded_newRF=-84.092443`, `DistanceRestraint=257.189539`, and perturbed/soft-core terms.
Porting these as unit tests gives per-term validation independent of the md binary — and is a
genuine second source of truth for the perturbed terms.

### Cross-cutting — differential audit (do continuously)

The reference suite is a bug oracle **only for wired paths.** Rules applied to every port:
1. **Reference test BEFORE wiring**, not after.
2. **Grep the C++ for self-flagged defects:** `grep -rniE 'bug|fixme|wrong|hack' interaction/ algorithm/ math/`
3. **Second-source uncertain physics** (RF self-terms, virial): derive from the GROMOS book; diff against C++.
4. **Reproduce genuine GROMOS quirks as named, documented decisions.**

### Deferred breadth (tracked, not scheduled)
- [ ] **PME / lattice-sum electrostatics** — RF stays the focus; PME needs investigation. Note `// wrong!!!` traps in `interaction/nonbonded/interaction/latticesum.{h,cc}` (FUTURE.md Dim 11).
- [ ] **Stochastic / Langevin dynamics** — `random_force` scaffolding exists; SD leap-frog unported.
- [ ] **Coarse-grained → Martini bridge** — FUTURE.md Dim 13; gated on nonbonded-conventions work.
- [ ] **Polarisable / charge-on-spring** — explicitly out of scope.

---

## How to Test

```sh
# Build
cargo build --release --bin md

# Run integration tests
cargo test -p gromos-md --test test_gromosXX_references

# Include ignored systems
cargo test -p gromos-md --test test_gromosXX_references -- --include-ignored

# Run a specific system
cargo test -p gromos-md --test test_gromosXX_references -- pair_lj --exact
```

**Differential audit:**
```sh
grep -rniE 'bug|fixme|wrong|hack' .local/gromosXX/md++/src/interaction/ .local/gromosXX/md++/src/algorithm/
```
