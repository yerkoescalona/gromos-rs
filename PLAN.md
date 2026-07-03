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

**37 of 39 tests pass.** (2 ignored: `aladip_vacuum_fep` — known FEP mismatch; `aladip_vacuum_em` — EM energy frame count off-by-one vs gromosXX)

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

**1.5 — Dim 9: Pairlist from O(N²) to cache-coherent O(N)** ← **NEXT (highest impact)**

> Maps to FUTURE.md §Dimension 9. Five sub-dimensions; 9a is the scaling blocker.
> Dispatch decision: **enum, not trait object** — `pairlist_algorithm.update()` is on the hot path
> every NSNB steps; an enum match inlines to a direct call with zero heap allocation. Two arms now;
> the spatial-index service (9e) is a separate struct, not a third arm.

> **⚠ The math reality that drives the whole test strategy.** Reference tests are NOT byte-identical:
> `test_gromosXX_references.rs` compares every MD step to GROMOS at `ENERGY_REL_TOL = 1e-8`
> (`FORCE_ABS_TOL = 1e-6`). The production pairlist is **not** order-normalized (`normalize_pairs`
> is a `#[cfg(test)]` helper only). So switching Standard→CellList keeps the pair *set* identical
> (already proven by set-equality tests) but changes pair *iteration order* → changes
> floating-point summation order → perturbs each step's energy by ~1e-13 relative → the perturbation
> grows under MD's positive Lyapunov exponent over a trajectory. `water_216_box` runs **100 steps**;
> `ch4_water_fep` runs perturbation on **2998 atoms**. Whether reordered sums stay under `1e-8` to
> step 100 is *empirical, not guaranteed*. Every sub-step below therefore states (a) its exact math
> invariant, (b) the isolated test that proves it *before* any MD wiring, and (c) the "math trap" it
> hides. **Discipline: prove the invariant in a unit test against a brute-force oracle first; only
> then touch the engine.** A sub-step that changes pair *order* must declare whether it targets
> bit-identical (sort to canonical order) or within-tolerance (measure the margin empirically).

**Done**
- [x] `CellListPairlistAlgorithm` — bins chargegroup COGs; O(N) for rectangular boxes, falls back to O(N²) for triclinic/vacuum
- [x] Martina bug NOT reproduced; solute/solvent classification by both atoms' roles
- [x] Set-equality validated vs `StandardPairlistAlgorithm` on all reference systems

---

**9a-0 — Plumbing only, zero float change** ✓ complete
> *Invariant held:* every reference system selects `Standard`; not one summation order changed.
> *Slot mapping:* faithful to gromosXX (`in_parameter.cc:1419-1422`): 0=Standard, 1=ExtendedGrid
> (not yet ported → fallback Standard), 2=CellList (Heinz & Hünenberger). IMD parser also
> recognises `"grid_cell"` keyword. 37/37 reference tests unchanged.
- [x] `PairlistAlgorithm { Standard, CellList }` enum in `gromos-core/src/pairlist.rs` with
  `from_imd()` and `update<BC>()` delegating via match.
- [x] Heuristic: auto-selects `CellList` only when Rectangular + chargegroups + `n_atoms > 5000`.
- [x] Re-exported from `gromos-core/src/lib.rs`.
- [x] 9 unit tests covering full dispatch matrix (algorithm × box_type × n_atoms × chargegroups).
- [x] Field type swapped in `forcefield.rs`, `gamd.rs`, `eds.rs`, `replica.rs`, `md.rs`,
  `simulation.rs`, `algorithm_sequence.rs`. `md.rs` calls `from_imd(imd.algorithm, …)`.
- [x] 37/37 reference tests pass unchanged.

**9a-1 — Turn CellList on; confront reordering head-on** ✓ complete
> *Result:* margin = **0.0 (bit-identical)** across all 100 steps of `water_216_box`.
> CellList iteration order happens to match Standard's on this system. Safe to lower threshold.
- [x] `test_pairlist_margin.rs`: runs `water_216_box` twice (Standard vs CellList via `"grid_cell"`),
  asserts per-step max |ΔE|/|E| < 1e-8. Observed margin = 0.000e0 (bit-identical).
- [x] No auto-heuristic: `ALGORITHM standard` (0) always means Standard. CellList is only
  activated when the input explicitly says `ALGORITHM grid_cell` (2), matching gromosXX
  semantics. All 37 reference files use `standard` → all still run Standard.
- [x] `water_1000_spc_gridcell` reference — 1000 equilibrated SPC (3.1057 nm, vol7 tutorial);
  gromosXX Grid_Cell_Pairlist (7×7×7), 10 steps; gromos-rs CellList passes. Proof chain:
  gromos-rs CellList == gromosXX grid_cell (this) + gromos-rs CellList == Standard (margin=0).
  38/38 reference tests pass.
- [ ] Benchmark `water_216_box` step time: Standard vs CellList — confirm O(N) gap is real.

**9d — Charge groups as first-class queryable primitive** (independent; can land alongside 9a)
> *Invariant:* `cg_table.cog(cg)` is **bit-identical** to `chargegroups[cg].center_of_geometry(pos)`
> (same atoms, same summation order). If true, 9d cannot change any energy — the pairlist sees
> identical COGs. This makes 9d provably energy-neutral.
> *Trap:* a stale cache. The COG cache must be recomputed when positions change; reading `cog()`
> before `update_cogs()` must be unambiguous (return `None`, not a stale/zero value).
- [ ] Add `ChargeGroupTable` to `gromos-core/src/topology.rs`: `groups: Vec<ChargeGroup>`,
  `atom_to_cg: Vec<usize>`, `n_solute: usize`, `cog_cache: Vec<Option<Vec3>>`. Methods:
  `cg_of_atom`, `atoms_of_cg`, `n_solute_cgs`, `is_solute_cg`, `update_cogs(positions)`, `cog(cg)`.
- [ ] Add `pub cg_table: Option<ChargeGroupTable>` to `Topology`; `Topology::build_cg_table()`;
  call it after topology finalization in `gromos-io/src/topology.rs`.
- [ ] **Unit tests (exact, fast, no MD):** (1) round-trip — `cg_of_atom(a)` ∋ `a` in `atoms_of_cg`;
  (2) **bit-exact** — `assert_eq!(table.cog(cg), chargegroups[cg].center_of_geometry(pos))` with `==`,
  not approx, for every CG; (3) staleness — `cog()` is `None` before `update_cogs`; after moving an
  atom + `update_cogs`, `cog` reflects the new position.
- [ ] Switch `CellList`/`Standard` to read `cg_table.cog(i)` (cached) instead of recomputing
  `center_of_geometry` each rebuild. **Validate:** all existing pairlist set-equality tests still pass.

**9e — Cell list as spatial-index service** (after 9a-1 is stable)
> *Invariant:* `pairs_within()` returns exactly the brute-force-N² set for the same cutoff; the
> refactored `CellList` produces the identical pair-set to the pre-refactor one (regression guard).
> *Trap (the big one):* `neighbors_within(query, radius)` with **radius > the grid build cutoff**.
> The grid is cells of size ≥ build-cutoff with a ±1 stencil; a larger query radius reaches beyond
> ±1 cells and **silently misses neighbors**. The API must either expand the stencil to
> `ceil(radius/cell_size)` or reject radius > cutoff. This is a "math won't land" landmine — test it.
- [ ] Add `SpatialIndex` to `gromos-core/src/pairlist.rs` (or `spatial_index.rs`):
  `build(positions, box_dims, cutoff)`, `pairs_within<BC>(...) -> impl Iterator<Item=(u32,u32)>`
  (i<j within cutoff), `neighbors_of<BC>(query_idx, ...)`, `neighbors_within<BC>(query_pos, radius, ...)`.
- [ ] **Unit tests vs brute-force oracle (no MD):** (1) `pairs_within` set == N² set on a random
  small box, periodic; (2) `neighbors_within` set == N² scan, **including a case with
  radius > build cutoff** and a case with radius < cell size; (3) query point not in the particle set.
- [ ] Refactor `CellListPairlistAlgorithm::update_cell_list()` to drive pair classification
  (solute/solvent, short/long, exclusion check) over `pairs_within()`. **Validate:** pre- vs
  post-refactor `CellList` pair-sets identical on all reference systems (regression).
- [ ] Re-export `SpatialIndex` from `lib.rs`.

**9c — Parallel cell build + displacement-triggered rebuild** (after 9e)
> *Parallel invariant:* parallel pair-set ≡ serial pair-set, every run (nondeterministic completion
> order must not change the *set*). *Trap:* if 9a-1 chose canonical sorting for bit-identity, the
> parallel merge must re-sort, or determinism is lost.
> *Trigger invariant (the subtle math):* rebuild when `2·max_disp > skin`. The factor **2** is because
> two atoms can each move `max_disp` toward each other, closing the gap by `2·max_disp`; a pair just
> outside the list at build time can enter the cutoff. Factor-1 silently misses pairs under motion.
- [ ] Parallelize the cell-pair enumeration loop with rayon `into_par_iter()` (already imported,
  `pairlist.rs:11`); each thread emits a partial pair vec, then concatenate (+ re-sort if 9a-1 needs canonical order).
- [ ] Add `use_displacement_trigger: bool` + `ref_positions: Vec<Vec3>` to `PairlistContainer`;
  in `Forcefield::apply()` compute `max_disp = max|pos − ref_pos|`, rebuild when `2·max_disp > skin`,
  refresh `ref_positions` on rebuild. NSNB step-counter stays the default (keeps reference tests fixed).
- [ ] **Tests:** (1) parallel vs serial set-equality, looped many times (nondeterminism guard);
  (2) **hand-built skin test** — one atom displaced just under `skin/2` (no rebuild needed: the pair
  was already listed thanks to skin → energy still exact) vs just over `skin/2` (rebuild fires);
  compare energy to an always-rebuild reference. This is the load-bearing correctness test for 9c.

**9b — Spatial reordering for cache coherence** (last; most invasive)
> *Invariant:* permute∘unpermute = identity ⇒ energies and forces **bit-identical** to the unsorted path.
> *Trap 1:* swapping `sorted_to_orig` and `orig_to_sorted` gives plausible-looking garbage forces.
> *Trap 2 (the one the first draft glossed):* the inner loop's **pair indices must also be remapped
> into sorted space** — positions in sorted order but pairs in original indices = total garbage. Pairs,
> positions, charges, IAC must all live in the same index space inside the loop; forces un-map on exit.
- [ ] Add `SpatialOrder { sorted_to_orig: Vec<u32>, orig_to_sorted: Vec<u32> }` to `pairlist.rs`,
  built from `SpatialIndex` cell assignments after each rebuild.
- [ ] **Unit tests first (no MD):** (1) `unpermute(permute(x)) == x` on random `Vec<Vec3>`;
  (2) `orig_to_sorted[sorted_to_orig[i]] == i` for all i (the two maps are true inverses).
- [ ] Confine permutation to `Forcefield`: working copies `sorted_pos`, `sorted_charges`, `sorted_iac`,
  remapped `sorted_pairs`; do NOT permute `Topology`/canonical `Configuration` (avoids exclusion-index
  rebuild + I/O complexity). Inner loops consume sorted slices; forces un-map via `sorted_to_orig`
  before accumulation into `conf.current_mut().force`. Trajectory writers read original-order `pos` — no I/O change.
- [ ] Gate on `n_atoms > 5000` (cache wins only dominate at scale).
- [ ] **Validate:** energies + forces bit-identical to unsorted on all reference tests (mathematical
  identity if maps are correct); benchmark 50k-atom synthetic box — expect >20% inner-loop speedup.

**Sequencing**
```
9a-0 (plumbing, zero float) → 9a-1 (CellList on, measure margin) → 9e (SpatialIndex) → 9c (parallel + trigger)
                          \                                                                    \
                           └→ 9d (CG table, exact COG cache, independent) ────────────────────→ 9b (spatial reorder, last)
```
> Each arrow is a commit boundary with green tests. 9a-0 and 9d carry zero numerical risk (prove-then-
> wire); 9a-1, 9c, 9b each isolate exactly one source of float divergence and test it against an oracle.

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

**2.4 — Code quality** ← **NEXT after Dim 9** (390 warnings hide real bugs)
- [ ] Clippy pass: `gromos-forces` (89), `gromos-integrators` (77), `gromos-io` (31), `gromos-core` (15)
- [ ] Replace bare `unwrap()` in non-test code with `.expect("msg")` or `?`
- [ ] Add missing `#[test]`: SHAKE constraints (0 today), improper dihedral unit test
- [ ] Split large files: `nonbonded.rs` (~1500 LOC), `gromos-io/topology.rs` (~1200 LOC)

**2.5 — Stub cleanup** — parked
- [x] `frameout`, `ener`, `rmsd`, `nhoparam`, `ext_ti_ana`, `bar`, `ext_ti_merge` — real implementations
- [ ] `visco`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` — stubs; parked

### Priority 3 — py-gromos API & education

> **Approach: small steps, no big design commitment yet.** The full vision (the `System`
> algebra, lazy vs eager `+`, solvation semantics, native topology building) is a *future*
> dimension — see FUTURE.md "Compositional topology construction in py-gromos" and the code
> sketch in `py-gromos/notebooks/00_api_design_mockup.ipynb`. Those open decisions (D1–D8)
> are **deliberately deferred**: we figure them out by walking the path, not by drawing the
> whole map first. P3 below is only the next concrete, decision-free steps.
>
> **The one settled idea:** mirror the two GROMOS files as two Rust objects —
> `Topology` (`.top`) and `Configuration` (`.cnf`) — and a `System = Topology + Configuration`
> that pairs them. Users normally hold a `System`. Start with loading from files; nothing else.

**3.1 — `System.from_files()` — load topology + coordinates as one object** ✓ done
> The smallest useful step. Commits to none of the deferred decisions. Mostly wires
> existing `PyTopology` / `PyConfiguration` together.
> *Invariant:* `System(topo, conf)` raises `ValueError` if `topo.n_atoms != conf.n_atoms`
> — earliest catch of the `.top`/`.cnf` mismatch that silently corrupted vsomm_modeler runs.
- [x] `PyTopology.charge` `#[getter]` — `self.inner.charge.iter().sum::<f64>().round() as i32`.
- [x] `Topology.from_file(path)` / `Configuration.from_file(path)` `#[staticmethod]`
  (keep the existing `__new__(path)` working — just an alias).
- [x] New `PySystem` in `pyo3-gromos/src/system.rs`:
  - Holds `Topology + Configuration`; validates atom-count match at construction.
  - Properties: `n_atoms`, `charge`, `positions`, `velocities`, `box`,
    `topology`, `configuration`.
  - `System.from_files(top, cnf)` `#[staticmethod]` — convenience loader.
  - `write(path)` — write the configuration side to `.cnf`.
- [x] Expose `System` in `__init__.py`; add to `.pyi` stub.
- [x] **Unit test:** load a reference system (`water_216_box`), check `n_atoms`/`charge`;
  mismatched topo/conf raises.

**3.2 — Run an existing `System` with native parameters** ✓ done
> Make a loaded `System` runnable without writing `.imd` files. Reuses the existing
> `PySimulation`; the new piece is parameter factories so Python need not author IMD text.
- [x] `impl Default for ImdParameters` in `gromos-io/src/imd.rs` (sensible GROMOS defaults).
- [x] Factory functions: `nve(dt, steps)`, `nvt(dt, steps, temperature)`,
  `npt(dt, steps, temperature, pressure)`, `steepest_descent(steps)`.
- [x] `PyInputParameters`: `from_file(path)` `#[staticmethod]` (alias of current `__new__`);
  `#[staticmethod]` factory wrappers; `write(path)` to serialise back to `.imd`.
- [x] `Simulation(system, params)` constructor accepting a `System` (in addition to the
  current topo/conf/params form).
- [x] **Unit test:** `nvt(...)` factory produces a runnable sim equivalent to loading the
  matching `.imd`; energies match the reference for `water_216_box`.

**3.3 — Energy out without files: reporters** ✓ done
> Stream energies from a run into Python instead of parsing a `.tre` file.
> *Trap:* reporter callbacks cross the GIL inside the Rust step loop — a panic in a callback
> must surface as a Python exception, not a process abort. *Invariant:* with no reporters
> attached, trajectories stay bit-identical to today (all reference tests pass).
- [x] Wire missing energy components into `PyEnergy`: `angle`, `dihedral`, `improper`
  (fields exist in the Rust `Energy` struct; zeroed today at `simulation.rs:572`).
  Also fixed the same zeroing bug in the `md` binary's `.tre` writer via a shared
  `EnergyFrame::from_energy()` conversion (`gromos-io/src/energy.rs`), now used by
  both `md.rs` and `PySimulation::run()`. Left two pre-existing instances of the same
  bug as documented `TODO(bug)` markers (REMD `replica.rs`, GaMD dihedral boost in
  `md.rs`) — out of scope for P3, not silently fixed.
- [x] `sim.run(steps, ene_freq=100)` — batch loop in Rust returning an
  `(n_frames × k)` numpy energy array (no Python-side per-step loop).
- [x] `EnergyTimeseries` (Python, `py-gromos/python/gromos/timeseries.py`): wraps that
  array; `to_dataframe()` (polars→pandas→dict), `plot(*components)`,
  `block_average(component, block_size)`.
- [x] **Test:** frame count correct; `run()` energies match a `step(1)` loop
  (`test_run_matches_step_loop`); `EnergyTimeseries` behavior (`test_energy_timeseries`).

**3.4 — Notebooks: teach the `from_files` path** (after 3.3) ✓ done
> Replaced the three existing notebooks (they reference phantom APIs like `gromos.State`)
> — deleted `01_understanding_pyo3_bindings.ipynb`, `02_molecular_systems_and_energy.ipynb`,
> `03_performance_deep_dive.ipynb`; none were referenced elsewhere (docs nav, CI-by-name,
> READMEs — confirmed before deleting).
- [x] `01_load_and_inspect.ipynb` — `System.from_files()`, topology/config inspection,
  mass histogram, O-O RDF of bulk water (plotly), round-trip write/reload, the
  atom-count `ValueError` guard, and manual `Topology.solvate()` + `System(topo, conf)`
  composition (needed for `aladip_solvated`, which `from_files()` can't load directly).
- [x] `02_short_md.ipynb` — native `InputParameters` factory vs `from_file()`, batch
  `sim.run()` → `EnergyTimeseries`, per-component energy breakdown plot, block-average
  convergence, `to_dataframe()`, NVE-vs-NVT temperature comparison (manual step loop,
  contrasting batch vs per-step APIs).
- [x] Deprecated `md_runners.py` (docstring notice, points at `Simulation`/`notebooks/02`).
- [x] `pyproject.toml` `notebooks` dependency group: `jupyter`, `ipykernel`, `matplotlib`,
  `plotly`, `polars`. CI (`ci.yml` test-notebooks job) and `Makefile` (`notebooks`,
  new `notebooks-test` target) now rely on `uv sync --all-groups` instead of ephemeral
  `--with` packages.
- [x] `EnergyTimeseries.to_dataframe()`/`.plot()` take a `backend` argument
  (`"polars"`/`"pandas"`/`"dict"`, `"plotly"`/`"matplotlib"`); defaults live on
  `gromos.timeseries.config` (polars + plotly).

**Bug found and fixed while building these notebooks:** the aladip energy-component
plot in `02_short_md.ipynb` came out with `angle`/`dihedral`/`improper` flat at zero —
3.3's fix only corrected the Python-binding *reporting* layer. The actual bug was one
level deeper: `calculate_bonded_forces_ntf` (`gromos-forces/src/bonded/mod.rs`) summed
all four bonded terms into a single `ForceEnergy.energy` scalar with no per-term
breakdown, and every call site (`algorithms/forcefield.rs` — the main MD path used by
every `Simulation.step()`, not just Python — plus `replica.rs`, `gamd.rs`) dumped that
combined total entirely into `bond_total`. Fixed by adding `bond_energy`/`angle_energy`/
`dihedral_energy`/`improper_energy` fields to `ForceEnergy`, populated at the point in
`calculate_bonded_forces_ntf` where each term is still a separate local value, before
combining. `water_216_box`'s existing 3.3 test didn't catch this because SPC water has
no angle/dihedral/improper terms — verified now with `aladip_solvated`, which has all
four. Total potential energy was always correct (it's the sum), so no reference test
regressed; 37/37 `gromos-md` reference tests and all `py-gromos` tests still pass.
Also fixed the matching GaMD dihedral-boost TODO in `md.rs` (now reads the real
`dihedral_total` instead of a hardcoded `0.0`) since it was the same root cause.

**Gap found, not fixed (out of scope here):** `InputParameters.nve()`/`.nvt()`/`.npt()`
factories don't expose SHAKE constraint settings (`ntc` is read-only) — running a
constrained system (anything with H-bonds, e.g. `aladip_solvated`) with factory-built
parameters silently produces an unstable, unconstrained simulation (confirmed: NaN
temperature within ~50-100 steps). `02_short_md.ipynb` documents this and uses
`from_file()` for any system that actually needs constraints. Also: the Python
bindings' `build_simulation()` never dispatches to a steepest-descent algorithm
sequence regardless of `InputParameters.steepest_descent()`/`ntem` — EM is only
reachable via the `md` binary today, not through `Simulation`. No `SteepestDescent`
algorithm-sequence building block is exposed to Python at all.

**3.5 — Minimal fix for the two convenience-path bugs** ✓ done (M1, M2, and SD)
> A composition-pattern redesign (constraints-on-`System`, unify the two builders, consolidate
> defaults, symmetric construction) was drafted and **audited down to its load-bearing core**. The
> audit's finding: only two things are actually broken, and both close in ~20 lines. The rest is
> pre-existing tech debt with no triggering need — deferred below, not scheduled. Ship these two
> first, each with its own verification; do NOT bundle them with the refactor.
>
> *Assumptions exposed by the audit (why the big refactor is not required here):*
> (1) The builder's SHAKE guard fires on `imd.ntc > 1` **alone** (`simulation.rs:107`), so a single
> constraint knob on the factories reaches SHAKE with **zero builder changes** — constraints-on-`System`
> is an aesthetic (OpenMM-orthogonality) choice, not a correctness requirement.
> (2) The temperature-DOF bug is independent of all of it and **unguarded by any test**.
> (3) The 1e-8 Python reference suite exercises `build_simulation` via the 3-arg
> `Simulation(topo, conf, params)` path (`_create_simulation`) — it does *not* A/B-compare the two
> builders, so a "unify builders" rewrite is guarded only insofar as `resolve_algorithm_sequence`
> still matches gromosXX on the 21 covered systems.

- [x] **M1 — constraint knob on the factories.** `ntc_from_constraints()` (`gromos-io/imd.rs`) maps
  `constraints="none"|"hbonds"|"allbonds"` → `ntc=1|2|3`; `ImdParameters::nve/nvt/npt` take it and
  return `Result` (invalid string → error, surfaced as `ValueError` in `pyo3-gromos/parameters.rs`).
  `InputParameters.nve/nvt/npt(..., constraints="none")` (Python default, backward compatible) plus a
  readable `.constraints` getter (inverse of the arg). No builder change — the existing
  `shake_enabled = imd.ntc > 1 || …` guard picks it up automatically.
  **Verified empirically**, not just by construction: factory-built `nvt(..., constraints="hbonds")` on
  `aladip_solvated` (solute H-bonds) stays finite and tracks the 300 K bath target over 300 steps at
  `dt=0.002`; `constraints="none"` on the same system still diverges to >2000 K within 50 steps
  (kept as the documented contrast, not a regression). Test: `test_constrained_system_stable_with_factory_params`
  (`py-gromos/tests/test_gromosXX_references.py`).
- [x] **M2 — fix the `temperature` getter DOF.** Extracted the duplicated DOF calc (was inline in both
  `build_simulation` and `resolve_algorithm_sequence`'s thermostat arm) into one shared
  `compute_total_dof(topo, imd)` (`pyo3-gromos/algorithm_sequence.rs`) — single source of truth,
  called by both builders. `PySimulation` now stores `total_dof` at build time and the `temperature`
  getter uses it instead of bare `n_atoms*3`, so it structurally cannot disagree with what the
  thermostat is coupling to. 37/37 Rust reference tests unchanged (this only touches the getter, not
  any energy path) — confirms the refactor didn't silently change physics.
- [x] **Steepest-descent via `Simulation`.** `build_simulation()` (`pyo3-gromos/simulation.rs`) now
  branches on `imd.ntem > 0`: pushes `Forcefield` → `SteepestDescentAlgorithm` (replacing leap-frog
  velocity+position) → SHAKE (if `constraints=`/`.in` calls for it — GROMOS constrains bonds during EM
  too) → `EnergyCalculation`, skipping COM removal/thermostat/barostat/`TemperatureCalculation`
  (GROMOS convention: EM has no velocities, `E_total = E_pot`) — mirrors the `md` binary's own EM
  sequence construction exactly. Extracted `shake_algorithm_from_imd()` as a shared helper so the EM
  and standard-MD branches build SHAKE identically (can't silently drift apart).
  New `SteepestDescent` algorithm-sequence building block exposed to Python
  (`algorithm_sequence.rs`, mirrors `LeapFrogVelocity`/`ShakeConstraints`), plus a new
  `AlgorithmSequence.minimize(topo, params)` preset and `from_parameters()` dispatching to it when
  `ntem > 0` — the composable path and the direct `Simulation(system, params)` path now agree
  (verified: both produce identical `algorithm_names`).
  **Verified empirically:** `aladip_vacuum` under `InputParameters.steepest_descent(steps=50)` drops
  from -13.5 to -47.98 kJ/mol then correctly plateaus once `dE < DELE` (GROMOS convergence, not a
  bug — the algorithm becomes a documented no-op after convergence). 37/37 Rust reference tests
  unchanged (EM was never reference-tested through Python before). Tests:
  `test_steepest_descent_via_simulation`, `test_steepest_descent_via_algorithm_sequence`.

- [x] **`sim.volume` / `sim.pressure` getters** (added while extending `02_short_md.ipynb`'s NVE/NVT
  comparison to include NPT — there was no per-step way to read either, only via `sim.run()`'s array
  columns). Mirrors the `temperature` getter's shape exactly. **Documented trap, not silently
  papered over:** `pressure` is only physically meaningful under NPT — `P = (2*KE - virial)/(3*V)`,
  and the virial term is only populated by `PressureCalculation`, which only NPT's algorithm sequence
  includes. Under NVE/NVT it still returns a number (the kinetic-only term, ~2000+ bar for
  `water_216_box`), not zero — verified this is genuinely misleading if untagged, so the getter's
  docstring, the API reference, and the notebook (faded/labeled trace) all say so explicitly. Test:
  `test_volume_and_pressure_getters` — cross-checks the getters against `run()`'s array columns, and
  confirms NVE/NVT volume is exactly fixed while NPT's responds to the barostat.

**Deferred — composition-pattern refactor** (tech debt; reappraise, do not schedule yet)
> Full audit prepended to `~/.claude/plans/golden-baking-liskov.md`. These are real but elective, and
> none is the bug the user hit. Revisit only when the FUTURE.md system-building algebra gives
> constraints-on-`System` a *second* caller — i.e. when there's a triggering need beyond elegance.
> - Constraints as a first-class `System` attribute (`constraints=`/`.constrain()`), orthogonal to
>   the ensemble factories (OpenMM-style: constraints on the model, ensemble on the integrator).
> - Unify the two imd→sequence builders: make `build_simulation` (`simulation.rs:66-295`) delegate to
>   `AlgorithmSequence.from_parameters` + `resolve_algorithm_sequence` (deletes ~150 dup lines). Risk:
>   the suite guards each builder against gromosXX but never against *each other*.
> - Consolidate scattered defaults (`300.0/0.1/4.575e-4/0.5/1.0`) into one `defaults` consts module in
>   `gromos-io/imd.rs`; no layer above gromos-io declares a default.
> - Symmetric construction: split `read_imd_file` into pure `parse_imd<R: BufRead>` + path wrapper; add
>   `write_imd_file`; curated `InputParameters(**kwargs)` ctor + get/set properties + `.save(path)`.
>   (Open A/B/C constructor-surface question lives here — premature to decide until this is scheduled.)

**Deferred to FUTURE.md (do NOT start until 3.1–3.4 are solid and usable)**
> The system-building algebra and everything that requires the open D1–D8 decisions:
> `ForceField`/`molecule(seq)` (subprocess `make_top` + coordinate pipeline), the `+`/`*`
> assembly algebra, `.solvate()`/`.pack()`, `.neutralize()`, native topology building, and the
> `system_builder.py` completions. These live in FUTURE.md and the design mockup; revisit with
> running code in hand. The traditional 8 binaries stay available throughout for back-compat.

### Priority 4 — Benchmarking
- [ ] Baseline `cargo bench --workspace -- --save-baseline v0.1`
- [ ] End-to-end MD step / pairlist / SHAKE / bonded benches
- [ ] Confirm O(N) scaling after Dim 9a wiring (see 1.5 benchmark tasks)

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
