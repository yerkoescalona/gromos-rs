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
`PairlistAlgorithm {Standard, CellList}` enum (`gromos-core/src/pairlist.rs`), `from_imd()`/
`update<BC>()` dispatch; slot mapping faithful to gromosXX (`in_parameter.cc:1419-1422`: 0=Standard,
1=ExtendedGrid unported→fallback Standard, 2=CellList). Auto-selects CellList only for
Rectangular+chargegroups+`n_atoms>5000`. Wired into all callers (`forcefield.rs`, `gamd.rs`,
`eds.rs`, `replica.rs`, `md.rs`, `simulation.rs`, `algorithm_sequence.rs`). 9 unit tests; 37/37
reference tests unchanged (every system still selects Standard, zero float change).

**9a-1 — Turn CellList on; confront reordering head-on** ✓ complete
Result: margin = **0.0 (bit-identical)** across all 100 steps of `water_216_box` — CellList's
iteration order happens to match Standard's here. Activated only by explicit `ALGORITHM grid_cell`
(2), not an auto-heuristic, matching gromosXX semantics; all 37 reference files still use `standard`.
New `water_1000_spc_gridcell` reference (1000 equilibrated SPC, real gromosXX Grid_Cell_Pairlist
7×7×7 grid, 10 steps) validates CellList against a real gromosXX oracle, not just against Standard.
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

**2.6 — Dim 12.5: provider-shape groundwork (QM/MM/ML seam)** ✓ classical + ML seam landed (2026-08-13)
- [x] `SpatialIndex` trait + `ConfigurationSpatialIndex` (`gromos-core/src/spatial_index.rs`) —
  query-based neighbor service independent of the MD pairlist's charge-group/twin-range shape,
  carrying periodic-image shift vectors per pair.
- [x] `PotentialProvider` trait (`gromos-forces/src/provider.rs`) — `contribute()` over read-only
  `Topology`/`Configuration` + an `AtomSelection` region, returns a scattered `Contribution` (energy
  + sparse per-atom forces + virial), not raw `&mut Configuration` access (rejected in review: whole-
  system mutable access defeats the "arbitrary but *scoped* atom subset" invariant, FUTURE.md P5).
  Concrete impls named `*Interaction` (matches gromosXX's own `Nonbonded_Interaction`/
  `QMMM_Interaction` convention).
- [x] `LjCrfInteraction` (`nonbonded/interaction.rs`) — first impl, wraps
  `lj_crf_innerloop_novirial` with zero math changes; solute-solute atom-level pairs only.
  `test_provider_reference.rs`: matches real gromosXX `pair_lj` energy to `1e-8` relative.
- [x] `SchNetInteraction` (`nonbonded/schnet.rs`, feature `ml`) — first ML provider: loads a
  TorchScript model via `tch`/libtorch, runs in-process (no Python subprocess/embedded interpreter),
  forces via autograd, neighbor graph from the same `SpatialIndex` classical terms use. Runtime
  choice (`tch`/TorchScript over ONNX) checked against 2026 prior art: GROMOS's own BuRNN embeds
  SchNetPack via `pybind11` (the GIL tax this design avoids); OpenMM-ML and GROMACS's NNPot both
  route through libtorch/TorchScript, not ONNX.
  Requires **libtorch 2.11.0 exactly** (`tch=0.24.0`'s hard-pinned version — PyPI CPU wheel works,
  a conda-forge MKL build failed to even `import torch`). Model I/O is `Dict[str, Tensor]` keyed by
  `schnetpack.properties` (`_positions`, `_atomic_numbers`, `_idx_i`/`_idx_j`, `_offsets`, `_n_atoms`,
  `_idx_m`) — `schnetpack.atomistic.PairwiseDistances` expects a host-supplied neighbor list, not a
  built-in one, confirming the `SpatialIndex` design independently. `scripts/export_toy_schnet.py`
  builds a real (untrained) SchNetPack 2.2.0 `NeuralNetworkPotential`; **SchNetPack v1 does not
  `torch.jit.script`** (`schnetpack/nn/neighbors.py`'s conditionally-typed return — confirmed against
  the real trained BuRNN tutorial model, `.local/gromos_tutorial_livecoms` t_06), v2's `Forces`
  module does once `calc_forces=True` populates `required_derivatives`, and computes forces inside
  the scripted graph (no separate `.backward()` call needed). Validated against central finite
  differences (`5e-3` absolute — measured float32 noise floor is ~1.4e-3 across 9 components on a
  3-atom fixture) since no gromosXX oracle exists for a neural-network potential (FUTURE.md P8).
  `schnet_burnn_reference.rs`: real 6-atom QM zone from `t_06`'s `eq_meoh_5.cnf`, real `QMZONE`/
  `RCUTQM=1.4` from `meoh.qmm`/`md.imd`. **Caught a real bug:** `contribute()` assumed every
  `neighbor_pairs()` result had both endpoints in `region` — false once the region sits inside a
  solvated system (`neighbor_pairs` deliberately also returns boundary-crossing pairs, for the future
  embedding case). Fixed by skipping outside-region pairs; documented as a current limitation (no
  embedding support yet). 67 `--features ml` tests green; classical 37/37 unaffected throughout.
- [ ] **Deferred, not started:** wiring any MD binary's force loop to iterate
  `Vec<Box<dyn PotentialProvider>>`; a real *trained* model (this pass's is untrained, seam-only);
  async/cancellable `contribute` for external-process providers; a shared typed-units boundary.
  `crates/gromos-forces/src/qmmm.rs` (older file-in/file-out QM scaffolding, forks `xtb`) predates
  this design and isn't its basis — left untouched, fate undecided.

**2.7 — Electrostatic embedding onto MM atoms (QM/MM + BuRNN)** ✓ Steps 1-5 landed (2026-08-13)

> 2.6 left `SchNetInteraction` skipping every cross-boundary neighbor pair, i.e. no environment
> coupling at all. Grounded in Poliak et al. 2025 (*J. Comput. Chem.* 46:e70053, mainline GROMOS
> QM/MM) and Gómez-Flores et al. 2022 (*J. Chem. Theory Comput.* 18:1213, BuRNN), plus direct
> inspection of the tutorial's training code and shipped model — two different mechanisms:

| | Classic EE (Poliak mainline) | BuRNN (what `t_06` runs) |
|---|---|---|
| MM charges → model | **yes**, polarizes the QM density | **no** |
| Model → charges | no | **no** (confirmed: shipped model is a single `Atomwise` output, no charge head) |
| Zones | QM + MM (+ link atoms) | inner + buffer + outer |
| Force path onto MM | QM forces / `f=qE` / pairwise Coulomb from QM charges | buffer gets NN forces; outer stays classical |
| Main risk | link-atom boundary (out of scope, matches gromosXX: no charge redistribution) | double-counting inner/buffer |

Confirmed by reading `train_dataset_tutorial/mopac.py` directly: BuRNN's training target is
`E_burnn = E_QM(inner+buffer) − E_QM(buffer) − E_solute_vac` (QM−QM via MOPAC, **not** QM−MM) — the
NN carries inner internal energy + inner↔buffer interaction + buffer polarisation response, not the
buffer's own internal energy. `get_burnn_forces()` agrees per-atom. Double-counting is therefore a
static training-time contract (A5) the orchestrator must reproduce as a fixed decomposition, not
something detectable at runtime.

- [x] **Step 1 — gather test proves the trait needs no signature change.**
  `tests/embedding_gather_reference.rs`: on the real `t_06` system, **1363 MM atoms** found within
  `RCUTQM`=1.4nm of the 6-atom QM zone via `neigh.neighbor_pairs`, charges/positions readable —
  `contribute()`'s existing `topo`/`conf`/`neigh` args and `Vec<(usize, Vec3)>` forces already carry
  everything EE needs.
- [x] **Step 2 — `Embedding` becomes a stated, type-level property.** `provider.rs` gained
  `Embedding {None, Mechanical, Electrostatic}` + a defaulted `PotentialProvider::embedding()`.
  `SchNetInteraction` declares `None` (its `Dict[str,Tensor]` contract has no environment channel)
  and now rejects an unsupported scheme with a clear error instead of silently dropping pairs.
- [x] **Step 3 — `ElectrostaticEmbedding` (`nonbonded/embedding.rs`).** Poliak path (c): QM charges +
  pairwise Coulomb computed in the MD code. First provider to put forces outside its own region
  (FUTURE.md P5). Validated three ways: closed-form two-point-charge oracle (exact), FD on an MM
  atom, and FD on a real `t_06` MM atom (E_embed = 508.06 kJ/mol, forces on all 1363 MM atoms,
  equal-and-opposite, no intra-region double-count).
- [x] **Step 4 — `zones.rs`: `Zone {Inner, Buffer, Outer}` + `PairOwner` + `ZonePartition::owner()`.**
  The full six-pair-class ownership table derived directly from `get_burnn_energy()`'s decomposition
  — an orchestration type, not logic hidden in a provider. Degrades to plain QM/MM with no buffer,
  fully classical with no zones. 6 unit tests.
- [x] **Step 5 — replaced: the planned oracle doesn't exist.** Originally scoped as "compare against
  the MM baseline the training pipeline used" — checked first and found there is no such baseline
  (training subtracts two MOPAC energies, no classical term appears anywhere in the pipeline,
  `get_reference_energies()` is never called). Building that test would have manufactured false
  confidence, so it was dropped rather than faked. **Replaced** with
  `zone_partition_reference.rs`: asserts Step 4's table *partitions* the real system — no pair
  claimed twice or zero times. On real `t_06` (`BUFFERZONE`=0.5nm, `RCUTQM`=1.4nm): 3513 atoms →
  inner 6 / buffer 99 / outer 3408; all 6659 QM-cutoff pairs split provider 609 / embedding 6050 /
  classical 0 cleanly.

**Deferred (named so they are not silently assumed)**
- [ ] Charge-output channel on `Contribution` — blocked on a model that actually predicts charges (A3).
  Note it would break the pure-additive provider model: one provider's output becomes another's input
  within a step, i.e. ordered evaluation, not summation. Do not add it speculatively.
- [ ] Link atoms / capping / boundary charge redistribution (A4).
- [ ] Mutual polarisation (MM polarised by QM *and* back) — Poliak: EE "includes electronic
  polarization of the QM zone by the MM atoms, but not vice versa". Out of scope for parity.

**2.8 — Wiring: from "the pieces exist" to "a QM/MM run happens"** ✓ P2.8-1..4a, 2b, 6 done;
P2.8-4b blocked on the user, P2.8-5 deliberately deferred

> After 2.7, every *piece* was built and validated in isolation but nothing was connected —
> `ZonePartition` had no non-test caller, `LjCrfInteraction` didn't filter by zone, no binary
> iterated providers. P2.8-1..4a below closed those one at a time. **Originally still open here:**
> no binary drove `ProviderOrchestrator`, and the QM/ML provider's own inner-zone energy was
> missing from every real-system test because this environment had no QM engine. Both closed: `xtb`
> is now `apt`-installable here (P2.8-2b), and `ProviderOrchestratorAlgorithm` wires the
> orchestrator into the real `AlgorithmSequence` step loop (P2.8-6). `SchNetInteraction`'s
> `--features ml` build was also unblocked later the same session: `libtorch` 2.11.0 (CPU wheel,
> exact pinned version) + `schnetpack` 2.2.0 installed into a local venv
> (`python3 -m venv --without-pip` + `get-pip.py`, since `python3.13-venv` needed root and wasn't
> available — see `schnet.rs` module docs for the full recipe); real `t_06` BuRNN reference test
> and all `nonbonded::schnet::tests` now run for real here, not skipped.

- [x] **P2.8-1 — `LjCrfInteraction` honours `ZonePartition`.** `pub zone_partition:
  Option<ZonePartition>` + `with_zone_partition()` builder (`interaction.rs`); `contribute()` filters
  through `classical_should_evaluate()`. `None` (default) preserves prior unpartitioned behaviour
  exactly. **Exit:** on real `t_06` (3513 atoms), partitioned LJ+CRF (1970995 kept pairs) + the sum
  over the 6659 excluded pairs reconstructs the unpartitioned full-system energy (1977654 pairs) to
  `1e-8` relative (`zone_partitioned_lj_crf_reference.rs`). 37/37 + provider reference unaffected.
  **Bug caught building the exit test:** ~0.5% energy mismatch traced to `Configuration::new()`
  defaulting `box_config` to vacuum — the test set `pos` but not `box_config`, so the provider
  computed unwrapped distances while the test's own spatial index correctly wrapped. Fixed by
  setting `box_config` from the real `.cnf`; documented as a trap for future bare-`Configuration`
  tests. **Open watch:** `Contribution.energy` is a single scalar, so per-term GROMOS accounting
  (lj/crf split) is lost at the provider boundary — not decided, not blocking yet.
- [x] **P2.8-2 — `ProviderOrchestrator` sums `Vec<Box<dyn PotentialProvider>>`.**
  `orchestrator.rs`: `register(region, provider)` + `evaluate()`, enforcing that an `Embedding::None`
  provider never places a force outside its own region (hard error, the P2.6 review finding this
  exists to check) — `Mechanical`/`Electrostatic` providers exempted by definition. 5 unit tests
  including a deliberately-lying mock provider. **Exit, revised before building:** the literal
  "reproduces 37/37" premise doesn't survive `LjCrfInteraction`'s solute-solute-only scope (33/37
  systems are unreachable by any orchestrator built from today's providers — a roster gap, not an
  orchestration bug). Checked instead: the 4 pure two-atom no-solvent systems
  (`pair_lj`/`pair_lj_mixed`/`nacl_pair`/`nacl_pair_box`) match both a direct non-orchestrated call
  and gromosXX's real energy to `1e-8` (`test_orchestrator_reference.rs`).
  **Bug caught:** `nacl_pair`/`nacl_pair_box` initially off by ~21% — gromosXX's reported energy
  includes the unconditional per-atom RF self-energy term that `LjCrfInteraction` deliberately
  excludes (separate GROMOS term); applied the same correction `interaction.rs`'s own oracle test
  already established.
  **Follow-up, same session as P2.8-4:** `qmmm_orchestrator_reference.rs` registers
  `LjCrfInteraction` (zone-partitioned) + `ElectrostaticEmbedding` together on real `t_06` and
  matches the two providers' direct calls added by hand — transparency extended to two providers on
  a real partitioned system. Still missing: the QM/ML provider's own inner-zone energy (no
  `libtorch` here).
- [x] **P2.8-2a — `SchNetInteraction` driving a real MD step loop (single provider, no
  orchestrator).** `schnet_nve_loop.rs` (feature `ml`): real leapfrog NVE dynamics, 500 steps
  (force eval → `LeapFrog::step` → force eval), not a single isolated call. Energy conservation as
  the check (architecture-blind to chemical accuracy, same FUTURE.md P8 tier as 2.6). `dt=1e-3`
  chosen empirically (scanned `1e-4..2e-3`; largest step giving real KE while keeping fluctuation
  at 0.0047% of the mean); tolerances set ~100x that measured margin. Proves `PotentialProvider` +
  `LeapFrog` compose correctly — narrower than P2.8-2's orchestrator, doesn't attempt multi-provider.
- [x] **P2.8-2b — `XtbInteraction`: a real `xtb` QM engine as a `PotentialProvider`.** Closes the
  other half of the gap the P2.8 intro note above flagged ("this environment has no QM engine
  (`xtb`/`mopac`)") — `xtb` 6.7.1-2 is now `apt`-installable here. `nonbonded/xtb.rs`: subprocess
  file-in/file-out, not the C-API linking `gromosXX`'s real `xtb_worker.cc` uses (that's a much
  bigger FFI-binding lift, deferred). **Grounded in real xtb runs, not docs alone** — the older
  `qmmm.rs::XTBEngine` scaffold (predates `PotentialProvider`, never wired up) parsed a `energy`/
  `gradient` file pair whose actual format doesn't match what it assumed; ran real GFN2-xTB on a
  toy water molecule and found: (a) xtb's xyz parser accepts a bare atomic number in the element
  column (verified identical energy either way) — so, like `SchNetInteraction`, atomic numbers are
  a caller-supplied table, no symbol lookup needed; (b) `<basename>.engrad` (ORCA-style: atom
  count, energy in Eh, flattened gradient in Eh/bohr, `#`-prefixed comments) is a cleaner parse
  target than the Turbomole `$grad`/bare `energy` files the old scaffold used. Added
  `HARTREE_TO_KJMOL`/`BOHR_TO_NM` to `gromos_core::units` (didn't exist) rather than reusing the
  old scaffold's unchecked magic constant. **Scope decision:** only `Embedding::None` (the region's
  own isolated-cluster energy) — `Mechanical`/`Electrostatic` rejected loudly, same policy
  `schnet.rs` already uses, so this can't double-count against `ElectrostaticEmbedding`'s own
  inner↔outer coupling (`zones.rs` A5 hazard). xtb natively supports point-charge embedding with
  forces returned on the MM side too (`pcharge`/`pcgrad` files, confirmed by reading
  `xtb_worker.cc` and testing the file format directly), but wiring that up was deferred at the
  time. **Done later this session as P2.8-2d** — the "harder charge-refresh problem" turned out to
  be smaller than expected; see that entry and the corrected P2.8-5.
  **Exit:** `nonbonded/xtb::tests` — a pinned real-oracle energy test (`-5.018180941704 Eh` for a
  specific water geometry on this machine's xtb 6.7.1-2), a central finite-difference
  self-consistency check on forces (validates this provider's own sign/unit-conversion pipeline,
  independent of xtb's physics — FUTURE.md P8 tier, same as `schnet.rs`'s own FD test), and the
  same loud-rejection-of-unsupported-embedding test `schnet.rs` has. `xtb_nve_loop.rs`
  (`gromos-md`): real leapfrog NVE loop driving `XtbInteraction` alone on a water monomer, 120
  steps at `dt=2e-4` (real GFN2 PES is far stiffer than the untrained SchNet toy model, needed a
  smaller step than P2.8-2a's `1e-3`) — 0.038% energy fluctuation, well inside the 0.5% bound.
  Not feature-gated (no compile-time linking like `ml`/`tch` needs), so it builds everywhere;
  tests skip cleanly if `xtb` isn't on `PATH`.
- [x] **P2.8-2c — the first real, combined QM/MM evaluation.** `qmmm_orchestrator_reference.rs`
  (earlier work) already proved `LjCrfInteraction` (zone-partitioned) + `ElectrostaticEmbedding`
  combine correctly on a real system, but its own module doc said outright: "this environment has
  no QM engine... not that the result is a complete QM/MM energy," and it used hard-coded stand-in
  charges for the embedding term. `xtb_qmmm_water_dimer_reference.rs` closes both: registers
  `XtbInteraction` + `ElectrostaticEmbedding` + zone-partitioned `LjCrfInteraction` together, with
  the embedding charges coming from a real xtb Mulliken-charge calculation, not a stand-in.
  **System:** `water_dimer` (Poliak dataset, already checked into
  `gromosXX_qmmm_references/water_dimer_mechst/`), not `t_06` — `t_06`'s 54a7 force field is
  united-atom, and its QM zone's true per-atom element identity lives in a `.qmm` file's
  `QMZONE`/`QMZ` column that nothing here parses yet; `water_dimer`'s two explicit-atom SPC waters
  need no such mapping (`[8,1,1]` per water, trivially).
  **Charges:** read directly from xtb's own `charges` side-effect file after a real
  `XtbInteraction::contribute()` call, rather than adding a charge-output channel to
  `Contribution` — that channel is explicitly deferred in `provider.rs` (P2.7: "would break the
  pure-additive provider model... do not add it speculatively"), so this reads xtb's own file
  instead of building it.
  **Exit:** real numbers, not just an assertion of correctness — for the real `water_dimer`
  starting geometry: `E_QM = -13307.5`, `E_classical = -735.5`, `E_embedding = -23.1`,
  `E_total = -14066.2` kJ/mol; real Mulliken charges `[-0.562, 0.280, 0.282]` (sums to ~0, a
  neutral water, as a sanity check); **3 QM-zone atoms and 3 MM atoms both carry nonzero force** —
  the concrete, checkable version of "QM and MM influence each other," not present in any earlier
  test. Orchestrator transparency (combined == three direct calls summed) still holds with all
  three providers, extending `qmmm_orchestrator_reference.rs`'s two-provider check.
  **Explicit non-goals, stated in the module doc so they aren't mistaken for gaps:** doesn't match
  gromosXX's real published number for this system (their reference run used AM1, xtb is a
  different method — same accepted mismatch as P2.8-2b); charges are static, one xtb call, not
  refreshed per step — this test deliberately keeps using path (c) (`ElectrostaticEmbedding`) to
  demonstrate three-provider orchestration; path (a) (real per-step-consistent embedding, no
  static charges at all) landed later this session as P2.8-2d and is demonstrated alongside it in
  the same test file; single-point, not a multi-step run (the waters are SHAKE-constrained rigid
  bodies in the real `.imd` — a dynamical version needs the constraint solver wired into the
  sequence too, a real separate follow-up, not attempted here).
- [x] **P2.8-2d — real electrostatic embedding in `XtbInteraction` (Poliak path (a)).** Grew out of
  a user question: "isn't the static-charge approximation and the per-step charge-refresh problem
  something other implementations already solved?" Read gromosXX's real `xtb_worker.cc` again with
  that question specifically — it already links xtb's own C API (`xtb_setExternalCharges`,
  `xtb_getPCGradient`), which is path (a): MM point charges polarize xtb's own SCF, and xtb hands
  back MM forces directly, so there's no charge-derivative term to solve at all (see the corrected
  P2.8-5 above). `nonbonded/xtb.rs`'s `XtbInteraction` now implements the same physics via xtb's
  CLI-level `pcharge`/`pcgrad` files (subprocess path, not the C API — consistent with P2.8-2b's
  original choice to avoid a native build dependency).
  **Units determined empirically, not assumed** — xtb's man page documents `pcharge`'s columns but
  not their length unit. Verified directly: placed a real point charge at two known distances
  (10 and 50 Bohr) from a bare Na⁺ ion and compared the measured energy shift to the analytic
  monopole-monopole Coulomb energy (`q₁q₂/r` in atomic units) under both Bohr and Ångström
  hypotheses. Bohr matched to 99.85% at r=50 (residual is real, expected SCF polarization —
  shrank correctly with distance); Ångström was off by the exact expected ~1.89× ratio. Also
  confirmed `pcgrad`'s row order matches `pcharge`'s input order and its sign convention
  (force = −gradient, same as `.engrad`) with a real two-point-charge run.
  **Scope:** `Embedding::Electrostatic` now real; `Embedding::Mechanical` still rejected (nothing
  implements it). This and `ElectrostaticEmbedding` (path (c)) are **alternatives, not additions**
  — documented prominently in both modules to prevent the double-counting `zones.rs` (A5) exists to
  catch. `ZonePartition::owner()` already unconditionally routes inner-outer electrostatics to
  `Embedding` regardless of which path handles it, so a zone-partitioned `LjCrfInteraction` composes
  correctly with either path with no extra wiring.
  **Exit:** `nonbonded/xtb::tests` — a real polarization test (`electrostatic_embedding_changes_
  qm_energy_and_puts_force_on_the_mm_atom`: an external −0.8e charge measurably shifts the QM
  energy and receives a real force back, which path (c) structurally cannot do since it can only
  ever add a fixed Coulomb term to a fixed QM energy); a finite-difference check covering *both*
  the QM atoms' and the MM atom's forces (the tier that validates this path's own sign/unit
  pipeline, independent of xtb's physics); loud rejection of `Mechanical` and of `Electrostatic`
  without a cutoff set. `xtb_qmmm_water_dimer_reference.rs` gained a second test,
  `path_a_real_electrostatic_embedding_on_the_real_water_dimer`, running the real water-dimer
  system through path (a) instead of path (c) — total energy agrees with the path (c) demo to
  within 0.09% (−14053.7 vs −14066.2 kJ/mol), the physically expected gap between real
  self-consistent polarization and a frozen-charge approximation of it.
- [x] **P2.8-2e — `MopacInteraction`: a real MOPAC engine, AM1/PM3/PM6/PM7 all for free.** Grew out
  of the user's original AM1 question — MOPAC (apt-installable, `mopac` 23.1.2) is the real engine
  behind AM1/PM3/PM6/PM7; `method` is a plain field on `MopacInteraction`, so all four come from one
  implementation, not four ports.
  **Architecture check, done before writing any code:** does a second real QM engine mean a new
  `QmEngine` trait (mirroring gromosXX's `QM_Worker` base class, ~20 virtuals, 9 subclasses)?
  `architecture.md` rule 3 says no — "a small, stable trait taxonomy — resist proliferation... about
  five contracts." `MopacInteraction` is just another `PotentialProvider`, like `XtbInteraction`.
  The real, mechanical duplication between the two (work-dir setup, stale-output removal, subprocess
  spawn/exit-check) moved to a new crate-private `qm_subprocess.rs`, shared by both (`XtbInteraction`
  lightly refactored to use it too, behavior unchanged — confirmed by its full existing test suite
  still passing unmodified).
  **Oracle checkpoint — the plan's own explicit gate, and it failed as originally scoped.** The
  BuRNN tutorial ships 1722 real archived MOPAC `.aux` files (`MOPAC_VERSION=MOPAC2016.22.067L`,
  PM7) — checked before building anything on top: re-ran one archived geometry through this
  machine's installed MOPAC 23.1.2. Single-point energy at the archived geometry: **−1427.7
  kcal/mol**; the archive itself reports **−1610.9 kcal/mol** — 165 kcal/mol (11%) apart. Diagnosed,
  not just noted: re-optimizing MOPAC 23 from that same starting geometry converges to **−1592.1
  kcal/mol** (close to the archive's −1610.9, ~1.2% — each version's *own* minimum agrees
  reasonably), and a fresh single-point exactly at that MOPAC-23-optimized geometry reproduces
  −1592.1 to 5 decimal places (MOPAC 23 is internally self-consistent — not a bug here). Conclusion:
  the archived geometry, real MOPAC2016 stationary point though it is, isn't close to a stationary
  point on MOPAC 23's own PM7 surface — the two versions' surfaces differ enough that cross-version
  *fixed-geometry* comparison is unreliable, even though independent optimization lands in the same
  ballpark. **Consequence:** the 1722 files aren't usable as tight single-point oracles. Pivoted to
  the pattern that already worked for `XtbInteraction` (P2.8-2b): pin real values from *this
  machine's own* MOPAC, not the 2016 archive.
  **`.mop`/`.aux` gotchas found against real output, not docs:** MOPAC's `GRAD` keyword reports
  nothing ("Keyword GRADIENTS used, but geometry has no variables") unless every atom's coordinate
  flags are `1` (optimization variables) — `1SCF` then prevents any actual optimization from moving
  the geometry, matching gromosXX's own tutorial `mopac.py` input convention exactly. `.aux` values
  use Fortran `D` exponents (`-0.161088149248928D+04`), rejected by `f64::from_str` — normalized to
  `E` first.
  **Embedding:** `Embedding::None` only — gromosXX's real `mopac_worker.cc` doesn't implement
  `parse_mm_gradients` (unlike `xtb_worker.cc`), so MOPAC never returns MM forces; pairing it with
  QM/MM coupling means path (c) (`ElectrostaticEmbedding`), documented in both modules.
  **Exit:** `nonbonded/mopac::tests` — a pinned water PM7 energy test (real MOPAC 23.1.2 run on this
  machine), loud rejection of unsupported embedding, a finite-difference force check (same tier as
  every other QM-engine test here). `tests/mopac_reference.rs` — a real multi-element (C, O, H)
  methanol molecule, energy *and* Mulliken charges (summed to `1e-14`, confirming neutrality) both
  matching a pinned real-MOPAC-23 run. Geometry is a standard, self-constructed methanol (typical
  bond lengths/angles), not lifted from the BuRNN archive — that repository has no explicit license
  file for its data, unlike the Poliak/Zenodo dataset (CC BY 4.0) already used elsewhere in this
  suite, so nothing from it is checked into this repo's git history.
- [x] **P2.8-3 — virial from `ElectrostaticEmbedding`.** Was `Mat3::ZERO` (NPT with an embedded QM
  zone would silently mispressure). Added per-pair atomic virial matching `ForceStorage`'s
  established convention. **Exit, three tiers:** closed-form exact (`trace(virial)==energy` via
  Euler's theorem for a `1/r` pair), FD on the scaling derivative (`1e-6`), and the real periodic
  `t_06` system (1363 pairs, `1e-8` relative) — the real-system check is the sharpest PBC-
  correctness test since a wrong minimum-image convention would break the identity even if energy
  and forces individually looked plausible. 82 lib tests green; 37/37 unaffected (unwired).
- [x] **P2.8-4 — resolved the `QMLJ` ambiguity.** Cloned real gromosXX (`.local/gromosXX`) and read
  `in_parameter.cc::read_QMMM` + `qmmm_interaction.cc::modify_exclusions` directly. **Answer:**
  `QMLJ` gates a classical LJ *supplement* for inner-inner/inner-buffer pairs only; inner-outer LJ
  is always classical regardless (only its electrostatics moves to embedding) — the opposite of the
  crate's prior undocumented assumption. Extended `ZonePartition` with `qm_lj: bool` (default
  `false`, matches `t_06`) + `lj_owner()` (separate from `owner()`) + `lj_only_should_evaluate()`;
  wired a second CRF-zeroed LJ pass into `LjCrfInteraction`. **Exit:** real `t_06` three-way
  decomposition (1970995 classical + 6050 lj-only-supplement + 609 fully-excluded = 1977654 total)
  reconstructs `full` to `1e-8`. 87 lib tests; 37/37 unaffected; release build clean.
- [x] **P2.8-4a — a real external gromosXX QM/MM oracle (not just source/synthetic).** Found the
  Poliak et al. 2025 dataset on [Zenodo 10.5281/zenodo.14549978](https://doi.org/10.5281/zenodo.14549978)
  (CC BY 4.0) — real production QM/MM runs (`WATER`/`WATER_DIMER`/`WATER_SINGLE`/`AA_SOLV`/
  `TRIPEPTIDES`, 4 semi-empirical methods, 3 embedding schemes). Most needs `mndo` (not
  apt-installable here, same blocker class as `xtb`/`libtorch`). **Reachable slice:**
  `WATER_DIMER` mechanical-embedding-constant-charge (`mechst`) needs no QM engine at all — gromosXX
  leaves QM-MM exclusions untouched for this scheme. `tests/gromosXX_qmmm_references/
  water_dimer_mechst/` (real topology/config/imd/energies, extracted via ranged HTTP reads from the
  980MB archive without downloading it whole; README carries provenance). `single_point_energy`
  reproduces the real gromosXX `nonb` to `1.3e-9` relative.
  **Deferred, not started:** everything needing an actual QM program (`elst`/`mech` runs,
  `AA_SOLV`/`TRIPEPTIDES`) — genuinely unreachable without `mndo`, not attempted speculatively.
- [ ] **P2.8-4b — unlock the full Poliak dataset via a real MNDO build.** Blocked on a step only
  the user can do: MNDO (`mndo.mpi-muelheim.mpg.de`) is "open source for non-commercial use," but
  distributed through a license form + Gitea instance requiring registration — not an anonymous
  `apt`/`git clone`, and not something to complete on the user's behalf. **Blocked, specifically:
  the user needs to be at a machine they don't currently have access to, to do that registration.**
  Waiting on them, not on any further investigation here.
  **Exit, once unblocked:** user hands over MNDO source (Gitea URL+token, a tarball, or a local
  path) → build MNDO (`gfortran`, confirmed not installed here but apt-installable) → build
  gromosXX from the `.local/gromosXX` clone wired to that binary (`NTQMSW=0`) → run the real
  Poliak `am1`/`om2`/`om3`/`pm3` `elst`/`mech` `WATER_DIMER` systems (and potentially
  `AA_SOLV`/`TRIPEPTIDES`) → reproduce their published energies as exact bit-for-bit oracles, the
  same rigor P2.8-4a's `mechst` slice already has, but for the actual electrostatic-embedding case
  `ElectrostaticEmbedding` was built for. `xtb` (apt-installable now, no registration) was offered
  as a faster alternative path — declined in favor of this one, since it wouldn't match Poliak's
  published numbers (different QM method).
- [x] **P2.8-5 — refresh QM charges per step. Superseded by P2.8-2d, not implemented as
  originally scoped.** Originally framed as "derive gromosXX's `MEDC` extra-force correction for
  fluctuating charges, or restrict to static charges." Read gromosXX's real source
  (`interaction/qmmm/qm_zone.cc`) directly rather than trusting the earlier paraphrase of Poliak's
  paper: **no such correction term exists anywhere in gromosXX.** Dynamic per-step charges exist
  only for *mechanical* embedding (`qm_zone.cc:177-196`) — the engine's charges are written
  straight into `topo.charge()` and picked up next step by the ordinary classical loop, `dq/dR`
  simply neglected. For *electrostatic* embedding, gromosXX doesn't refresh `ElectrostaticEmbedding`-
  style static charges at all — it uses path (a) instead (P2.8-2d): the QM program returns MM
  forces directly, so there is no charge-derivative problem to solve, because nothing is
  differentiating a *fixed*-charge approximation in the first place. That's the real fix, and it's
  done for `xtb` (`XtbInteraction`'s `Embedding::Electrostatic`, `nonbonded/xtb.rs`).
  **What's left, correctly scoped down:** engines that only support path (c) (no `parse_mm_gradients`
  — confirmed only `xtb_worker.cc`/`orca_worker.cc`/`turbomole_worker.cc` implement it; MOPAC/MNDO/
  DFTB don't) still only have `ElectrostaticEmbedding`'s static charges available. Refreshing those
  per step would be a real, much simpler feature than originally scoped (no derivative term needed —
  just rewrite `region_charges` from a fresh engine call each step, mirroring gromosXX's own
  mechanical-embedding behavior) — a legitimate but smaller future item, not started.
- [x] **P2.8-6 — wire `ProviderOrchestrator` into a running binary's step loop.** Named as the real
  remaining gap the P2.8 intro note above already flags. Every provider/orchestrator/zone piece up
  to P2.8-2b was validated by test harnesses calling `orchestrator.evaluate()` directly, or (P2.8-2a/
  P2.8-2b) a bespoke `LeapFrog` loop — nothing went through the real step machinery.
  **Premise correction made while starting this:** there is no single Rust `Simulation` struct.
  `md.rs::main()` and `pyo3-gromos/src/simulation.rs::build_simulation()` each independently build a
  `gromos_core::algorithm::AlgorithmSequence` (`Vec<Box<dyn Algorithm>>`) and drive it with
  `run_step()` — *that* is the real shared entry point, not a `Simulation` type. The standard
  sequence both build: `RemoveCOMMotion → Forcefield → LeapFrogVelocity → thermostat →
  LeapFrogPosition → constraints → TemperatureCalculation → PressureCalculation/Barostat →
  EnergyCalculation`.
  **What was built:** `gromos-forces/src/orchestrator_algorithm.rs` — `ProviderOrchestratorAlgorithm`
  implements `Algorithm` directly (no new gromos-integrators↔gromos-forces dependency edge needed;
  gromos-forces already depends on gromos-core, where `Algorithm` lives), wrapping a
  `ProviderOrchestrator` + `Periodicity` (refreshed from the box each step, mirroring `Forcefield`'s
  own NPT-safe refresh block verbatim). Pushed immediately after `Forcefield` in a sequence, it adds
  (`+=`) to `Forcefield`'s already-computed force and virial.
  **Energy-bookkeeping trap found and fixed:** the summed contribution can't go into
  `potential_total` — `update_potential_total()` (called by both `Forcefield` and the final
  `EnergyCalculation`) only sums bond/angle/dihedral/improper/cross_dihedral/lj/crf/ls/sasa, so a
  direct write there would be silently wiped out by that later recompute. `Energy::total()` adds
  `special_total` on top separately (same slot position restraints use), so that's where it goes
  instead. A second trap: `special_total` isn't unconditionally reset to a fresh per-step baseline
  by anything else (`Forcefield` only *sets* it, conditionally, when position restraints are
  configured) — a naive `+=` here would silently compound across steps whenever they aren't. Fixed
  by overwriting `special_total` rather than accumulating; documented as a known limitation
  (`ProviderOrchestratorAlgorithm` currently claims exclusive ownership of that field — running it
  alongside active position restraints in the same sequence isn't supported yet).
  **Exit:** `gromos-md/tests/xtb_orchestrator_sequence.rs` — a real `AlgorithmSequence` (`Forcefield`
  + `ProviderOrchestratorAlgorithm` registered with a real `XtbInteraction` + `LeapFrogVelocity` +
  `LeapFrogPosition` + `TemperatureCalculation` + `EnergyCalculation`) run for 120 steps via
  `run_step()` on a water monomer whose classical `Forcefield` contributes exactly zero (zero
  charge/LJ/bonds) — deliberately narrower than "both nonzero simultaneously" (same incremental
  discipline P2.8-2a used), isolating what needed proving: the orchestrator's energy/force/virial
  actually flow through the real pipeline. NVE energy conservation via `conf.old().energies.total()`
  (GROMOS convention — results land in `old()` after leapfrog's state exchange): 0.0032% drift over
  120 steps, tighter than P2.8-2b's bespoke-loop equivalent (0.038%). A synthetic-provider unit test
  in `orchestrator_algorithm.rs` separately checks the `+=`-vs-overwrite bookkeeping in isolation
  from any real physics.
  **Deferred, not attempted:** wiring this into `md.rs`'s actual CLI (`@qmmm` flag parsing — still
  just help text, unused, `md.rs:98`) or into `pyo3-gromos`'s `build_simulation()` — this proves the
  *machinery* composes, not end-user CLI/Python wiring (that's P3.7, now unblocked); a dedicated
  `qmmm_total` `Energy` field plus `EnergyFrame`/`.tre`/`PyEnergy` reporting (`special_total` reuse
  is coarse, same accepted-coarseness precedent as P2.8-1's `Contribution.energy` note); P2.8-5
  (refreshing QM charges per step) — untouched, independent of this task.

**2.9 — The full QM/MM → ML/MM pipeline: generate, train, compare**

> Grew directly out of a user question: with real QM/MM (P2.8) and a real, working `--features ml`
> build now both in hand, is the whole "use QM/MM to generate training data, train an ML potential,
> run it as ML/MM, compare against QM/MM" pipeline actually reachable? Answer: yes — every piece
> existed in isolation (`XtbInteraction` for real QM data, `SchNetInteraction` for loading *any*
> real trained TorchScript SchNetPack 2 model, `scripts/export_toy_schnet.py` for the architecture
> and export recipe), the missing part was the glue.

- [x] **`libtorch`/`schnetpack` installed and verified.** `libtorch` 2.11.0 CPU (the exact version
  `tch=0.24.0` pins) + `schnetpack` 2.2.0, into `/tmp/torch_venv`. `python3.13-venv` needed root
  (not available non-interactively here) — worked around with `python3 -m venv --without-pip` +
  manually bootstrapping pip via `get-pip.py`, then the documented recipe (`schnet.rs` module docs)
  otherwise unchanged. `--features ml` builds; `nonbonded::schnet::tests` (3/3) and
  `schnet_burnn_reference.rs` (the real `t_06` methanol QM-zone test) now run for real, not skipped
  — the first time ML/MM has actually executed in this environment.
- [x] **Physical-correctness decision, made before writing any code.** Training target is the QM
  zone's *isolated* energy (`XtbInteraction` with `Embedding::None`), not an electrostatically-
  embedded one. Reasoning: `SchNetInteraction`'s architecture has no environment/charge input
  channel, so training on embedded energies would bake a specific environment's electrostatic
  contribution into a model structurally unable to represent it — exactly the same reasoning behind
  the real BuRNN training target (`zones.rs` docs: `E_burnn = E_QM(inner+buffer) − ...`, environment-
  blind, with embedding handled as a *separate* additive term at inference time). Real ML/MM
  embedding needs a charge-output channel on `Contribution` that P2.7 already deliberately defers —
  a genuine, separate follow-up, not attempted here.
- [x] **Phase 1 — `crates/gromos-md/examples/generate_qm_training_data.rs`.** 10 independent short
  real `xtb`-driven NVE trajectories (91 frames each, same proven-stable setup as `xtb_nve_loop.rs`
  — no classical bond potential or SHAKE needed; xtb's own intramolecular chemistry keeps geometry
  bounded) from randomly perturbed starting geometries/velocities — one long trajectory would only
  explore a single energy shell. 910 real frames written, energies ranging ~−13100 to −13300 kJ/mol
  (matches the scale of the P2.8-2b pinned water oracle). Plain self-describing text output
  (GROMOS-native units throughout — nm/kJ/mol/kJ·mol⁻¹·nm⁻¹ — no unit conversion needed downstream,
  matching `schnet.rs`'s own unconverted pass-through), no new `serde_json` dependency needed.
- [x] **Phase 2 — `scripts/train_qmmm_schnet.py`.** Imports `build_model` directly from
  `export_toy_schnet.py` (same architecture, no duplication) with a modest capacity bump
  (`n_atom_basis=32`, `n_interactions=2` vs the toy's 16/1) since this one needs to actually fit
  data. Manual training loop (not SchNetPack's `Task`/`Trainer`/hydra machinery — simpler to read
  and debug at this scale): combined energy + force MSE loss, `Adam`, train/val split, 400 epochs.
  Exports via the identical `torch.jit.script`/`save` call `export_toy_schnet.py` already proved
  works with this architecture.
- [x] **Phase 3 — `crates/gromos-md/tests/qm_vs_ml_comparison.rs`** (`--features ml`). Real
  `XtbInteraction` vs the trained `SchNetInteraction` on a held-out trajectory (a deliberately
  different, larger starting perturbation than any generator trajectory — not a reserved split).
  Energy and force RMSE checked against generous, explicitly-documented tolerances — the goal is
  proving the *pipeline* works end-to-end, not a chemical-accuracy claim (≈900 frames, a small
  toy-scale model, CPU training on one water molecule).
  **Results, real run:** training loss dropped from 6.2M to ~3-6K over 400 epochs (noisy —
  unbatched single-sample SGD-like updates — but a genuine, large reduction, not a plateau at
  initialization). On 41 held-out frames from a deliberately-different starting perturbation:
  energy RMSE 64.0 kJ/mol against a ~−12800 kJ/mol total-energy scale (≈0.5% relative), force
  component RMSE 668 kJ/mol/nm — both comfortably inside the tolerances (200 / 5000) and, more to
  the point, small enough relative to the energy/force scales actually present to say the trained
  model is really tracking xtb's PES, not just passing a loose bound. Full pipeline (generate →
  train → compare) ran for real, end to end, in this environment.
- **Explicitly not attempted:** real embedding for the trained net (needs `Contribution`'s deferred
  charge channel); training on a QM/MM *zone* rather than a bare vacuum molecule (would need the
  same isolated-target reasoning applied to a real inner/buffer/outer system, e.g. `t_06` — natural
  next step, bigger scope); anything resembling chemical-accuracy validation (would need a much
  larger dataset/model/training budget than this CPU-only sandbox check calls for).

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

**3.1 — `System.from_files()`** ✓ done — `System = Topology + Configuration`, atom-count
validated at construction; `.from_files()`/`.write()`; exposed in `__init__.py`/`.pyi`.

**3.2 — Run a `System` with native parameters** ✓ done — `ImdParameters::default()` +
factories (`nve`/`nvt`/`npt`/`steepest_descent`); `PyInputParameters.from_file()`/`.write()`;
`Simulation(system, params)` two-arg constructor.

**3.3 — Reporters (energy out without files)** ✓ done — `angle`/`dihedral`/`improper` wired
into `PyEnergy` (also fixed the same zeroing bug in the `md` binary's `.tre` writer via a
shared `EnergyFrame::from_energy()`); `sim.run(steps, ene_freq)` batch loop → numpy array;
`EnergyTimeseries` (`to_dataframe()`, `plot()`, `block_average()`).
**Bug found+fixed along the way:** `calculate_bonded_forces_ntf` (`gromos-forces/bonded/mod.rs`)
summed all four bonded terms into one scalar with no per-term breakdown — affected every
`Simulation.step()`, not just Python (plus REMD/GaMD). Fixed by adding per-term energy fields,
populated before combining. Total potential energy was always correct, so no reference test
regressed; caught because `aladip_solvated` (unlike `water_216_box`) has angle/dihedral/improper
terms to show zero.

**3.4 — Notebooks** ✓ done — replaced 3 stale notebooks referencing phantom APIs with
`01_load_and_inspect.ipynb` and `02_short_md.ipynb`; deprecated `md_runners.py`.
**Gap found, not fixed then (closed by 3.6):** factory params didn't expose SHAKE
(`aladip_solvated` diverged to NaN) and `build_simulation` never dispatched to steepest-descent.

**3.5 — M1/M2/SD convenience-path fixes** ✓ done — a full composition-pattern redesign was
audited down to its load-bearing core: only two things were actually broken.
- **M1:** `constraints="none"|"hbonds"|"allbonds"` knob on the `nve`/`nvt`/`npt` factories
  (→ `ntc=1|2|3`); verified empirically (`aladip_solvated` stays finite with `"hbonds"`, still
  diverges with `"none"` as documented contrast).
- **M2:** `temperature` getter's DOF calc was duplicated and could silently disagree with the
  thermostat's; extracted to one shared `compute_total_dof()`.
- **SD:** `build_simulation` now branches on `ntem > 0` to a steepest-descent sequence
  (mirrors `md.rs`'s EM sequence exactly); new composable `SteepestDescent`/`AlgorithmSequence.minimize()`
  building blocks agree with the direct path. Verified: `aladip_vacuum` EM converges and plateaus.
- Also added `sim.volume`/`sim.pressure` getters (documented trap: NVE/NVT `pressure` returns a
  real but not-physically-meaningful kinetic-only number, not zero).

**3.6 — Rust↔Python reference-test parity gap** ✓ done (FEP deferred)
> `py-gromos`'s `REFERENCE_SYSTEMS` (20) vs the full Rust suite (37) had 18 systems missing — real
> Python-API defects (`build_simulation` never wired up features already implemented lower in the
> stack), not a test-authoring gap. Fixed one each: **GENVEL** (`ntivel` parsed but never read,
> also explained `water_1000_spc_gridcell`'s 39% mismatch); **SETTLE/LINCS** (implemented, never
> dispatched — added `ConstraintSelection::from_imd`; regression caught: gating on `imd.nsm` broke
> a 3.5 test since the compositional path never sets it, fixed by gating on actual solvent-atom
> count instead); **Nosé-Hoover(+chain)** (same story, added `push_thermostat()`); **distance/
> position restraints** (`distrest=`/`posresspec=`/`refpos=` kwargs on `Simulation(...)`);
> **triclinic/truncated-octahedron box** (`NTB=-1` transform; also fixed `sim.forces` frame
> rotation — found but left out of scope: gromos-rs's own `.trc` position output already disagreed
> with gromosXX pre-existing, tracked via `POSITION_MISMATCH_SYSTEMS`). Deferred on explicit user
> direction: FEP topology (`ch4_water_fep`, `aladip_vacuum_fep`); the composable
> `AlgorithmSequence`/`resolve_algorithm_sequence` path was untouched. Rust suite 37/37 unchanged
> throughout; Python suite grew 73→100.
> **Follow-up:** `REFERENCE_SYSTEMS` was still missing 7 of 37 active Rust systems (6 steepest-
> descent EM systems + `water_216_box_com_rot`) — the features were already wired, just never added
> to the test list. Added all 7 (6 pass cleanly; `aladip_vacuum_em` has the same known EM
> frame-count off-by-one the Rust suite already `ignore:`s — energies/forces still validate, kept in
> `POSITION_MISMATCH_SYSTEMS` rather than dropped). Python suite: 100→121 (118 passed, 2 documented
> position-mismatch skips). Remaining gap vs Rust is exactly the 2 deferred FEP systems.

**3.7 — ML potential binding for `py-gromos`: name zones, don't count atoms** ✓ done

> Grew out of a design conversation about `SchNetInteraction` (2.6/2.7) and `py-gromos`
> (`crates/pyo3-gromos`): both already existed, but nothing connected a Python user's model to a
> `Simulation`'s zone definition, and `ZonePartition::new(n_atoms, inner, buffer)`
> (`gromos-forces/src/zones.rs`) took raw `&[usize]` — exactly the "renumber the ligand every time
> the topology changes" ergonomics `AtomSelection::from_string` (`gromos-core/src/selection.rs`,
> 2.1) already solved for atom selection generally. This section is that fix applied to zones, plus
> the binding that makes it reachable from Python.

**Two real corrections found while implementing, not just style — the original sketch above (now
struck through in spirit, kept for history) would have shipped a subtle bug or an impossible
requirement:**

1. **`Simulation(..., ml_potential=, ml_region=, ml_buffer=)` kwargs, not a post-hoc
   `add_ml_potential()` method.** Read `build_simulation`'s actual internals: it eagerly calls
   `md_sequence.init()` and runs a full `run_step()` for step 0 before ever returning to the
   caller. A post-construction method would either miss step 0's contribution (already primed
   without it) or need to re-run `init()`/step 0 after insertion — fragile (`AlgorithmSequence`'s
   `push()` only appends, can't insert at the required position right after `Forcefield` anyway,
   and re-`init()`-ing already-initialized algorithms isn't obviously safe). The existing
   `distrest=`/`posresspec=`/`refpos=` kwargs already solve exactly this shape of problem — this
   follows that precedent instead.
2. **The orchestrator's ML term is additive on top of `Forcefield`, not a zone-partitioned
   *replacement* of its classical treatment.** Checked before wiring anything: `Forcefield` (the
   real classical algorithm in every `Simulation`'s `AlgorithmSequence`) has **no**
   `ZonePartition` field at all — unlike the provider-pattern `LjCrfInteraction`, it computes
   classical LJ+CRF for the whole system unconditionally. The original sketch's "constructs a
   `ProviderOrchestrator` entry (classical `LjCrfInteraction` + the ML provider, zone-partitioned)"
   would have registered a *second* classical term through the orchestrator on top of
   `Forcefield`'s already-unconditional one — double-counting every inner-zone pair, silently,
   exactly the failure mode `zones.rs` (assumption A5) exists to prevent. Caught by checking
   `Forcefield`'s struct definition directly rather than assuming the sketch's premise. Fixed by
   registering *only* the ML potential through the orchestrator — an honest, documented **additive
   ML correction term**, not (yet) a rigorous "replace classical with ML for the inner zone"
   scheme. Giving `Forcefield` itself zone-partition awareness is real, separate follow-up work (a
   change to the production classical algorithm), not attempted here — see
   `ml_potential.rs::build_ml_orchestrator_algorithm`'s own doc comment for the full reasoning.
3. **`elements` is caller-supplied, not derived from `Topology`.** The original sketch assumed a
   `Topology` → atomic-number derivation that doesn't exist: `Topology` has only `iac` (a
   force-field type index) and `mass`, no element field. Every real QM/ML provider this session
   got atomic numbers from an external source, never from `Topology` — building that inference is
   real, separate, speculative work. `SchNetPotential(model_path, cutoff, elements: list[int])`
   takes them explicitly.

**What landed:**
- `ZonePartition::from_selections(topo, inner, buffer: Option<&str>)` (`zones.rs`) — thin wrapper
  over `AtomSelection::from_string`, 3 unit tests (name match, with-buffer, bad-selector error).
- `SchNetPotential` (`crates/pyo3-gromos/src/ml_potential.rs`, `#[cfg(feature = "ml")]`) — a
  *recipe* (model_path/cutoff/elements), not eagerly loaded (no topology available yet at Python
  construction time to validate against).
- `resolve_zone_partition(topology, inner, buffer=None) -> (inner, buffer, outer)` — standalone
  pyfunction, always available (doesn't need `ml`), used both by the exit-criterion test and
  internally by `Simulation`'s own construction path (one code path, no duplication).
- `Simulation`'s three new kwargs, threaded through `build_simulation` via a plain
  `MlPotentialSpec` struct (not itself a pyo3 type, so the function signature doesn't change
  across `ml`/non-`ml` builds — only the code that *acts* on `Some` is feature-gated). Pushed
  immediately after `Forcefield` in both the minimization and standard-MD branches, before step 0
  runs.
- `pyo3-gromos`/`py-gromos` both gained `ml = ["...//ml"]` feature forwarding.
  `__init__.py` imports the two new names inside `try/except ImportError` (first time this
  workspace's Python layer needed to handle an optional compiled member); confirmed for real —
  built both with and without `--features ml`, non-`ml` build imports cleanly with `SchNetPotential`
  simply absent (`AttributeError`, not a crash), matching test suite (3 tests) shows 2 skipped, 1
  passed (the kwarg-pairing validation, which works regardless of `ml`).
- **Exit criterion, met on the real `t_06` fixture:** `test_ml_potential.py`'s
  `resolve_zone_partition(topo, "m:a", None)` on the real BuRNN tutorial topology reproduces the
  exact inner zone `zone_partition_reference.rs` already validates (indices 0-5, the 6-atom
  methanol solute — `t_06`'s entire solute *is* the QM zone, so `"m:a"`, "all solute atoms",
  resolves to it exactly). That fixture's *buffer* zone is computed geometrically each step (atoms
  within a radius), not by a static selector string, so it's out of scope for `from_selections`
  and not compared — documented in the test, not silently skipped. A second test builds a real
  `Simulation` (`water_single`, the Rust suite's own Level-1 reference system) with
  `ml_potential=`/`ml_region="1:res(SOL:a)"` and confirms it constructs and steps — the untrained
  toy model from `export_toy_schnet.py`, same honesty tier as every other `ml` test here.
- ~~No libtorch in this environment~~ — done: `libtorch` 2.11.0 CPU + `schnetpack` 2.2.0 in
  `/tmp/torch_venv` (see `schnet.rs` module docs for the recipe).
- **Still open, unrelated to this section's own scope:** a real trained (not toy) model — that's
  exactly what P2.9's QM/MM→ML/MM pipeline produces, usable with this binding once trained on a
  real QM/MM zone rather than a bare vacuum molecule (P2.9's own stated next step).

**Follow-up — `XtbPotential` + `.evaluate()`: a real QM reference from Python, not just Rust.**
`qm_vs_ml_comparison.rs` (P2.9) already runs real `XtbInteraction` vs a trained `SchNetInteraction`
and checks RMSE — but only in Rust. `SchNetPotential` above could be attached to a `Simulation`,
yet there was no way to get a real QM value from Python at all to compare it against. Closed by:
- `PyXtbPotential` (`crates/pyo3-gromos/src/qm_potential.rs`, new module, **not** feature-gated —
  `XtbInteraction` is a subprocess wrapper around the real `xtb` binary, no `libtorch` needed).
  `XtbPotential(work_dir, elements, gfn=2, charge=0, multiplicity=1)`, `.evaluate(positions)` — same
  `Embedding::None` isolated-cluster scope `qm_vs_ml_comparison.rs` itself uses.
- `SchNetPotential.evaluate(positions) -> (energy, forces)` — a standalone, direct call (not wired
  into a `Simulation`'s step loop), so both potentials can be called on the same positions array
  for comparison. Both `.evaluate()` methods build a throwaway `Configuration`/`Topology`/
  `AtomSelection::all`/`Vacuum` `ConfigurationSpatialIndex` internally and call the shared
  `PotentialProvider::contribute` trait method directly — this is a reference/comparison utility,
  not production QM/MM.
- `py-gromos/tests/test_qm_vs_ml_comparison.py` (new) — real `xtb` vs the real trained model from
  P2.9's pipeline (`/tmp/trained_water_schnet.pt`, reused if present) on a held-out water geometry
  distinct from any training trajectory frame; falls back to a fresh untrained toy model (weaker
  claim: finite, right-shaped output only) if no trained model is on disk. **Real result, from
  Python, matching the Rust-side pipeline's own honesty tier:** energy |diff| = 23.8 kJ/mol, force
  RMSE = 1146.6 kJ/mol/nm (both well inside the same generous tolerances `qm_vs_ml_comparison.rs`
  uses — 200 kJ/mol / 5000 kJ/mol/nm — this is a pipeline-correctness check, not a chemical-accuracy
  claim).
- Verified both ways: `maturin develop --features ml` (4/4 relevant tests pass, real RMSE above) and
  plain `maturin develop` (non-`ml` build: `XtbPotential` present and importable, `SchNetPotential`
  correctly absent, `_HAS_ML` is `False`, the 3 `ml`-only tests skip cleanly with a clear reason
  rather than erroring). Full `cargo test -p gromos-forces -p gromos-md -p gromos-core -p
  pyo3-gromos --features ml` and `pytest py-gromos/tests/` both green (one unrelated pre-existing
  failure: `test_energy_timeseries` needs `polars`, not installed in this venv — untouched by this
  work). `cargo clippy -p pyo3-gromos --features ml` clean for both new files (`qm_potential.rs`,
  the `evaluate()` addition to `ml_potential.rs`) — the large pre-existing warning set elsewhere in
  the workspace (loop-index idioms, `.map_or`, etc.) predates this session and isn't reintroduced.

**Stretch, not scheduled — proximity-based zone definition.** `AtomSelection`'s grammar
(`selection.rs`) has no distance operator; the common way people actually define an ML/QM inner
region is "ligand + everything within a cutoff shell," not just by residue name/number. Two ways to
get there, neither started: (a) a manual two-step using the existing `SpatialIndex::neighbor_pairs`
to expand a residue selection into a shell at setup time, entirely in Rust/Python without touching
the grammar; or (b) a `within(radius, spec)` grammar extension backed by the same `SpatialIndex`.
(a) needs no parser change and is the more likely first move if this becomes real.

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
