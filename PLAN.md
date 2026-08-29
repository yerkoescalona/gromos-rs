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
- **`PLAN_ARCHIVE.md`** — the full working notes of finished items (Dim 9a; 2.6–2.9 QM/MM + ML;
  3.5–3.7 py-gromos). PLAN.md keeps one paragraph per finished item with its load-bearing decisions;
  the archive keeps the traps, oracle numbers and rejected alternatives verbatim. `BENCHMARKING.md`
  plays the same role for Priority 4; `CHANGELOG.md` is the per-version record.

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
| 3   | water_216_nve_nobath | 648 | absent MULTIBATH block ⇒ no bath (IMD parser defaults, 3.9 A18) | **PASS** |
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
| 4   | meoh_water_fep | 2862 | charged CH3OH→dummy in 953 SPC water, λ=0.5 (RF self/excluded-pair terms) | **FAIL, ignored** — CRF off by 0.16 kJ/mol at frame 0, LJ exact (2026-08-29) |

**38 of 41 tests pass.** (3 ignored: `aladip_vacuum_fep` — FEP mismatch, bisected 2026-08-29 (see 1.7); `meoh_water_fep` — perturbed RF term for charged atoms; `aladip_vacuum_em` — EM energy frame count off-by-one vs gromosXX)

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

**1.5 — Dim 9: Pairlist from O(N²) to cache-coherent O(N)** — 9a done and O(N) confirmed (P4); 9d/9e/9c/9b open

> Maps to FUTURE.md §Dimension 9. Dispatch is an **enum, not a trait object** (`PairlistAlgorithm
> {Standard, CellList}`, `gromos-core/src/pairlist.rs`): `update()` runs every NSNB steps and an enum
> match inlines to a direct call with no heap allocation. The spatial-index *service* (9e) is a separate
> struct, not a third arm.

> **The math reality that drives the whole test strategy.** Reference tests are NOT byte-identical:
> `test_gromosXX_references.rs` compares every MD step to GROMOS at `ENERGY_REL_TOL = 1e-8`
> (`FORCE_ABS_TOL = 1e-6`). The production pairlist is not order-normalised (`normalize_pairs` is a
> `#[cfg(test)]` helper), so any change to pair *iteration order* changes floating-point summation order
> (~1e-13 relative per step) and grows under MD's positive Lyapunov exponent over a 100-step trajectory.
> **Discipline: prove the invariant in a unit test against a brute-force oracle first; only then touch
> the engine.** A sub-step that changes pair *order* must declare whether it targets bit-identical
> (sort to canonical order) or within-tolerance (measure the margin empirically).

**9a — CellList plumbed, switched on, and fast** ✓ (full notes: `PLAN_ARCHIVE.md` §9a)
- [x] 9a-0 plumbing, zero float change: enum dispatch; `from_imd()` slot mapping faithful to
      `in_parameter.cc:1419-1422` (0=Standard, 1=ExtendedGrid unported→Standard, 2=CellList); wired into
      every caller; Martina bug not reproduced; set-equality vs Standard on all reference systems.
- [x] 9a-1 CellList on, reorder margin measured: **0.0 (bit-identical)** over 100 steps of
      `water_216_box`. Activated only by explicit `ALGORITHM grid_cell`, never a size heuristic; real
      gromosXX 7×7×7 oracle `water_1000_spc_gridcell` added.
- [x] O(N) benchmarked and repaired (P4, `BENCHMARKING.md`): the grid pruned nothing until the cell size
      came from IMD `SIZE` (`grid_cell_pairlist.cc:136-141`) with inter-cell distance pruning; cell-pair
      enumeration is now parallel with order-preserving concatenation. Closes the old "confirm the O(N)
      gap is real" item.

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

**9e — Cell list as spatial-index service** (9a-1 is stable ✓; the `SpatialIndex` *service* already
exists since 2.6 — `gromos-core/src/spatial_index.rs`, used by the QM/ML providers — so this item is now
"rebuild the MD pairlist on it", not "create it")
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
- [x] Parallelize the cell-pair enumeration loop — done in P4 (`BENCHMARKING.md`): ranges over
  `par_chunks` into private `PairlistContainer`s, appended in order (order-preserving, so no re-sort
  was needed; 9a-1's bit-identical margin held). The charge-group-based path got the same treatment.
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
  **2026-08-29: now exercised and failing** — `meoh_water_fep` (charged 54a7 CH3OH → dummies, λ=0.5):
  LJ exact (3e-11), CRF off by 0.16 kJ/mol at frame 0. Reference in place, `ignore`d until fixed.
- **`aladip_vacuum_fep` bisected (2026-08-29)** with single-block `.ptp` variants against the native
  gromosXX (scratch, not committed): unperturbed run exact; `PERTBONDSTRETCH` exact; perturbed angle /
  improper / proper-dihedral **energies are booked into the bond slot** (angle case: angle −0.715, bond
  +0.715 — totals right, `.tre` columns wrong); soft bond/angle/improper energies differ slightly (real
  formula differences); `PERTATOMPARAM` (charged, λ-mixed masses) differs in LJ, CRF *and* kinetic
  (mass mixing convention); `PERTATOMPAIR` differs in LJ. Fix order: energy slots → charged perturbed
  RF (also closes `meoh_water_fep`) → atom pairs → soft bonded → masses.

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

**2.4 — Code quality** ← after 3.9 (390 warnings hide real bugs)
- [ ] Clippy pass: `gromos-forces` (89), `gromos-integrators` (77), `gromos-io` (31), `gromos-core` (15)
- [ ] Replace bare `unwrap()` in non-test code with `.expect("msg")` or `?`
- [ ] Add missing `#[test]`: SHAKE constraints (0 today), improper dihedral unit test
- [ ] Split large files: `nonbonded.rs` (~1500 LOC), `gromos-io/topology.rs` (~1200 LOC)

**2.5 — Stub cleanup** — parked
- [x] `frameout`, `ener`, `rmsd`, `nhoparam`, `ext_ti_ana`, `bar`, `ext_ti_merge` — real implementations
- [ ] `visco`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` — stubs; parked

**2.6 — Dim 12.5: the provider seam (QM/MM/ML)** ✓ landed 2026-08-13 (`PLAN_ARCHIVE.md` §2.6)
- [x] `SpatialIndex` trait + `ConfigurationSpatialIndex` (`gromos-core/src/spatial_index.rs`): a query-based
      neighbour service with per-pair periodic shift vectors, independent of the MD pairlist's
      charge-group/twin-range shape. (The MD pairlist itself is *not* yet built on it — that is 9e.)
- [x] `PotentialProvider` (`gromos-forces/src/provider.rs`): `contribute()` over read-only
      `Topology`/`Configuration` + an `AtomSelection` region → a scattered `Contribution` (energy, sparse
      per-atom forces, virial). Whole-system `&mut Configuration` access was rejected in review because it
      defeats the "arbitrary but *scoped* subset" invariant (FUTURE.md P5). Impls are named `*Interaction`
      (gromosXX's own convention).
- [x] `LjCrfInteraction` — first classical provider; matches real gromosXX `pair_lj` to 1e-8
      (`test_provider_reference.rs`). `SchNetInteraction` (feature `ml`) — TorchScript through `tch`,
      in-process, forces computed inside the scripted graph. **libtorch 2.11.0 exactly** (`tch=0.24.0`'s
      pin; recipe in `schnet.rs` docs). SchNetPack v1 does not `torch.jit.script`, v2 does. Validated by
      finite differences — no gromosXX oracle exists for a neural potential (FUTURE.md P8). Bug caught:
      `neighbor_pairs()` deliberately returns boundary-crossing pairs; the provider must skip them.
- [ ] Deferred: async/cancellable `contribute` for external engines; a typed-units boundary;
      `crates/gromos-forces/src/qmmm.rs` (older file-in/file-out scaffold) untouched, fate undecided.

**2.7 — Electrostatic embedding onto MM atoms (QM/MM + BuRNN)** ✓ Steps 1–5 (`PLAN_ARCHIVE.md` §2.7)
Grounded in Poliak et al. 2025 (mainline GROMOS QM/MM: MM charges polarise the QM density) and
Gómez-Flores et al. 2022 (BuRNN: environment-blind NN over inner/buffer/outer zones). Read from the
tutorial's `mopac.py`: `E_burnn = E_QM(inner+buffer) − E_QM(buffer) − E_solute_vac` (QM−QM, never QM−MM),
so inner/buffer double counting is a *training-time* contract (assumption A5 in `zones.rs`) that the
orchestrator reproduces as a fixed decomposition — it cannot be detected at runtime.
- [x] `Embedding {None, Mechanical, Electrostatic}` as a type-level provider property; an unsupported
      scheme is rejected with an error, never silently dropped.
- [x] `ElectrostaticEmbedding` (`nonbonded/embedding.rs`, Poliak path (c): QM charges + pairwise Coulomb in
      the MD code) — first provider placing forces outside its region. Closed-form, finite-difference and
      real `t_06` checks (1363 MM atoms within `RCUTQM`, E_embed = 508.06 kJ/mol, equal-and-opposite forces).
- [x] `zones.rs`: `Zone {Inner, Buffer, Outer}`, `PairOwner`, `ZonePartition::owner()` — the six-pair-class
      ownership table as an *orchestration* type. `zone_partition_reference.rs`: on real `t_06` every one of
      the 6659 QM-cutoff pairs is claimed exactly once (provider 609 / embedding 6050 / classical 0). The
      originally planned "MM baseline" oracle was found not to exist and was dropped rather than faked.
- [ ] Deferred, named so they are not assumed: a charge-output channel on `Contribution` (it would break
      the additive model — ordered evaluation, not summation; do not add speculatively); link atoms /
      boundary charge redistribution; mutual polarisation (Poliak: QM by MM only, not vice versa).

**2.8 — Wiring: from "the pieces exist" to "a QM/MM run happens"** ✓ done except 4b (`PLAN_ARCHIVE.md` §2.8)
- [x] `LjCrfInteraction` honours `ZonePartition` (`with_zone_partition()`); partitioned + excluded pairs
      reconstruct the full energy to 1e-8 on `t_06`. **Trap recorded:** `Configuration::new()` defaults
      `box_config` to vacuum — set it from the `.cnf` in every bare-`Configuration` test.
- [x] `ProviderOrchestrator` (`orchestrator.rs`): `register(region, provider)` + `evaluate()`; an
      `Embedding::None` provider that places a force outside its region is a hard error (the P2.6 review
      finding, now checked). Transparency (combined == sum of direct calls) holds with three providers on
      the real Poliak water dimer. Note: gromosXX's reported energy includes the per-atom RF self-term
      that `LjCrfInteraction` deliberately excludes.
- [x] Real QM engines as providers, subprocess file-in/file-out (not gromosXX's C-API linking, deferred):
      `XtbInteraction` (`nonbonded/xtb.rs`; `.engrad` parse; `HARTREE_TO_KJMOL`/`BOHR_TO_NM` added to
      `units`) with **real electrostatic embedding, Poliak path (a)** via `pcharge`/`pcgrad` (length unit
      proved empirically to be Bohr, 99.85% match to the analytic monopole term). Path (a) and
      `ElectrostaticEmbedding` (path (c)) are **alternatives, not additions** — documented in both modules.
      `MopacInteraction` (`nonbonded/mopac.rs`): AM1/PM3/PM6/PM7 from one impl, `Embedding::None` only
      (gromosXX's `mopac_worker.cc` returns no MM forces), shared `qm_subprocess.rs`. Oracle lesson: the
      BuRNN archive's MOPAC2016 energies are 11% off MOPAC 23 at fixed geometry → pin oracles from *this
      machine's* engine, never from an archive.
- [x] Virial from `ElectrostaticEmbedding` (Euler identity, FD, real periodic `t_06` at 1e-8).
- [x] `QMLJ` resolved from `in_parameter.cc::read_QMMM` + `modify_exclusions`: it gates a classical LJ
      *supplement* for inner-inner/inner-buffer pairs only; inner-outer LJ is always classical.
      `ZonePartition::qm_lj` + `lj_owner()`.
- [x] External gromosXX QM/MM oracle: Poliak dataset (Zenodo 10.5281/zenodo.14549978, CC BY 4.0);
      `water_dimer_mechst` reproduced to 1.3e-9 (`tests/gromosXX_qmmm_references/`). P2.8-5 (per-step
      charge refresh) superseded: gromosXX has no `MEDC` correction; path (a) makes it moot for xtb.
- [x] `ProviderOrchestratorAlgorithm` (`orchestrator_algorithm.rs`) implements `Algorithm`, pushed right
      after `Forcefield`. **Energy-bookkeeping trap:** contributions go to `special_total` and *overwrite*
      it (a `+=` compounds across steps; `potential_total` is recomputed by `update_potential_total()`), so
      the algorithm claims `special_total` exclusively — coexistence with active position restraints is
      unsupported. NVE drift 0.0032% over 120 real steps (`xtb_orchestrator_sequence.rs`).
- [ ] **P2.8-4b — real MNDO build** to unlock the rest of the Poliak dataset: **blocked on the user's
      registration** (license form + Gitea; not something to do on their behalf). Exit once unblocked:
      build MNDO (`gfortran`), wire gromosXX (`NTQMSW=0`), reproduce the `elst`/`mech` `WATER_DIMER`
      energies as exact oracles.
- [ ] Per-step refresh of static charges for path-(c)-only engines (MOPAC/MNDO/DFTB) — small, not started.
- [ ] Zone-aware `Forcefield`, so a term can *replace* classical physics in a zone instead of only adding
      to it (needed by 3.9 terms that are not purely additive).

**2.9 — The QM/MM → ML/MM pipeline: generate, train, compare** ✓ ran end to end (`PLAN_ARCHIVE.md` §2.9)
- [x] Training target = the QM zone's *isolated* energy (`Embedding::None`): `SchNetInteraction` has no
      environment channel, so an embedded target would bake one environment into a model that cannot
      represent it — the same reasoning as BuRNN's environment-blind target.
- [x] `examples/generate_qm_training_data.rs` (10 short xtb NVE trajectories, 910 frames, GROMOS units)
      → `scripts/train_qmmm_schnet.py` (same architecture as `export_toy_schnet.py`, energy + force MSE)
      → `tests/qm_vs_ml_comparison.rs` (held-out trajectory: energy RMSE 64 kJ/mol on a −12 800 kJ/mol
      scale, force RMSE 668 kJ/mol/nm — a pipeline proof, not a chemical-accuracy claim).
- [ ] Not attempted: embedding for the trained net (needs the deferred charge channel); training on a real
      inner/buffer/outer zone (`t_06`); anything resembling chemical-accuracy validation.

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

**3.5 — M1/M2/SD convenience-path fixes** ✓ done (`PLAN_ARCHIVE.md` §3.5)
`constraints="none"|"hbonds"|"allbonds"` on the `nve/nvt/npt` factories (→ `ntc`, reaching the existing
SHAKE guard with no builder change); one shared `compute_total_dof()` for the thermostat and the
`temperature` getter; `build_simulation` branches on `ntem > 0` to steepest descent (mirrors `md.rs`);
`sim.volume`/`sim.pressure` (pressure is the kinetic-only term outside NPT — documented trap). A full
composition redesign was drafted then parked — see FUTURE.md Addendum 2026-07 and 3.9 below, which
supersedes it.

**3.6 — Rust↔Python reference-test parity gap** ✓ done, FEP deferred (`PLAN_ARCHIVE.md` §3.6)
The 18 systems missing from the Python suite were real binding defects, not test gaps — features
implemented lower in the stack that `build_simulation` never dispatched: NTIVEL=1, SETTLE/LINCS
(`ConstraintSelection::from_imd`, gated on the topology's *actual* solvent atoms because the object path
never sets `imd.nsm`), Nosé-Hoover(+chain), distance/position restraints (`distrest=`/`posresspec=`/
`refpos=` kwargs), truncated-octahedron box + `sim.forces` frame rotation. Python suite 73→121. Known,
tracked divergence: `POSITION_MISMATCH_SYSTEMS` (`aladip_trunc_oct`, `aladip_vacuum_em`) — Rust-side
`.trc` output, not a binding bug. Deferred on user direction: FEP from Python (`ch4_water_fep`,
`aladip_vacuum_fep`) — closed by 3.9 step 1.

**3.7 — ML potential binding for py-gromos: name zones, don't count atoms** ✓ done (`PLAN_ARCHIVE.md` §3.7)
- [x] `ZonePartition::from_selections(topo, inner, buffer)` over `AtomSelection::from_string`.
- [x] `SchNetPotential(model_path, cutoff, elements)` (feature `ml`; a recipe, not eagerly loaded — no
      topology exists at construction time), `resolve_zone_partition()` (always available),
      `XtbPotential(...)`; `.evaluate(positions)` on both gives a real QM reference from Python
      (`test_qm_vs_ml_comparison.py`: |ΔE| = 23.8 kJ/mol vs the P2.9-trained model).
- [x] `Simulation(..., ml_potential=, ml_region=, ml_buffer=)` as kwargs, not a post-hoc method:
      `build_simulation` runs step 0 eagerly and `AlgorithmSequence::push` can only append. **Physics
      caveat kept:** the ML term is *additive* on top of the whole-system `Forcefield` — an honest
      correction term, not a zone-replacing ML/MM scheme (needs 2.8's zone-aware `Forcefield`).
- [x] `elements` are caller-supplied: `Topology` has no element field (only `iac`/`mass`).
- [ ] Stretch: proximity-based zone definition — most likely a `SpatialIndex`-backed two-step expansion
      of a residue selection, not a grammar change.

**3.8 — py-gromos composition audit (2026-08-29)** — findings, each checked in the code

> Question: the Python API keeps getting "customised to the occasion" — a kwarg here, a descriptor
> there — and new terms (QM, ML, restraints, FEP, GaMD/EDS/REMD) keep arriving. What is the scalable
> shape, and how should IMD files and in-memory Python/Rust objects relate?

1. **Three copies of "IMD → algorithm sequence", drifting.** (a) `crates/gromos-md/src/bin/md.rs:987-1298`
   — the most complete (SETTLE, LINCS, Nosé-Hoover(+chain), FEP, restraints, truncated octahedron; GaMD/EDS
   are applied in the loop *after* `run_step`, not as algorithms). (b) `crates/pyo3-gromos/src/
   simulation.rs::build_simulation` — a hand-ported subset of (a) plus feature kwargs. (c) `crates/pyo3-
   gromos/src/algorithm_sequence.rs` presets + `resolve_algorithm_sequence` — Berendsen + SHAKE only: no
   SETTLE, LINCS, Nosé-Hoover, restraints, ML, vacuum or triclinic box (`build_simulation_from_sequence`
   hard-codes `SimBox::rectangular`). **`gromos-md` has no library** (`lib.rs` is a doc comment), so (a)
   is unreachable — that is *why* (b) and (c) exist. Verified: `AlgorithmSequence::new()` has exactly
   these three non-test callers; the other bins (`gamd`, `eds`, `remd`, `mdf`, `md_mpi*`) do not use the
   sequence at all.
2. **Concrete (a)↔(b) divergences today** — the Python engine is not the `md` engine:
   - `four_pi_eps_i` from the topology's PHYSICALCONSTANTS block: honoured in (a), default constant in
     (b). A no-op on the current suite (every reference `.topo` carries 138.9354) but a semantic gap.
   - `parallel_nonbonded = n_atoms > 100` in (a); always serial in (b). The parallel kernels are
     reference-verified at 1e-8 relative, not proven bit-identical to the serial ones.
   - NSM: (a) recomputes it from the coordinate file (and warns); (b) trusts `imd.nsm`.
   - Thermostat DOF: (a) gates solvent-constraint DOF on the *selected* constraint algorithms, (b)/(c)
     use `compute_total_dof` gated on `imd.ntc/ntcs/nsm`; (a) still carries `TODO: solute constraint DOF`.
   - FEP (`@pttopo`, λ, soft-core, `.trg`) exists only in (a).
3. **The descriptor model is closed.** `AlgorithmDescriptor` has 12 fixed variants; adding one term
   touches six places (pyclass, enum variant, `extract_descriptor`, resolver arm, the presets,
   `__init__.py`/`.pyi`) and still leaves (a) and (b) unaware of it.
4. **Force terms arrive as constructor kwargs** (`distrest=`, `posresspec=`, `refpos=`, `ml_potential=`,
   `ml_region=`, `ml_buffer=`): expressible on the file path only, wired one-off inside
   `build_simulation`. The right home exists in Rust — `PotentialProvider` + `ProviderOrchestrator`,
   architecture.md's keystone — but nothing lets a user *compose* providers.
5. **Defaults live in three places** (`300.0`, `-1.0`, `1.0`, `4.575e-4`, `0.5`, `1000` in `md.rs`,
   `simulation.rs`, `algorithm_sequence.rs`), none next to `ImdParameters::default()`.
6. **"Objects" still need an IMD underneath.** Every descriptor field is `Option<T>` resolved with
   `unwrap_or(imd.x)`; `Simulation.from_sequence(topo, conf, params, seq)` requires `params`. The object
   path is a layer *over* the file path, not an alternative to it.
7. **There is no IMD writer.** `gromos-io/imd.rs` has `read_imd_file` only; a run composed in Python
   cannot be handed to gromosXX or reproduced from a file. (3.2's mention of `.write()` on
   `InputParameters` is stale — only `System.write` exists.)
8. **No test that the builders agree.** The Rust suite drives the *binary* (`Command::new(md_bin())`), the
   Python suite drives `build_simulation`; nothing compares (a)↔(b)↔(c). FUTURE.md (Addendum 2026-07)
   named exactly this trap and parked unification "until there's a test that can prove they agreed".
9. Leftovers: `test_basic.py`/`test_advanced_features.py` are placeholders for a past API (31 tests, all
   skipped, referencing types that do not exist); `Vec3`/`Frame` are `f32` while the core and every array
   getter are `f64`; `md_runners.py` (deprecated subprocess runners) and `analysis.py` (stubs) are still
   imported by `__init__.py`; `gromos.pyi` is hand-written with no coverage check.

Worth keeping — right ideas, wrong place: `System = Topology + Configuration`; sequence-as-data with
edit verbs (`insert_after`/`remove`/`replace`) as the escape hatch; the shared `compute_total_dof` and
`resolve_zone_partition`; the `constraints=` knob; `ConstraintSelection::from_imd` gating on actual
solvent atoms.

**3.9 — The composition model: one recipe, one builder, three front-ends** ← **NEXT**

> **Goal, stated as a property the tests can check:** a Python user and the `md` binary can never run
> different physics for the same input, because there is exactly one place that turns "what to run"
> into an `AlgorithmSequence`, and every front-end is only a way to *produce* that "what to run".

> **Reviewed 2026-08-29 from five framework perspectives** — OpenMM; GROMACS (mdp → grompp → tpr,
> gmxapi); HOOMD-blue v3; ASE / i-PI / PLUMED / OpenMM-ML; Polars / pyo3. Full texts:
> `PLAN_ARCHIVE.md` §3.9-review. Every code claim the reviews made was verified in the repo before
> being adopted (list under "Verified while reviewing"). The consensus is folded into the model below
> (changes marked ▲) and summarised in "External review consensus".

**The model**

```
   .imd ──parse──▶                               stage 1              stage 2                 stage 3
   Recipe(...) ───▶  RunRecipe ── build_plan ──▶ Vec<AlgorithmSpec> ── instantiate ──▶ AlgorithmSequence ── start ──▶ Simulation / md loop
   dict ──────────▶  (plain, versioned data)   (self-contained data,   (Box<dyn Algorithm>;         (init + step 0)
                          │                     validated, editable)    reads ONLY plan + RunControl)
                          ├──write──▶ .imd (+ run bundle)      ├──dump──▶ recipe + plan as text (`md @dump`)
                          └──echo───▶ <run>.recipe.toml + Diagnostics (what was defaulted / passed through / coerced)
```

- **`RunRecipe`** (`crates/gromos-run/src/recipe.rs` — ▲ a new *library crate*, see A4): plain
  `Serialize + Deserialize + PartialEq` data with ▲ `version: u32`, grouped by orthogonal concern —
  `control` (dt, steps, t0, seed, initial velocities/temperature, COM removal);
  `boundary` (vacuum / rectangular / triclinic / truncated octahedron);
  `forcefield` (cutoffs, RF, pairlist algorithm/size/frequency, NTF flags, `four_pi_eps_i` policy,
  ▲ `restraints` — position/distance, GROMOS-native, **owned by `Forcefield` as today** (A16) —
  and ▲ `terms: Vec<TermSpec>` — only `PotentialProvider` implementors);
  `constraints` (solute none/hbonds/allbonds × SHAKE/LINCS; solvent none/SHAKE/SETTLE/LINCS;
  tolerance, NTISHK); `ensemble` (▲ `thermostat: Option<…>` with **all** baths, not `temp_bath[0]`;
  `barostat: Option<…>` with the tensor form; no `-1.0`-means-off sentinels); `integrator`
  (leap-frog | steepest descent + EM parameters); `perturbation` (λ, soft-core, pttopo); `outputs`
  (write frequencies — kept for the IMD round trip and `md.rs`; `Simulation` never writes files
  implicitly); ▲ `execution { parallel: Auto|Serial|Parallel, threads }` — a platform concern, not
  physics, **not written to IMD** (documented non-IMD group); ▲ `passthrough` — an *allowlist* of inert
  raw blocks: a physics-bearing block gromos-rs does not model (GAMD, EDS, LOCALELEV, QMMM, …) is an
  error unless `allow_passthrough=[…]` names it (A17). ▲ **Defaults:** `RunRecipe::default()` is
  *derived* from `ImdParameters::default()` (one table, not two kept equal by a test) and means "what
  gromosXX does when the block is absent" — verified by running gromosXX on `to_imd(default())` (G7).
- **`TermSpec`** — the providers: `SchNet {model, cutoff, elements, region, buffer, …}`, `Xtb {…}`,
  `Mopac {…}`, later `LocalElevation`, … `#[serde(tag = "kind", deny_unknown_fields)]` so Python's
  `Term(kind, **params)` *is* this data and a typo'd parameter is an error, never a dropped key.
  ▲ Every engine term carries `coupling: Delta | Replace { qm_lj }` — `Delta` (additive, what exists)
  is the only value `instantiate` accepts until 2.8's zone-aware `Forcefield` lands; `Replace` is data
  now and is rejected with "requires zone-aware Forcefield (2.8)", never silently treated as additive
  (A8) — plus `units` (per-engine length/energy factors, gromosXX `QMUNIT` defaults) and provenance
  (model checksum, engine version) in the bundle. ▲ `Contribution`/the orchestrator gain a per-term
  energy slot keyed by registry name, so `.tre`, `run()` and parity tests see each term (G10).
- **`AlgorithmSpec`** — the sequence-as-data the escape hatch edits: `RemoveCOM {…}`,
  `Forcefield {…, restraints}`, `Orchestrator {terms}`, `LeapFrogVelocity`, `Thermostat {…}`,
  `LeapFrogPosition`, `Shake {…}`, `Settle`, `Lincs {…}`, `SteepestDescent {…}`,
  `TemperatureCalculation`, `PressureCalculation {…}`, `Barostat {…}`, `EnergyCalculation`.
  ▲ One spec = one algorithm (no combined `LeapFrogIntegrator`); ▲ **fully resolved** — DOF, virial
  type, NSM, `four_pi_eps_i`, restraint paths and the `parallel` decision are *values* written by
  `build_plan`, never `Option` + "look it up in the recipe"; serde + `PartialEq` + `version` like the
  recipe. Replaces pyo3's closed `AlgorithmDescriptor`.
- **One builder, ▲ three stages** (`crates/gromos-run/src/{plan,build}.rs`):
  `build_plan(&RunRecipe, &Topology) -> Result<Vec<AlgorithmSpec>, RunError>` encodes the GROMOS
  step order once; ▲ `validate_plan(&[AlgorithmSpec]) -> Result<(), RunError>` checks the ordering
  invariants, declared as registry data (`Orchestrator.must_follow = Forcefield`;
  `Barostat.requires = PressureCalculation`; `EnergyCalculation` last; `SteepestDescent.excludes =
  {LeapFrog*, Thermostat, TemperatureCalculation}`; exactly one `Forcefield`; at most one
  `Orchestrator`; `Barostat` + a provider that reports no virial is rejected) and runs on *every* plan,
  edited or not; ▲ `instantiate(&[AlgorithmSpec], &RunControl, &Topology, &Configuration) ->
  Result<AlgorithmSequence, RunError>` **takes no `&RunRecipe`** — an edited plan is exactly what runs
  (G8); ▲ `start(&mut seq, …)` = `init()` + step 0 (NTISHK moves positions here, so it is physics),
  shared by `md.rs` and `Simulation` and reusable by a later `Simulation.set_positions()`.
- **`prepare_system(&RunRecipe, TopologyData, CoordinateData, &RunInputs) -> Result<(Topology,
  Configuration, Diagnostics), RunError>`** — everything both `md.rs` and `build_simulation` do
  *before* the sequence (NSM from coordinates, perturbation-topology merge, truncated-octahedron
  transform, initial velocities, validation). Lifted once. ▲ `Diagnostics` (blocks absent → defaulted,
  passed through, values coerced) is returned by `from_imd`, `prepare_system` and `build_plan`;
  `md` writes `<run>.recipe.toml` next to `.tre` (GROMACS's `mdout.mdp`) and `Simulation.recipe` /
  `.diagnostics` expose the same; ▲ `md @input x.imd @dump` prints recipe + plan and exits — the
  Rust-side A/B.
- **Front-ends.** `RunRecipe::from_imd(&ImdParameters) -> (RunRecipe, Diagnostics)` /
  `to_imd(&Topology)` (▲ needs the topology to synthesise MULTIBATH `LAST` and FORCE `NRE` for a
  Python-born recipe) — the file path with GROMOS's bundled `ntc`/`ntcs`/`ntf` semantics kept
  (FUTURE.md's parked ensemble⟂constraints split is *not* reopened). ▲ Factories `RunRecipe::nve/nvt/
  npt/minimize(...)` live in the library, not the binding (the binding owning defaults is how Polars'
  eager/lazy defaults drifted). `gromos.Recipe` — a pyo3 wrapper, ▲ immutable by construction
  (`with_term(...)`, `with_control(dt=…)` return new objects; a `#[getter]` returning a `Vec` makes
  `recipe.terms.append` a silent no-op), `from_dict/to_dict` through ▲ `pythonize` (not JSON strings:
  NaN/inf, errors name the kwarg), pickling via serde, `to_dict` same-version-only. `AlgorithmSequence
  .from_recipe(system, recipe)` — the plan, editable ▲ by index/identity as well as by name (a name is
  ambiguous once a kind repeats).
- **Escape hatch stays.** `Simulation(system, recipe, sequence=edited)` validates and instantiates an
  edited plan with the recipe's `RunControl`; ▲ `Simulation.plan` is a frozen snapshot, never a live
  handle (OpenMM's `reinitialize()` FAQ, HOOMD's `SyncedList` bugs).
- ▲ **Errors:** one `RunError { MissingFeature{term, feature}, UnknownKind, UnknownParameter{kind,
  field}, InvalidPlan{reason}, Diagnostics, Io }` in the library; `gromos.exceptions` mirrors it with
  one `From` in the binding (today: ~50 ad-hoc `PyErr::new::<…>` sites and `Result<(), String>`).
- **Providers compose** through `ProviderOrchestrator` exactly as P3.7's ML term does today.

**External review consensus (2026-08-29)** — verdict per decision; `Modify` means "keep, with the
condition folded in above"

| Decision | OpenMM | GROMACS | HOOMD | ASE/i-PI | Polars | Outcome |
|---|---|---|---|---|---|---|
| D1 recipe as grouped data, defaults once | Modify | Modify | Agree | Modify | Modify | kept; defaults *derived*, presence-aware parsing, `execution` group outside IMD, restraints ≠ terms |
| D2 builder as plan → instantiate | Agree | Agree | Modify | — | Modify | kept; three stages, plan self-contained, `instantiate` takes no recipe |
| D3 serde bridge + exhaustive registries | Modify | Agree | Modify | Agree | Modify | kept; `deny_unknown_fields`, `pythonize`, `version`, registry carries (name, type, default, available) |
| D4 terms as additive providers | Modify | Modify | Modify | Disagree | — | changed: `coupling` declared per term, per-term energy slot, restraints stay `Forcefield`-owned |
| D5 IMD front-end + writer + bundle | Modify | Agree | Agree | Modify | Modify | kept; `to_imd(&Topology)`, allowlisted passthrough, provenance, factory outputs also gromosXX-checked |
| D6 one-cycle deprecation shims | Agree | Agree | Agree | — | Modify | kept; shims in Python, `pytest.warns`, migration table in the same PR |
| D7 editable plan | Modify | Modify | Modify | Agree | Modify | kept; `validate_plan` mandatory, index/identity edits, frozen snapshot |
| D8 guards + step order | Agree | Agree | Agree | Modify | Modify | kept; determinism baseline first, G2 pins `Serial`, Rust-side A/B at step 2, step 5 needs a physics oracle |

*Unanimous:* the skeleton — each reviewer mapped it onto their framework's *mature* shape
(`createSystem → System → Context`; `grompp → tpr → mdrun`; HOOMD v3 declarative → attach;
`DslPlan → IR → Executor`); `#[serde(deny_unknown_fields)]`; parity test before any refactor.
*Four of five independently found the same two gaps:* "term" was not one contract (restraints are
`Forcefield` fields; the orchestrator overwrites `special_total`, so two terms could not coexist), and
an editable plan had no validator. *Splits, resolved:* where restraints live → `Forcefield`-owned
(zero physics-path change; the restraint oracle is a single system, `nacl_1water_distres`); whether
`outputs` belong in the recipe (OpenMM: not the model; HOOMD: instantiate the writers) → kept for the
IMD round trip and `md.rs`, `Simulation` never writes implicitly; crate placement (Polars: split now)
→ adopted; step 5's metric (ASE: a file count proves wiring, not physics) → both required.

*Changes adopted, each traceable to a step:* C1 plan self-containment (G8, step 2) · C2 `validate_plan`
(G9, step 2) · C3 restraints ≠ terms + per-term energy (G10, step 2, before any second `TermSpec`
variant) · C4 `coupling`/`units`/provenance on engine terms (steps 3, 5) · C5 presence-aware parser +
`Diagnostics` + derived defaults (step 2; fixes the verified absent-MULTIBATH bug) · C6
`deny_unknown_fields` + `pythonize` + `version` (steps 2–3) · C7 `execution` group; G2 pins `Serial`;
same-path rerun baseline (step 0) · C8 allowlisted passthrough (step 2) · C9 `gromos-run` crate
(step 1) · C10 `RunError` / `gromos.exceptions` (step 1) · C11 shared `start` stage (step 1) · C12
shims in Python, `mypy.stubtest` for G5, `md @dump` (steps 2–3) · C13 step 5 = wiring count **and** a
physics oracle.

*Verified while reviewing* (claims checked in the code, not trusted): an `.imd` **without MULTIBATH
silently runs a Berendsen bath** (300 K, τ = 0.1) in *both* engines — `read_imd_file` starts from
`ImdParameters::default()` (`imd.rs:329`) and only replaces `temp_bath` inside the MULTIBATH arm
(`imd.rs:422-500`); the only bath-less references are the two EM systems, which skip the thermostat,
so the suite cannot see it. `parse_f64`/`parse_i32`/`parse_usize` coerce garbage to 0 (`imd.rs:1072-
1082`). Every path reads `temp_bath[0]` only. `XtbInteraction.virial == Mat3::ZERO` (`xtb.rs:359`).
The orchestrator sums one scalar (`orchestrator.rs:96`). xtb uses fixed filenames in `work_dir` and
no timeout (`xtb.rs:292-295`). SchNet computes in f32 (`schnet.rs:201`). The parallel kernels reduce
through rayon `fold/reduce` whose tree depends on work-stealing (`innerloops.rs:576-607`), so
`np.array_equal` under `parallel` is not guaranteed even for the same path twice. The only existing
builder-vs-builder test compares `.names` (`test_gromosXX_references.py:598`). `HARTREE_TO_KJMOL`
(8 s.f.) and `BOHR_TO_NM` (12 s.f.) are of different vintages. `dlamt` is parsed and never applied.
`pyo3-gromos`'s `ml` feature forwards only `gromos-forces/ml`. `__init__.py` hardcodes `0.1.0`
while the workspace is `0.0.26`.

**Assumptions — exposed, each with what happens if it is wrong and how it is checked** (▲ = revised
or added after review)

| # | Assumption | If wrong | Check |
|---|---|---|---|
| A1 | The three paths *intend* identical physics for the same IMD; every A/B difference is a defect or a missing feature, never a "different but fine" convention. | The difference becomes an explicit recipe option, never an implicit per-path behaviour. | Step 0's divergence table: every row gets a resolution. |
| A2 | Serial and parallel nonbonded kernels are **not** assumed bit-identical (reference-verified at 1e-8 only). | Python inherits `parallel: Auto` and shifts within 1e-8. | The serial-vs-parallel delta is measured once (n ≥ 3, std) and recorded in BENCHMARKING.md. |
| A3 ▲ | Between front-ends of the *same* build, "bit-identical" means `==` on energies/forces/positions after N steps **with `execution.parallel = Serial`** — the parallel reduction tree is work-stealing-dependent (verified), so parallel results are compared at the reference tolerances only. A same-path rerun baseline is measured before any cross-path comparison. | Serial not deterministic → a real bug, not a tolerance. | `np.array_equal`, never `allclose`, in `test_front_end_parity.py`; baseline test `same_path_twice`. |
| A4 ▲ | The recipe/builder live in a **new `gromos-run` library crate** (serde non-optional; no clap/env_logger/mpi/cudarc), depended on by the `gromos-md` binaries and `pyo3-gromos`; the `gromos` facade does not depend on it. | Fold back into `gromos-md` — but only if the split fails, not for convenience. | `cargo build --workspace`; `cargo tree -p pyo3-gromos -i gromos-run`; no `process::exit` in `gromos-run`. |
| A5 | An `.imd` describes a run only together with its auxiliary files (topo, cnf, pttopo, posresspec, refpos, distrest, model path). The recipe carries those as `RunInputs`; the full serialisation is a bundle in the reference systems' existing `input.toml` shape, ▲ plus provenance (model checksum, engine version). | — the shape already exists in 37 directories. | `from_bundle(to_bundle(dir))` round trip on every reference directory. |
| A6 | GaMD/EDS/REMD/MPI/CUDA binaries stay out of scope: they don't use `AlgorithmSequence`, aren't reachable from Python, and GaMD/EDS in `md.rs` act *after* `run_step` (and re-parse the input out-of-band, `md.rs:646-668`). ▲ Their blocks are physics-bearing and therefore *not* silently passed through (A17). | They become `AlgorithmSpec` variants when P1.9 turns them into algorithms. | `grep -L AlgorithmSequence crates/gromos-md/src/bin/*.rs` lists exactly those bins. |
| A7 | The IMD writer must emit what `read_imd_file` *and* gromosXX accept. | A round trip passes while gromosXX rejects the file. | Step 2 runs the gromosXX native build on every `to_imd` output of the suite ▲ and on every factory output. |
| A8 ▲ | Terms are additive **and say so**: `coupling: Delta` is the only accepted value until 2.8; `Replace` exists as data and is rejected at `instantiate`. A `Delta` engine term is an honest correction, not QM/MM. | A term needs to remove classical pairs → that *is* 2.8, not a 3.9 change. | `instantiate` rejection test; `TermSpec` docs. |
| A9 | Existing Python users (`InputParameters` factories, `Simulation(..., distrest=…)`, `AlgorithmSequence.nvt(...)`) get **one deprecation cycle** with warning shims, then removal; `__version__` is 0.1.0, no stability promise yet. ▲ Shims live in Python, print their `Recipe(...)` equivalent, and ship with a migration table. | Keep the shims longer. | `pytest.warns(DeprecationWarning)` tests (the pyproject ignores the category globally); notebooks and examples migrated in the same PR. |
| A10 | Adopting `md.rs`'s rules in the shared code (`four_pi_eps_i` from the topology, NSM from coordinates) moves Python *toward* gromosXX and changes no reference result (all topologies carry 138.9354; NSM is consistent in all 37). | A Python-suite failure in step 1 means that system was wrong before — fix forward, document. | Python suite green after step 1. |
| A11 | One DOF formula (constraint-aware from *actual* solvent + solute constraints) replaces the three, ▲ and it is written into the plan by `build_plan` (from the plan's `Shake/Settle/Lincs` entries), never re-derived at `instantiate`. Thermostat DOF changes energies, so the NVT references are the guard. | `water_216_nvt*` fail → the formula is wrong; revert to `md.rs`'s and record why. | Rust 37/37 and Python suite after step 1; G8. |
| A12 ▲ | serde tagged enums are the Rust↔Python data bridge, through `pythonize` (not JSON strings), with `deny_unknown_fields` on every variant and per-struct `#[serde(default)]` decisions (lenient on recipe groups, strict on `TermSpec`); `to_dict` is same-version-only. | Generate per-field getters with a macro instead. | `Recipe.from_dict(r.to_dict()) == r` on all 37; a typo'd kwarg fails. |
| A13 ▲ | `gromos.pyi` stays hand-written; drift is caught by `python -m mypy.stubtest gromos` (needs `#[pyo3(signature = …)]` on every method) plus the registry coverage test. | Adopt `pyo3-stub-gen` later. | `stubtest` clean in CI (step 3). |
| A14 | The ML feature stays compile-time: `TermSpec::SchNet` always exists as data; `instantiate` returns `RunError::MissingFeature` without `--features ml`; ▲ registry entries carry `available: bool` and `gromos.build_info()` says which features were compiled; the `ml` feature is forwarded through `gromos-run`, `pyo3-gromos` and `py-gromos`. | — | Non-`ml` build test asserts the error variant; CI builds both feature sets. |
| A15 | The `md` binary's console output and file layout are not part of parity; only energies, forces, positions and dH/dλ are. | — | Reference suite unchanged. |
| A16 ▲ | Position/distance restraints stay `Forcefield`-owned parameters (`AlgorithmSpec::Forcefield { restraints }`), **not** `TermSpec`s: moving them is a physics-path change guarded by one reference system. | A second restraint oracle appears → revisit; until then "one variant + one arm" is a claim about providers only. | Restraint results bit-identical before/after step 2 on `nacl_1water_distres` and the posres EM systems. |
| A17 ▲ | `passthrough` is an allowlist of inert blocks; an unknown or physics-bearing block (GAMD, EDS, LOCALELEV, QMMM, …) is an error unless explicitly allowed. A7 cannot catch this class (gromosXX reproduces *its own* energies while the recipe runs something else). | — | `from_imd` on each such block errors; the allowlist is a test fixture. |
| A18 ▲ | The IMD parser becomes presence-aware (`Option` per optional block: MULTIBATH, PRESSURESCALE, PERTURBATION, DISTANCERES, POSITIONRES, ENERGYMIN) and unparseable numbers are diagnostics, not zeros. "Absent block" ≠ "default values": the recipe default is *defined* as gromosXX's absent-block behaviour. | A reference regresses → that system relied on the silent default; fix forward. | New `water_216_nve_nobath` reference (step 0) fails before C5 and passes after. |
| A19 ▲ | ML terms are not bit-reproducible (f32 tensors, libtorch thread pool); parity for an ML term is at a stated per-term tolerance, recorded in G5's table, never `==`. | — | `test_front_end_parity.py` marks ML rows with their tolerance. |

**Drift guards — what keeps the model honest after the refactor** (all run under `just test`, `just
lint` and `cd py-gromos && uv run pytest`; ▲ = revised or added after review)

- **G1 — one code path under both oracles.** The `md` binary calls `prepare_system` + `build_plan` +
  `validate_plan` + `instantiate` + ▲ `start`; so does `Simulation`. The Rust suite (drives the
  binary) and the Python suite (drives the binding) then test the *same* builder from two sides. This
  is the structural fix; everything below is belt-and-braces.
- **G2 — front-end parity** (`py-gromos/tests/test_front_end_parity.py`): every reference system
  through `Simulation(system, Recipe.from_imd(...))`, through `Simulation(system, Recipe.<factory>(...))`
  where a factory can express it, and through `Simulation(system, recipe,
  sequence=AlgorithmSequence.from_recipe(...))` → `np.array_equal` on energies, forces and positions
  after NSTLIM steps, ▲ with `execution.parallel = Serial` pinned and a `same_path_twice` baseline.
- **G3 — round trips**: `from_imd(parse(to_imd(r))) == r` on all 37 imds ▲ plus a byte-level check
  that passthrough blocks are re-emitted verbatim; `from_dict(to_dict(r)) == r`;
  `from_bundle(to_bundle(dir))` on every reference directory; the gromosXX native build accepts every
  `to_imd` output — round trips *and* factory outputs — and reproduces its own `expected/` (A7).
- **G4 — exhaustive registries at compile time**: `TermSpec::registry()` and `AlgorithmSpec::
  registry()` are built with an exhaustive `match` over `Self` (`#![deny(clippy::wildcard_enum_match_
  arm)]` in the crate) — adding a variant without an entry ▲ `(name, parameters with types and
  defaults, available, ordering constraints, one *runnable* example)` is a compile error.
  `gromos.terms()` / `gromos.algorithms()` expose them.
- **G5 ▲ — stubs and coverage**: `python -m mypy.stubtest gromos` replaces the substring test;
  `test_every_kind_has_a_parity_case` asserts every kind is exercised by at least one reference system
  or unit fixture built from its registry example (the `restraints + SchNet` combination is a required
  row; ML rows carry their tolerance, A19).
- **G6 ▲ — no second builder / no second reader**: `just lint` greps for `AlgorithmSequence::new()`
  outside `crates/gromos-run/src/` (excluding `#[cfg(test)]` modules), for `.push(Box::new(` under
  `crates/pyo3-gromos/`, for `parse_file(`/`read_imd_file(` outside `gromos-run` and `gromos-io`, and
  for `process::exit` inside `gromos-run` — all must return nothing.
- **G7 ▲ — defaults in one place**: `RunRecipe::default()` is derived from `ImdParameters::default()`
  (no second table); gromosXX runs `to_imd(RunRecipe::default())` and every factory output and its
  energies are what the recipe produces; the literal grep for `300.0`, `-1.0`, `4.575e-4`, `0.5` is
  scoped to `crates/pyo3-gromos/src` and `crates/gromos-md/src/bin`.
- **G8 ▲ — plan self-containment**: `instantiate` has no `&RunRecipe` parameter (compile-time), and
  for every reference recipe `r` and every recipe field `f` that `build_plan` reads,
  `instantiate(build_plan(r))` is bit-identical to `instantiate(build_plan(r))` after `r.f` is changed
  *and the plan re-derived by hand to the same values* — i.e. nothing at stage 2 reaches past the plan.
  DOF, virial, `four_pi_eps_i`, NSM and `parallel` are asserted present in the plan.
- **G9 ▲ — plan validation**: one test per ordering invariant that a violating plan (built by
  `insert_after`/`remove` from a valid one) is rejected by `validate_plan` with a named reason; the
  invariants are registry data, so a new kind without them fails G4.
- **G10 ▲ — per-term energies are visible**: each term's energy appears in `.tre`, in `run()`'s
  columns and in G2's comparison; a compensating-error test (two terms whose totals cancel) must still
  fail parity.

**Order of work** — each step is a commit boundary with green suites; nothing later starts until the
gate before it holds. Sizes are estimates of focused work (S ≈ ½ day, M ≈ 1–2 days, L ≈ 3+ days);
▲ marks review-driven additions.

- [x] **Step 0 — Measure before moving (S).** ✓ 2026-08-29. `py-gromos/tests/test_front_end_parity.py`
      v0: every `REFERENCE_SYSTEMS` entry through `Simulation(topo, conf, params, **kw)` vs
      `Simulation.from_sequence(topo, conf, params, AlgorithmSequence.from_parameters(topo, params))`,
      `np.array_equal` after NSTLIM steps. **Result: 27/37 bit-identical; 10 divergences**, all
      `xfail(strict=True)` entries in `EXPECTED_DIVERGENCE` naming the missing feature (rows below) —
      strict, so that when the builder lands they must be *deleted* (a silent pass is drift the other way).
      - [x] ▲ `same_path_twice` baseline: bit-identical on all 37 systems (the Python path is serial).
            Binary, `scripts/kernel_determinism.py` (n = 3 per thread count): run-to-run bit-identical
            at 1 and at 8 threads on `water_216_box` and `ch4_water_fep`; 1 vs 8 threads differ in the
            last printed digit (max rel Δ 1.6e-13 / 2.7e-13). Pinning the thread count therefore already
            gives exact same-build parity; `Serial` stays the documented stronger setting (A3 stands).
      - [x] Divergence table filled (below): 10 A-vs-C rows + 2 engine-vs-gromosXX rows.
      - [x] Serial-vs-parallel delta recorded in BENCHMARKING.md ("Kernel determinism") with n.
      - [x] ▲ `water_216_nve_nobath` generated with the CODATA-patched `BUILD/program/md` — the only
            gromosXX build that reproduces the committed references bit-for-bit (`BUILD_native`/
            `BUILD_omp` differ in the last digit; `scripts/regen_gromosXX_references.py` already points
            at `BUILD`). gromosXX logs "Adding a bath, no temperature coupling" and its energies equal
            `water_216_box`'s to the digit; gromos-rs prints "Setting up thermostat: Berendsen" and
            fails at frame 0 (E_total −4270.69 vs −4270.19). Ignored in the Rust suite, strict-xfail
            (`EXPECTED_ENGINE_FAILURES`) in Python, until step 2.
      - [x] ▲ `test_steepest_descent_via_algorithm_sequence` compares energies and positions — and
            passes: the two EM paths agree exactly.
      - [x] Rust-side A/B is not possible yet (no library) — arrives at step 2 (`@dump`,
            `cargo test -p gromos-run`).
      **Gate met 2026-08-29:** Python suite 185 passed / 2 skipped / 13 xfailed; every parity failure is
      a strict xfail with a named feature and a table row; `water_216_nve_nobath` fails for the
      documented reason in both suites.

- [x] **Step 1 — Lift, zero behaviour change (M).** ✓ 2026-08-29. New crate `crates/gromos-run`
      (A4): `prepare_system`, `build_sequence_from_imd(&ImdParameters, &Prepared, &RunInputs,
      &RunOptions)`, `start`, `RunInputs {pttopo, posresspec, refpos, distrest}`, `RunOptions
      {parallel}`, `RunError`, `total_dof`, `periodicity_of`. Bodies moved from `md.rs:414-572` +
      `987-1298`; the binary's `println!`/`process::exit` sites became `RunError` variants and
      `PrepareNote`/`BuildSummary` data that `md.rs` prints (`md.rs` 2060 → 1359 lines). `md.rs` and
      `build_simulation` both call it; `ml` forwarded through `gromos-run` → `pyo3-gromos` → `py-gromos`.
      - [x] Divergences resolved *toward `md.rs`* (A10/A11): `four_pi_eps_i` from PHYSICALCONSTANTS
            (no oracle — every reference topology carries the default; propagation is by inspection);
            NSM from coordinates (`nsm_comes_from_the_coordinate_file`); `ParallelPolicy::Auto`
            (`parallel_policy_resolves_like_the_binary`); one DOF formula with solute constraints
            (`total_dof_counts_solute_constraints`; the survey found no reference with NTC>1 *and* a
            live thermostat, so no reference moved). Tests: `crates/gromos-run/tests/prepare_and_build.rs`.
      - [x] `build_simulation_from_sequence` goes through `prepare_system` (vacuum/triclinic
            preparation; its `Forcefield` still comes from the descriptor resolver, so
            `aladip_trunc_oct` stays an xfail until step 2).
      - [x] FEP reaches Python — not through a kwarg but through `Topology.apply_perturbation(path)`
            (the `.ptp` merge became `Topology::apply_perturbation` in gromos-core, shared with
            `md @pttopo`). `ch4_water_fep` added and **passes** from Python; `aladip_vacuum_fep` strict
            xfail (the Rust suite's known mismatch); `perturbation_topology_follows_ntg` covers the
            NTG=0 rules.
      - [~] `gromos.exceptions` **deferred to step 3**: step 1 keeps the builtin exception types the
            binding raised before (`run_err` maps `RunError` variants onto IOError/ValueError/
            RuntimeError) — zero behaviour change was the point of this step; the typed hierarchy
            arrives with the recipe's new error kinds.
      - [x] Found and fixed on the way: the descriptor path built a *serial* `Forcefield` while the
            binary is parallel above 100 atoms — the six >100-atom systems lost exact A/C parity the
            moment path A inherited the binary's policy (A2, exactly as predicted). The resolver now
            applies `ParallelPolicy::Auto` and the topology's `four_pi_eps_i` too.
      **Gate met 2026-08-29:** Rust 37 passed / 3 ignored; Python 208 passed / 16 skipped / 18 xfailed
      (no XPASS: the step-0 xfails name features the *descriptor* path still lacks — they close in
      steps 2–4, not here); `grep -c 'md_sequence.push' crates/pyo3-gromos/src/simulation.rs` = 0;
      `cargo tree -p pyo3-gromos -i gromos-run` shows no cycle; no `process::exit`/`println!` in
      `gromos-run`; clippy clean on the new crate.

- [x] **Step 2 — `RunRecipe` + `AlgorithmSpec`, IMD both ways (L).** ✓ 2026-08-29. Deviations from the
      sketch are marked ✗ and carried, not hidden.
      - [x] `ImdParameters` **lossless and strict** (gromos-io): every field of every modelled block
            is now parsed — `ntirtc`, `ntisti`, `ntwse/ntwg/ntwb`, `ntpp`, `ncyc`, `ntcs0`, MULTIBATH
            DOF sets, PRESSURESCALE `couple`/SEMI, NONBONDED lines 3–7, multi-line TITLE, NTWV/NTWF as
            *frequencies* (they were `bool`s — a force-write frequency of 5 came back as 1); garbage
            numbers are errors naming block and token (`imd.rs` used to coerce to 0); `parse_imd_str`
            for in-memory text. **Absent means absent** (A18): `temp_bath` defaults to empty (no bath)
            and `nscm` to 0 (gromosXX `parameter.h`: `skip_step 0`) — the second silent-default bug,
            found the same way as the first. ✗ Presence is read from `raw_blocks` rather than
            `Option<Block>` fields, and `raw_blocks` stays a `HashMap` (the writer sorts passthrough)
            — same guarantees, no struct upheaval.
      - [x] `gromos_io::imd_write::{write_imd, write_imd_file}` — every modelled block regenerated in
            GROMOS layout, shortest-round-trip floats, passthrough verbatim. **Two defects only a real
            gromosXX read caught (A7 earned its keep):** fixed-width columns fused `1e-10` with its
            neighbour (`NA2CLC: 20.0000000001`), and an "off" DISTANCERES block invented for a run
            without one is *validated* by gromosXX (`DIR0 > 0`, `TAUDIR > 0`) — optional blocks are now
            written only when read or non-default.
      - [x] `gromos-run/src/recipe.rs`: `RunRecipe` (`version = 1`; groups `system / control /
            boundary / forcefield (+restraints, +terms) / constraints / ensemble / minimisation /
            perturbation / outputs / execution / passthrough`), serde with `deny_unknown_fields`
            everywhere and `default` on groups, `Default` *derived* from `ImdParameters::default()`
            (G7), `from_imd_with(imd, &PassthroughPolicy) -> (RunRecipe, Diagnostics)`, `to_imd()`,
            `to_imd_string(n_atoms)`, TOML/JSON, factories `nve/nvt/npt/minimize` delegating to
            `gromos-io`'s. `Diagnostics` names every absent optional block and its meaning, every
            passthrough, and required blocks gromosXX would refuse to run without. ✗ `to_imd` takes an
            atom count, not a topology (MULTIBATH DOFSET / FORCE NRE are the only atom-dependent
            lines). ✗ `RunBundle` + provenance (A5) not started — step 3, with the Python `Recipe`.
      - [x] Field table: "What the recipe models" below — every field of the 18 modelled blocks is
            carried; unknown blocks are rejected unless allowed (`md` allows GAMD/EDS, Python nothing).
      - [x] `gromos-run/src/plan.rs`: `AlgorithmSpec` (14 kinds; one spec = one algorithm; fully
            resolved — DOF, virial, `four_pi_eps_i`, `atoms_per_solvent`, restraint paths, `parallel`),
            `TermSpec` (`schnet`, `coupling: delta | replace`), `build_plan`, `validate_plan` driven by
            per-kind `KindRules` (G9), registries `KINDS` / `name()` / `examples()` with a completeness
            test (G4). `build.rs`: `instantiate(plan, topo, conf, periodicity)` — **no recipe parameter**
            (G8, compile-time); `build_sequence_from_recipe`; `build_sequence_from_imd` is a one-line
            front-end. `md @dump` prints recipe + plan (JSON) and exits; `md` writes
            `<tre>.recipe.toml` with the diagnostics as comments; `Simulation.recipe_toml` /
            `.plan_json` / `.diagnostics` expose the same (`py-gromos/tests/test_recipe.py`).
      - [~] `Barostat` + a provider without a virial is rejected by `validate_plan`; `coupling=replace`
            is rejected with the 2.8 message. ✗ **The per-term energy slot on `Contribution` (G10) is
            not done** — it changes `gromos-forces` and stays the precondition for a second `TermSpec`
            variant (step 5 does it first).
      - [x] Tests: `gromos-io/tests/imd_roundtrip.rs` (read → write → read exact on all 41 inputs,
            write idempotent, absent-block defaults, garbage → error, passthrough verbatim, DOFSET/NRE
            synthesis); `gromos-run/tests/recipe_roundtrip.rs` (imd → recipe → imd → recipe, text,
            TOML, JSON, derived defaults, diagnostics, `UnknownBlock`, factories, unknown-field error);
            `plan.rs` unit tests (registry complete, JSON round trip, typo → error, canonical plans
            validate, five violations rejected by name). **Rust-side A/B:**
            `build::tests::plan_matches_legacy_builder_bit_for_bit` — the frozen step-1 builder
            (`legacy_builder.rs`, `cfg(test)`) vs the recipe path, energies/positions/forces compared
            with `to_bits()` after up to 10 steps on all 41 references. **A7:**
            `scripts/roundtrip_imd_gromosXX.py` — **40/40** rewritten inputs accepted by the real
            gromosXX and reproducing a fresh run of the original to the digit (31 also match the
            committed `expected/` exactly; the other 9 — vacuum, COM-removal, NPT, small boxes — differ
            from `expected/` only in ~1e-32 COM-kinetic terms that a fresh run of the *original* input
            shows too: environment noise, invisible at 1e-8). `water_216_nve_nobath` is a regular
            passing reference again; its `ignore`/xfail are gone.
      **Gate met 2026-08-29:** Rust reference 38 passed / 2 ignored; Python green (211 + `test_recipe.py`,
      16 skipped, 15 xfailed); G3/G4/G7/G8/G9 green; A7 40/40. ✗ Not met: `ImdParameters` is still the
      argument of `build_simulation` in `pyo3-gromos` — it becomes `Recipe` in step 3.

**What the recipe models** (every `ImdParameters` field; "kept" = parsed only since step 2)

| Block | Fields | Recipe group | Engine use |
|---|---|---|---|
| TITLE | multi-line title (kept as lines) | `title` | echoed |
| SYSTEM | NPM, NSM | `system` | NSM is a hint; the coordinate file decides |
| INITIALISE | NTIVEL NTISHK NTINHT NTINHB / NTISHI NTIRTC NTICOM / NTISTI / IG TEMPI (NTIRTC, NTISTI kept) | `control` | NTIVEL, NTISHK, NTICOM, IG, TEMPI |
| STEP | NSTLIM T DT | `control` | all |
| BOUNDCOND | NTB NDFMIN | `boundary` | all |
| MULTIBATH | algorithm, NBATHS, TEMP0/TAU per bath, DOFSET, LAST/COM-BATH/IR-BATH (kept) | `ensemble.thermostat` (`None` when absent) | one bath; NBATHS > 1 is a named error, not a truncation |
| PRESSURESCALE | COUPLE SCALE COMP TAUP VIRIAL, SEMI (kept), 3×3 pressure | `ensemble.barostat` | `[0][0]` of the tensors, VIRIAL |
| COMTRANSROT | NSCM | `control.com_motion` | absent = 0 = off |
| PRINTOUT | NTPR, NTPP (kept) | `outputs` | NTPR |
| WRITETRAJ | NTWX NTWSE NTWV NTWF NTWE NTWG NTWB (NTWSE/NTWG/NTWB kept; NTWV/NTWF now frequencies) | `outputs` | NTWX, NTWE, NTWF ≠ 0 |
| CONSTRAINT | NTC, NTCP, NTCP0, NTCS, NTCS0 (kept for SHAKE) | `constraints` | all |
| FORCE | NTF×6, NEGR, NRE | `forcefield.{bonds,…}`, `energy_groups` | NTF; energy groups carried, not used |
| PAIRLIST | ALGORITHM NSNB RCUTP RCUTL SIZE TYPE | `forcefield.pairlist` | all |
| NONBONDED | NLRELE APPAK RCRF EPSRF NSLFEXCL + lines 3–7 (kept) | `forcefield.electrostatics` (+ `lattice_sum`) | APPAK, EPSRF, RCUTL; lattice-sum lines inert |
| POSITIONRES | NTPOR NTPORB NTPORS CPOR | `forcefield.restraints.position` | all but NTPORS |
| DISTANCERES | NTDIR NTDIRA CDIR DIR0 TAUDIR FORCESCALE VDIR NTWDIR | `forcefield.restraints.distance` | NTDIR, CDIR, DIR0, TAUDIR |
| PERTURBATION | NTG NRDGL RLAM DLAMT ALPHLJ ALPHC NLAM NSCALE | `perturbation` | NTG, RLAM, ALPHLJ, ALPHC, NLAM (DLAMT parsed, never applied — open) |
| ENERGYMIN | NTEM NCYC (kept) DELE DX0 DXM NMIN FLIM | `minimisation` | all but NCYC |
| anything else | GAMD, EDS, QMMM, … | `passthrough` | **rejected** unless the caller allows it (`md`: GAMD/EDS, applied out-of-band; Python: none) |
| — | `force_groups` (never populated), `pme_order`/`pme_alpha` (never in a file) | carried | none |

- [x] **Step 3 — Python `Recipe`, `Term`, `Algorithm`; migrate the kwargs (M).** ✓ 2026-08-29.
      - [x] `crates/pyo3-gromos/src/recipe.rs`: `Recipe` (wraps `RunRecipe`; immutable — `update(**groups)`
            deep-merges and returns a copy, `with_term/without_terms/with_inputs/with_execution`; factories
            forwarded from the library; `from_imd/to_imd/save_imd/to_bundle/from_bundle`; `from_dict/
            to_dict/to_toml/to_json` via `pythonize`; `==`; pickling), `Term(kind, **params)` and
            `Algorithm(kind, **params)` validated by serde against the tagged enums — unknown kind or
            field → `gromos.exceptions.RecipeError` (not a dedicated `UnknownParameter`: one class per
            failure kind, four in all, each subclassing the builtin the binding raised before);
            `gromos.terms()` → `{"kind", "params" (example), "feature", "available"}`,
            `gromos.algorithms()` → `{"kind", "params", "rules"}`, `gromos.build_info()`.
      - [x] `Simulation(system, recipe)`; `Simulation(system, recipe, plan=plan)` (re-validated);
            `recipe.plan(system)` returns the plan as a `Plan` (by index and by kind: `insert/insert_after/
            insert_before/remove/replace/validate`, JSON) — a `Plan` class instead of
            `AlgorithmSequence.from_recipe`, so the descriptor type can go in step 4; `Simulation.plan`,
            `.recipe`, `.diagnostics`. The twelve descriptor pyclasses were *not* rewritten as `Algorithm`
            constructors: they die with the descriptor path in step 4.
      - [~] `coupling: "replace"` rejected with the 2.8 message (A8, since step 2). `units` and provenance
            (model checksum in the bundle): ✗ not done — folded into step 5 with the first engine term.
      - [x] Shims (A9) — in Rust (`parameters::deprecated`: one `DeprecationWarning` naming the
            replacement), not in a `_deprecation.py`: `InputParameters` → `Recipe.from_imd`/factories; the
            six `Simulation` kwargs → `recipe.with_inputs(...)` / `with_term(Term("schnet", …))`;
            `AlgorithmSequence.nve/nvt/npt/minimize/from_parameters` + `Simulation.from_sequence` →
            `recipe.plan(system)` + `plan=`. Every shim *translates into a recipe* — no second path
            (`test_recipe.py::test_deprecated_forms_warn_and_are_translations`, parity A vs B below).
            Migration table: `py-gromos/docs/user-guide/quick-start.md`; `__version__` from package
            metadata.
      - [x] `gromos.pyi` `stubtest`-clean (`MYPYPATH=python python -m mypy.stubtest gromos.gromos
            --allowlist stubtest_allowlist_no_ml.txt`; no allowlist on an `ml` build): pyo3 classes are
            `@final` and construct through `__new__`, slot methods are positional-only, `_AlgorithmKind`/
            `_TermKind` literals; docs and notebooks `01`/`02` on `Recipe` (`examples/06` drives the
            binaries and never used `InputParameters` — untouched); G5 tests `test_pyi_lists_every_kind`
            and `test_every_kind_has_a_parity_case` (`orchestrator` exempt until step 5, by name).
      - [x] Non-`ml` build: `Term("schnet", …)` constructs with `available=False`; `Simulation` raises
            `MissingFeatureError` (A14) — `test_terms_registry_and_the_schnet_term`.
      **Gate met 2026-08-29:** parity A (`InputParameters` + kwargs) vs B (`Recipe`) vs D (`Recipe` + its
      own `Plan`) exact (`array_equal`) on all 40 systems in `REFERENCE_SYSTEMS`, no xfail. `Serial`
      was not needed: all three go through one builder with `ParallelPolicy::Auto`, i.e. the same kernels
      at the same thread count (run-to-run determinism at fixed thread count: BENCHMARKING.md). `maturin
      develop` with and without `--features ml` (ml: build, full suite 311 passed / 13 skipped / 15
      xfailed, strict stubtest). Python suite 308 passed / 16 skipped / 15 xfailed (was 215). Rust:
      `gromos-run`/`gromos-io` tests and the 38/40 reference suite unchanged; `pyo3-gromos` gains three
      unit tests. Found on the way: `System.from_files` rejected a solute
      topology with a solvated coordinate file (the reference suite never hit it — it passed `Topology`
      and `Configuration` separately); fixed (`system::atom_count_ok`, unit-tested).

- [x] **Step 4 — Delete the copies (S).** ✓ 2026-08-29. Removed: `resolve_algorithm_sequence`,
      `AlgorithmDescriptor`, the twelve descriptor pyclasses and the preset bodies (`algorithm_sequence.rs`
      1473 → 101 lines: only the deprecated `AlgorithmSequence.nve/nvt/npt/minimize/from_parameters` names,
      each returning the `Plan` of the parameters via `gromos_run::build_plan`), `RestraintFiles`/
      `MlPotentialSpec` (gone since step 3), `build_simulation_from_sequence` (the `AlgorithmSequence::new()`
      in `simulation.rs`), the ad-hoc `PyValueError` sites for recipe/atom-count errors (→ `RecipeError`),
      `test_basic.py`/`test_advanced_features.py`, `md_runners.py`; `analysis` out of the default namespace;
      `Vec3`/`Frame`/`rmsd`/`rdf` to `f64` (`gromos_core::math::Vec3`); the `Simulation.add_ml_potential`
      mention; `py-gromos`'s `glam` dependency. Added: `gromos_run::read_imd` — the `md` binary,
      `Recipe.from_imd`, `InputParameters` and the bundle loader all open parameter files through it — and
      the G6 gates in `just lint` (four greps after clippy; `just` is not installed on this machine, the
      recipe body was run directly). `Simulation.recipe`/`.plan` are always present.
      ✗ Narrowed: the gate greps `read_imd_file(`/`parse_imd_str(`, not `parse_file(` — `GamdBlock`/
      `EdsBlock::parse_file` in `md.rs` re-read the file for the passthrough blocks they own (parsing them
      from `recipe.passthrough` is 1.9's job, GaMD/EDS as algorithms).
      **Gate met 2026-08-29:** G6 clean; `algorithm_sequence.rs` 101 lines; Rust `gromos-run`/`gromos-io`/
      `pyo3-gromos` tests and the 38/40 reference suite green; Python suite 300 passed / 8 skipped / 3 xfailed
      (was 308 / 16 / 15: the eleven placeholder skips and the twelve descriptor xfails are gone); stubtest clean;
      CHANGELOG 0.0.31; the four CONTEXT.md files name the recipe as the only entry point.
      `test_front_end_parity.py` lost path C and its twelve `xfail` rows; the deprecated names are now
      checked as a translation (`from_sequence` shim vs `Recipe`, exact) on every non-restraint system.

- [x] **Step 5 — First new term through the new door (S), as the proof.** ✓ 2026-08-29.
      `Term("xtb", region=…, elements=[…], gfn=2, charge=0, multiplicity=1, work_dir=None, timeout_s=600,
      coupling="delta")` over `XtbInteraction`, no feature gate.
      (i) wiring, measured: **three files in `gromos-run`** — the `TermSpec::Xtb` variant (`recipe.rs`),
      five registry match arms (`plan.rs`: `KINDS`/`name`/`examples`/`feature`/`provides_virial`/
      `coupling`) and one `instantiate` arm (`build.rs`) — plus the `_TermKind` literal in `gromos.pyi`
      and the test. One more than the target of two: the registry is a separate file from the
      variant by design (G5's exhaustive matches), not a leak — nothing else knew about the term
      (no binding change, no md.rs change, no `.pyi` class). The orchestrator is no longer behind
      `--features ml`; only the SchNet arm is.
      (ii) physics — `py-gromos/tests/test_xtb_term.py` on the water-dimer fixture (solute water
      carrying the term, unconstrained, classical bonded terms off; rigid SPC solvent water): adding
      the term adds exactly `XtbPotential`'s direct energy (`rel 1e-9`) and forces of the region and
      nothing outside it; two terms report separately and add up; NVE through the real step loop
      within the Rust oracle's thresholds (fluctuation < 0.5 %, half drift < 0.2 %; measured
      4 × 10⁻⁵ % and 1 × 10⁻⁶ % over 120 steps at dt = 0.1 fs, T ≈ 297 K); `coupling="replace"` rejected with the 2.8 message (`PlanError`); `Barostat` + xtb
      rejected (no virial); an `elements` table not covering the region → `RecipeError`.
      ▲ `work_dir` per term index (`<tmp>/gromos-rs-xtb-term-<i>` unless given) and a subprocess
      timeout (`XtbInteraction::with_timeout`; `qm_subprocess::run_subprocess_with_timeout` kills the
      child) — `timeout_s=0` is a `RunError`, not a hang.
      **G10, in memory:** `Energy.term_energies` (registry name → energy, `xtb:0`/`xtb:1` for a
      repeated kind) filled by `ProviderOrchestrator::evaluate_with_terms`; `Simulation.term_energies`.
      ✗ not in the `.tre` file nor in `run()`'s columns yet — the `.tre` format and the 12-column
      array are gromosXX-shaped; a per-term block/column set is a separate decision (FUTURE.md).
      ✗ provenance (model checksum in the bundle) still not done; `units` on terms not modelled.

**Definition of done for 3.9**
- [x] `AlgorithmSequence::new()` appears in exactly one non-test file, `crates/gromos-run/src/build.rs`,
      and `instantiate` has no `&RunRecipe` parameter. ✓ step 4 (G6 gate in `just lint`).
- [ ] `Simulation`'s constructor has no per-feature kwargs; `gromos-run` has no `process::exit`.
      `process::exit` ✓ (G6). The `distrest`/`posresspec`/`refpos`/`ml_*` kwargs survive as deprecation
      shims (warn, translate into the recipe) until the next release — remove them then.
- [~] Every reference system passes through all three front-ends bit-identically (G2, `Serial`), every
      `to_imd` output and every factory output is accepted by gromosXX (G3/G7/A7), and every plan an
      `insert_after`/`remove` can produce is either valid or rejected with a named reason (G9).
      G2 ✓ (steps 3–4, with `Auto` on every path rather than `Serial`); A7 ✓ for the 40 reference
      `.imd` files (step 2); factory outputs through gromosXX: not re-verified since step 2; G9 ✓
      (`validate_plan` rules, step 2; every `Plan` edit is re-validated, step 3).
- [x] Adding a term = one variant + one arm (step 5 measured, number recorded here) **and** a physics
      oracle passed. ✓ step 5: variant + registry arms + `instantiate` arm (three `gromos-run` files),
      `test_xtb_term.py` as the oracle.
- [x] `water_216_nve_nobath` passes; the divergence table has no open rows. ✓ (step 2; the A-vs-C rows
      closed by deletion in step 4).

**Divergence table** (rows open as they are found; close as steps land — keep the closed ones)

| System | Paths | First differing quantity / step | Cause | Resolution | Closed in |
|---|---|---|---|---|---|
| `water_216_nve_nobath` (any `.imd` without MULTIBATH) | all gromos-rs paths vs gromosXX | E_total at frame 0: −4270.69 vs −4270.19 (Berendsen scaling already acts in step 0) | parser kept `TempBathParameters::default()` (Berendsen 300 K, τ 0.1) — verified `imd.rs:306-317,329,422-500`; gromosXX: "Adding a bath, no temperature coupling" | defaults now mean "absent" (A18); the reference passes in both suites | **closed, step 2** |
| any `.imd` without COMTRANSROT | all gromos-rs paths vs gromosXX | COM velocity removed at step 0 | `nscm` defaulted to 100000 (removal at step 0); gromosXX `parameter.h`: `skip_step 0` | default `nscm = 0`; `imd_roundtrip.rs::absent_optional_blocks_mean_gromosxx_defaults` | **closed, step 2** |
| any `.imd` with > 1 bath | all gromos-rs paths vs gromosXX | thermostat scaling, step 1 | every builder read `temp_bath[0]` (`simulation.rs:483-492`, `md.rs:357-367`) | `ensemble.thermostat` carries every bath; `build_plan` rejects NBATHS > 1 with a named error instead of truncating (one-bath thermostat) | **closed (rejected loudly), step 2** — multi-bath support itself is a physics item |
| `nacl_1water_settle` | A vs C | kinetic, frame 0 (1.68259628e-3 vs 1.68252320e-3); positions ~7e-8 | C runs SHAKE where NTCS=3 asks SETTLE (only `ShakeConstraints` exists) | one builder (step 1); descriptor path deleted (step 4) | open → step 1 |
| `nacl_1water_lincs` | A vs C | kinetic, frame 0 | C runs SHAKE where NTCS=2 asks LINCS (solvent) | as above | open → step 1 |
| `aladip_vacuum_lincs` | A vs C | kinetic, frame 0 (32.6256 vs 32.6709) | C runs SHAKE where NTCP=2 asks LINCS (solute) | as above | open → step 1 |
| `water_216_nvt_nosehoover` | A vs C | kinetic, frame 0 (2606.29 vs 2605.82) | C runs Berendsen where MULTIBATH algorithm = 1 asks Nosé-Hoover | as above | open → step 1 |
| `water_216_nvt_nhc_chain` | A vs C | kinetic, frame 0 (2606.29 vs 2605.82) | C runs Berendsen where MULTIBATH algorithm ≥ 2 asks an NH chain | as above | open → step 1 |
| `water_216_box_com_rot` | A vs C | kinetic, frame 0 (2603.98 vs 2605.10) | `RemoveCOMMotion(initial: bool, nscm)` collapses NTICOM=2 to 1 — rotation removal lost (not in the 3.8 audit) | as above | open → step 1 |
| `aladip_trunc_oct` | A vs C | kinetic, frame 0 (165.756 vs 165.807); positions diverge | `build_simulation_from_sequence` hard-codes `SimBox::rectangular`; no NTB=−1 transform | `prepare_system` (step 1) | open → step 1 |
| `nacl_1water_distres` | A vs C | not runnable on C | restraints are `Simulation` kwargs; `from_sequence` has none | `RunInputs` (step 1) → `forcefield.restraints` (step 2) | open → step 1 |
| `aladip_solvated_em_posres` | A vs C | not runnable on C | same (posresspec/refpos) | as above | open → step 1 |
| `aladip_solvated_em` | A vs C | not runnable on C | same (posresspec/refpos) | as above | open → step 1 |
| `ch4_water_fep`, `aladip_vacuum_fep` | A vs C | kinetic, frame 0 | C never resolves the PERTURBATION block onto `Forcefield` (λ, soft-core, NLAM) | recipe `perturbation` group (step 2) | open → step 2 |
| six >100-atom systems (`water_216_*`, `water_1000_spc_gridcell`) | A vs C, after step 1 | last-digit energies from frame 0 | A inherited the binary's `ParallelPolicy::Auto`; C's resolver built a serial `Forcefield` (A2) | resolver applies the same policy + `four_pi_eps_i` | **closed, step 1** |
| every reference system | A vs B vs D (`InputParameters` vs `Recipe` vs `Recipe` + `Plan`) | none — `array_equal` on energies, positions, forces | — | the shims translate into the recipe; `plan=` is stage 1 of the same builder | **closed, step 3** |
| (rows above marked "open → step 1/2") | | | | path C — the descriptor resolver — deleted with its xfail table; the names that survive are translations into the recipe, checked exact against path B | **closed, step 4** (deleted, not fixed) |

**Explicitly out of scope here** (tracked elsewhere): the `System` building algebra (FUTURE.md,
D1–D8); zone-aware `Forcefield` and `coupling: Replace` (2.8); GaMD/EDS/REMD as algorithms (1.9); PME
(deferred breadth); `Simulation.set_parameter(λ, T0)` runtime knobs and `set_positions()` context
reuse (OpenMM's `Context.setParameter` lesson — a follow-up once `start` exists); `dlamt` λ scheduling
(parsed, never applied); per-step charge refresh for path-(c) engines (2.8).

**Superseded:** the July 2026 "composition-pattern refactor" audit (`PLAN_ARCHIVE.md` §Deferred,
`~/.claude/plans/golden-baking-liskov.md`) — its trigger has arrived (the ML/QM/restraint terms are the
second caller it waited for) and its A/B condition is step 0.

**Deferred to FUTURE.md (do NOT start until 3.1–3.4 are solid and usable)**
> The system-building algebra and everything that requires the open D1–D8 decisions:
> `ForceField`/`molecule(seq)` (subprocess `make_top` + coordinate pipeline), the `+`/`*`
> assembly algebra, `.solvate()`/`.pack()`, `.neutralize()`, native topology building, and the
> `system_builder.py` completions. These live in FUTURE.md and the design mockup; revisit with
> running code in hand. The traditional 8 binaries stay available throughout for back-compat.

### Priority 4 — Benchmarking ✓ Phases 0–3 done (2026-08-29); see `BENCHMARKING.md`
Engine-vs-engine against gromosXX 1.6.0 (both `-march=native`), protocol + assumptions +
every result in `BENCHMARKING.md`; harness `scripts/bench_engines.py` (n, mean ± std, interleaved
repeats), `scripts/make_solvated_box.py` (24k-atom production box via `sim_box`), `just bench`.
- [x] Single core: 0.48×/0.77× → **0.97×/1.11×** on the 648/3 000-atom references, **1.14×** on
      the 24k production box (sixteen changes, all reference-suite-verified; per-phase breakdowns).
- [x] Threads: pairlist, long-range and CG grouping were serial → now parallel; 8 cores **1.38×**
      vs gromosXX-OpenMP on the 24k box (from 0.53×), S(8) = 4.1×; SMT point measured (no gain).
- [x] O(N) cell list confirmed and fixed (was pruning nothing: grid from IMD `SIZE`, distance-pruned).
- [ ] 81 000-atom box (`make_solvated_box.py --replicate 3`); remaining serial phases: SHAKE/SETTLE
      over molecules, integration, thermostat.
- [ ] Phase 4 MPI — blocked on four build-plumbing items (BENCHMARKING.md §6.2); single node can
      only prove correctness, not speed. Phase X (integrated GPU) recorded as a future spin-off.
- [ ] Criterion micro-benches (`cargo bench`) remain the per-kernel regression guard; not yet baselined.

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

**Python (build the extension first; add `--features ml` for the ML term):**
```sh
cd py-gromos
uv run maturin develop --release
uv run pytest tests/ -v
uv run pytest tests/test_front_end_parity.py -v      # 3.9 G2 — exists after step 0
```
