# PLAN_ARCHIVE.md — full working notes of finished PLAN.md items

Moved out of `PLAN.md` on 2026-08-29 to keep the roadmap readable. Each section below is the
*verbatim* text PLAN.md carried when the item was closed — the traps, oracle numbers, rejected
alternatives and "why" that a one-paragraph summary cannot hold. Nothing here is a task: PLAN.md keeps
the checkable state, CHANGELOG.md keeps the per-version record, this file keeps the reasoning.

Sections: §9a (Dim 9 plumbing and CellList activation), §2.6–§2.9 (provider seam, electrostatic
embedding, QM/MM wiring, ML pipeline), §3.5–§3.7 (py-gromos convenience path, parity gap, ML binding),
§Deferred (the July 2026 composition-pattern audit, superseded by PLAN.md 3.9), §2.10 (energy
library profile).

Line references inside the archived text (`md.rs:98`, `simulation.rs:66-295`, …) are as of the date the
item closed and may have moved since.

---

## §9a — Dim 9 plumbing and CellList activation (PLAN.md 1.5, archived 2026-08-29)

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

---

## §2.6–§2.9 — provider seam, embedding, QM/MM wiring, ML pipeline (PLAN.md Priority 2, archived 2026-08-29)

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

---

## §2.10 — Energy library profile: same grammar, stricter rules, versioned (PLAN.md Priority 2, archived 2026-09-03)

Spec written 2026-09-03, all four steps landed in 0.0.47 the same day. Gate results are appended to each step.

> The `ene_ana` library is the right design for a research code — a researcher can re-partition a
> system, recompute with a new topology and describe the resulting file without recompiling — and
> it has one hole: nothing ties a library to the trajectory it is applied to. `ENER 26 2` instead of
> `ENER 52 1` reads every value into the wrong slot and prints `totkin = −52.8` with no error;
> `ENEVERSION`, the only guard, is hand-typed in both files (`out_configuration.cc:58`), only a
> *warning* in gromos++ (`ene_ana.cc:463`), and the tutorial's stale library carries the same string
> as md++'s current one. The fix is not a new format — that would forfeit the oracle and the
> tutorials — but a **profile** of the existing one: every profile addition is a column-0 `#` line
> (both upstream readers strip them; a new block would break gromos++'s strict block order,
> `EnergyTraj.cc:454`), the layout is defined in code and everything else is derived from it, and a
> library that disagrees with the file's fingerprint is refused with a diff.
>
> **The spec is `docs/src/reference/energy-library.md`** (profile v1): grammar, canonical form and
> SHA-256 fingerprint (names excluded, so renaming stays free), trajectory and library
> self-description lines, the three reader tiers, the structural checks that hold for any
> partitioning, the official library from `gromos-io`, provenance, versioning, refusal catalogue.
> Steps below implement it; the spec's Status table is updated as each lands.

Options weighed for "who owns the official library" (recorded in the spec §5): ship the file and
locate it at runtime (gromos++'s model — editable, unverified, the status quo); embed the file text
(immutable, but writer and text can still disagree); **layout in code with library, fingerprints and
a writer round-trip test derived from it** (chosen); signed distribution (wrong threat model — the
danger is accident and drift, and a fingerprint line already makes any edit visible).

- [x] **Step 1 — Layout in code, fingerprint, official library (M).** `EnergyLayout` in
      `gromos_io::energy_traj`: the 2023-04-15 `ENERTRJ`/`FRENERTRJ` declarations as data;
      `canonical_form()` and `fingerprint()` per spec R2; `official_library()` generates the schema
      sections and appends md++'s `VARIABLES` (embedded text, GPL-2 like ours). Checked-in copy at
      `crates/gromos-io/data/ene_ana.md++.lib`, a test asserting it equals the generated one and
      that its fingerprint equals md++'s `data/ene_ana.md++.lib`'s (so we prove we ship *their*
      layout). `size` slots validated before allocation (R6, first two bullets — closes the
      `f64 as usize` hazard). Mutation tests on the library text: reshape (`ENER 26 2`), shift
      (`ENER 53`), dropped `size`, renamed subblock and size (must give the *same* fingerprint).
      **Gate:** fingerprint of the checked-in library == fingerprint of md++'s; rename passes,
      every other mutation changes the fingerprint; the ten gromos++ reference cases unchanged.
      **Result:** landed as `Schema` (not `EnergyLayout`) — the same type parses a library's
      section and holds the built-in layout, so there is one canonical form. ENERTRJ
      `sha256:aeede256…fbda96`, FRENERTRJ `sha256:2cacf7dd…3c3f`; equal to md++'s library by
      test and to `sha256sum` of the canonical text by hand. Rename → same fingerprint;
      reshape/shift/dropped size/the tutorial's stale library → different. Official library
      byte-equal to the generator (`UPDATE_OFFICIAL_LIBRARY=1` rewrites); gromos++ `ene_ana`
      gives identical output with it and with md++'s own. Reference cases unchanged.
- [x] **Step 2 — Reader tiers and refusals (M).** `EnergyTraj::open` establishes the tier (R5)
      and applies R6 on every frame, with the exact first lines of spec §3. `ene_ana @library`
      becomes optional (official library when absent); a `VARIABLES`-only library is accepted.
      `ext_ti_ana` reads `.trg` through `energy_traj` instead of its own parser, so it gets the
      checks for free.
      **Gate:** the tutorial's stale `ene_ana.md++.lib` against a gromosXX `.tre` fails at open
      with the tier-**b** diff (`SPECIAL: library … x 12, built-in … x 13`) instead of
      `BONDED[1][1] outside`; `ENER 26 2` fails on frame 1 by the invariant with no fingerprint and
      at open with one; the reference cases still byte-identical to gromos++.
      **Result:** `read_preamble` + `EnergyTraj::bind` (the reader is a `TextReader` the caller
      owns, so there is no `open`). Stale tutorial library against gromosXX's `.tre`: refused at
      open, `ENERGY03 ENER:   library 50 x 1,   built-in 52 x 1` / `SPECIAL: library
      NUM_ENERGY_GROUPS x 12, built-in NUM_ENERGY_GROUPS x 13`; against our self-described
      `.tre` the tier-**a** diff with both provenance lines. `ENER 26 2` without a fingerprint:
      first-frame identity warning naming the tier; with one: refused at open. The identity is
      `ENER[1] = ENER[2] + ENER[3] + ENER[21]` (`energy.cc:599`, special total at slot 21 —
      the spec said `[2]+[3]`; corrected) and is a **warning**, not a refusal: a verified layout
      is a fact about slots, not about the writer's sums, and our own pre-0.0.47 files fail it
      (see Step 3). The `ext_ti_ana` item is dropped: ours is a trapezoid over `⟨dH/dλ⟩`
      through the fixed-slot `read_free_energy_trajectory`, not a port of gromos++'s
      library-driven program; it should go through `energy_traj` if it ever becomes that port.
- [x] **Step 3 — Writers describe what they wrote (S).** `EnergyWriter` and `FreeEnergyWriter`
      emit the R3 lines from `EnergyLayout` and the R8 provenance (writer version, topology hash
      and path, energy-group ranges); a round-trip test writes a frame and reads it back under
      tier **a**, so the hand-written writer is checked against the layout it claims.
      **Gate:** gromos++ `ene_ana` (`.local/gromosPlsPls/gromos++/BUILD/programs`) on our `.tre`
      gives byte-identical output before and after; `md`'s `.omd` header prints the fingerprint.
      **Result:** both met. Pinning the first-frame identity to gromosXX's slot numbers found
      what L4 warned of, the other way round: the writer filled slots the layout has, but the
      wrong ones — special at `ENER[19]`, SASA `[20]`, constraints `[22]`, distance restraints
      `[23]` instead of `[21]/[22]/[24]/[25]` (`out_configuration.cc:3074-3125`); `dhdl_special`
      likewise in `.trg`. Every `ene_ana`
      `totspecial`/`totdisres`/`totconstraint` on a gromos-rs file had read 0. Fixed;
      `EnergyFrame::special` is now gromosXX's `special_total` (restraints included).
- [x] **Step 4 — `ene_ana @print_library [trajectory]` (S).** Prints the official library, or the
      library a self-describing trajectory says it was written with (spec §4). The tutorial-1 notes
      and `crates/gromos-analysis/tests/data/README.md` switch from the checkout's library to it.
      **Gate:** `@print_library md.tre | ene_ana @library /dev/stdin …` reproduces the default run.
      **Result:** met; `@print_library` alone equals the checked-in file; editing a `subblock`
      line in the printed library makes it refuse (R4). The reference test and README use the
      built-in library; the `.local` md++ copy is no longer needed to run the suite.
- [ ] Later, Priority 3: expose `EnergyTraj` to py-gromos as a typed table (frame × term × group
      pair) — the re-partition workflow is a `group_by` there, which is where a richer derived
      language belongs rather than in the `.lib` (spec §5, last paragraph).

**Assumptions** — all five held; L1 by the Step 3 gate, L2 by the rename test, L3 by the tier-**b** gate, L4 by the round trip, L5 by `cargo tree` (`sha2` is new to the workspace and brings seven small RustCrypto crates: `digest`, `block-buffer`, `crypto-common`, `generic-array`, `typenum`, `cpufeatures`, `cfg-if`).

| # | Assumption | If wrong | Check |
|---|---|---|---|
| L1 | Column-0 `#` lines are dropped by both upstream readers wherever our files reach them (`Ginstream.cc:165-172`, `blockinput.cc:53-58`). | Our `.tre` breaks gromos++ → R3 lines move to the `TITLE` block text. | Step 3 gate runs the gromos++ build on our output. |
| L2 | Subblock and size *names* carry no information about where a byte lands; only shapes, order and block names do. | A name that changes reading → it enters the canonical form (profile v2). | Step 1 rename-mutation test: same fingerprint, same values. |
| L3 | md++'s `ENEVERSION` string does **not** identify a layout uniquely in the wild (the tutorial proves it). | — this is why tier **b** compares fingerprints, not strings, and why R6 exists for files with no self-description. | Tier-**b** gate in step 2 uses the stale tutorial library. |
| L4 | The hand-written `EnergyWriter` sequence equals the 2023-04-15 layout; it is checked by round trip, not made layout-driven (a value→slot map would be more machinery than it saves). | A slot the writer fills and the layout lacks → the round-trip test fails, not the user. | Step 3 round-trip test under tier **a**. |
| L5 | `sha2` (pure Rust, default features) is an acceptable `gromos-io` dependency. | FNV-64 in-tree, spec R2 amended to say so. | `cargo tree -p gromos-io`. |

---

## §3.5–§3.7 — py-gromos convenience path, parity gap, ML binding (PLAN.md Priority 3, archived 2026-08-29)

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

---

## §Deferred — the July 2026 composition-pattern audit (superseded by PLAN.md 3.9)

**Deferred — composition-pattern refactor** (superseded by 3.9 — the trigger it waited for has
arrived: the ML, QM and restraint terms are the "second caller")
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

---

## §3.9-review — five framework-perspective reviews of the composition model (2026-08-29)

Five independent reviews of PLAN.md 3.8/3.9 as it stood before the consensus was folded in (that pre-review text is what the reviews cite by line number; the post-review 3.9 in PLAN.md differs). Each reviewer was given the same non-negotiable constraints (GROMOS-faithfulness at 1e-8, the `.imd` as first-class input, Rust owns physics, explicit step order) and the same output format. The consensus table and the adopted changes C1–C13 are in PLAN.md 3.9; the code claims below were verified before adoption.

### Claims verified in the code

- **Absent MULTIBATH ⇒ silent Berendsen thermostat.** `read_imd_file` starts from
  `ImdParameters::default()` (`crates/gromos-io/src/imd.rs:329`), whose `temp_bath` is
  `vec![TempBathParameters::default()]` = Berendsen, 300 K, τ = 0.1 (imd.rs:306-317). `temp_bath` is
  only replaced inside the `"MULTIBATH" =>` arm (imd.rs:422-500). Both `md.rs:362-367` and
  `simulation.rs:488-493` then compute `thermostat_on = tau > 0.0` → true. gromosXX runs no bath
  when the block is absent. Only two reference inputs lack MULTIBATH (`aladip_solvated_em`,
  `aladip_solvated_em_shake`) and both are EM, which skips the thermostat — the suite cannot detect it.
- **Silent numeric coercion.** `parse_f64`/`parse_i32`/`parse_usize` (imd.rs:1072-1082) map any
  unparseable token to `0.0`/`0`.
- **Restraints are `Forcefield` fields, not providers** (`forcefield.rs:124-130`; set in
  `simulation.rs:329-410` and `md.rs:1140-1147`); only LjCrf/SchNet/Xtb/Mopac/Embedding implement
  `PotentialProvider`. `ProviderOrchestratorAlgorithm` overwrites `special_total`
  (`orchestrator_algorithm.rs`), so a second orchestrator or coexistence with position restraints is
  unsupported today.
- **GAMD/EDS are parsed out-of-band** (`md.rs:646-668`, `GamdBlock::parse_file` / `EdsBlock::parse_file`
  re-read the input file) and applied after `run_step` (`md.rs:1591-1700`) — a second IMD reader.
- **`dlamt` is parsed (imd.rs:170) but never applied** by `md.rs`.
- **`ntc_from_constraints` + the `nve/nvt/npt/steepest_descent` factories on `ImdParameters`**
  (imd.rs:943-1053) are a defaults site the 3.8 audit did not list.
- **Unit constants of different vintages:** `HARTREE_TO_KJMOL = 2625.4996` (units.rs:96, 8 s.f.) next to
  `BOHR_TO_NM = 0.0529177210903` (units.rs:101, 12 s.f.); no reference system exercises them.
- **`XtbInteraction` virial is `Mat3::ZERO`** (`xtb.rs:359`, "not computed yet") — an NPT recipe with an xtb
  term would mispressure silently.
- **Orchestrator energy is one scalar** (`orchestrator.rs:96`, `energy += contribution.energy`) — no per-term
  channel; parity tests on totals cannot see compensating errors between terms.
- **xtb subprocess: fixed filenames in `work_dir`** (`region.xyz`, `region.engrad`, `pcharge`, `pcgrad`,
  `xtb.rs:292-295`) and no timeout in `qm_subprocess.rs` — two xtb terms or two threads would collide.
- **SchNet casts positions to f32** (`schnet.rs:201,210`, `Kind::Float`) — bit-reproducibility of an ML
  term across runs is not guaranteed (libtorch thread pool).

### Review 1 — OpenMM perspective

## Verdict (3 sentences)
This is `ForceField.createSystem → System → Context` done correctly for GROMOS: `build_plan` is `createSystem`, `Vec<AlgorithmSpec>` is the inspectable `System`+`Integrator`, `instantiate` is `Context`, and making the plan pure data sidesteps OpenMM's worst wart (`Context.reinitialize()`) instead of inheriting it. Its weakest seam is the one OpenMM never fully closed: at `Simulation(system, recipe, sequence=…)` (PLAN.md:542-543) two objects can describe the same dynamics and nothing yet says which fields win. Write that contract, add plan validation in `instantiate`, and the rest is execution.

## Scorecard
D1 `Modify` — grouping and single-site defaults beat OpenMM's defaults scattered across `Force` constructors (`ewaldErrorTolerance`, `constraintTolerance`), but `parallel: Auto|Serial|Parallel` (l.512) is a `Platform` concern, not physics: it has no IMD field, so `to_imd` cannot be total; pass it to `instantiate`/`Simulation(platform=…)` as OpenMM does.
D2 `Agree` — exactly `createSystem()` (needs Topology, l.527) vs `Context` (needs positions, l.528); an editable `System` *before* `Context` is what made openmmtools/perses possible.
D3 `Modify` — serde is `XmlSerializer` and that aged well only because every `Force` carries a `version` and proxies migrate; add `version` to `RunRecipe` now, and l.651's "unknown parameter → ValueError" needs `#[serde(deny_unknown_fields)]` on every variant — serde silently ignores unknown keys by default.
D4 `Modify` — restraints are not providers: `PositionRestraints`/`DistanceRestraints` are fields of `Forcefield` (forcefield.rs:124-130); only LjCrf/SchNet/Xtb/Mopac/Embedding implement `PotentialProvider`. `TermSpec::PositionRestraints` (l.518) therefore either moves restraint physics into a provider (a physics-path change that must stay bit-identical, and the orchestrator *overwrites* `special_total`, orchestrator_algorithm.rs:23) or gets an `instantiate` arm that mutates `Forcefield`. Say which; step 5's "two files" count depends on it.
D5 `Modify` — OpenMM only ever *read* prmtop/top/psf; writers live in ParmEd because a foreign writer needs model context. `to_imd()` (l.536) takes no `&Topology`, yet a Python-born recipe must synthesise MULTIBATH `LAST` and FORCE `NRE`. Give `to_imd` the topology.
D6 `Agree` — at 0.1.0 one cycle is right; the `simtk.openmm`→`openmm` shim lasted years only because of install base.
D7 `Modify` — the right hatch (it is `CustomIntegrator`, serialisable), but the docstring's own `seq.remove(...)` pattern applied to `Shake` leaves thermostat DOF computed from `imd.ntc` (algorithm_sequence.rs:1155); DOF must derive from the plan's constraint entries at `instantiate`.
D8 `Agree` — strict-xfail then lift-verbatim is the platform-vs-Reference discipline OpenMM runs on every commit; add a same-path-rerun determinism baseline before any cross-path `array_equal` (rayon kernels, l.456).

## What OpenMM learned the hard way that this plan should adopt
1. **Validate at binding, not at step.** `Context()` throws for PME without box vectors; `instantiate` must reject: no `Forcefield`, more than one, `Orchestrator` not immediately after `Forcefield` (it *adds* to forces `Forcefield` overwrote, orchestrator_algorithm.rs:6-11), `EnergyCalculation` not last, `TemperatureCalculation` before constraints. `insert_after` makes every one reachable from Python.
2. **Separate runtime knobs from structure.** `Context.setParameter` (λ, k, T0) vs `updateParametersInContext` vs `reinitialize()`. `Forcefield.lambda` exists (forcefield.rs:130); `dlamt` is parsed (imd.rs:170) and never applied by `md.rs`. Reserve `Simulation.set_parameter(name, value)` routed to an `Algorithm::set_parameter` (Rust validates) so TI/annealing never rebuilds.
3. **Context reuse.** openmmtools' `ContextCache` exists because creation is expensive. `build_simulation` rebuilds the pairlist and LJ matrix and reloads the SchNet model per construction; add `Simulation.set_positions()` that re-runs `init` + step 0 — the training loops in `scripts/` will want it.
4. **Outputs are not the model.** `Reporter`s attach to `Simulation`, never to `System`. Keep `outputs` (l.516) for the IMD round trip only; Python writes files through explicit reporters, never because `recipe.outputs.ntwx` was set (`md.rs:352-353,1381` honours it; `Simulation` cannot — a front-end divergence A15 hides).
5. **Determinism is a platform property.** CUDA needed `DeterministicForces`. Step 0 should prove path (a) == path (a) on rerun before comparing (a)↔(b), or the divergence table fills with noise.

## What NOT to copy from OpenMM
1. `CustomNonbondedForce`/Lepton energy strings and runtime plugins — physics-in-strings from Python violates "Rust owns all physics" and cannot be verified at 1e-8 against gromosXX; the closed exhaustive enum (G4) is correct.
2. The unordered `Force` set + `setForceGroup` execution model — the explicit ordered sequence is the constraint *and* the debuggability win; don't let force groups sneak in as ordering.
3. A mutable `System` after `Context` plus `reinitialize(preserveState=True)` — OpenMM's number-one FAQ. Expose `sim.plan` as a frozen snapshot, not a handle.

## Risks the plan misses
1. Two orchestrators (ML term today, restraint terms tomorrow) break the `special_total` overwrite semantics (orchestrator_algorithm.rs:23); l.523 shows one `Orchestrator {terms}` — enforce exactly one.
2. G3's `==` on `f64` after `to_imd` fails unless the writer emits shortest-round-trip floats (`{:?}`), not `{:.4}`.
3. G7 has two default sites by construction (`ImdParameters::default()` and `RunRecipe::default()`) and only tests agreement; derive one from the other.
4. `refpos` defaults to initial positions (simulation.rs:357-358), so `TermSpec::PositionRestraints` resolves at `instantiate` from `Configuration`; spec it as `refpos: Option<Path>` with `None = initial` or `to_dict` of the same plan differs run-to-run.
5. `passthrough` GAMD/EDS blocks are inert in Python but active in `md.rs` (md.rs:1591-1700); `instantiate` should error when passthrough carries a block the front-end will not honour.

## The one change you would insist on before implementation starts
Write the precedence contract before step 2: `instantiate(plan, recipe, …)` reads **only** `control`, `boundary`, `perturbation` inputs and `passthrough` from the recipe; `ensemble`, `constraints`, `integrator`, `forcefield` and `terms` are consumed by `build_plan` and cannot influence `instantiate` — including DOF, which derives from the plan's `Shake/Settle/Lincs` entries. Guard it as G8: for every reference recipe `r`, `instantiate(build_plan(r), r.with_ensemble(other))` is bit-identical to `instantiate(build_plan(r), r)`. This is the wall OpenMM never built between `Integrator.getTemperature()` and `AndersenThermostat`, and every later term will lean on it.

### Review 2 — GROMACS (mdp → grompp → tpr, gmxapi) perspective

## Verdict (3 sentences)
This plan is, in substance, "build `grompp` for GROMOS": one preprocessing stage (`build_plan`) resolves everything into a complete, inspectable plan the engine executes without re-deriving, and A7 (gromosXX must accept every `to_imd`) is a stronger guard than GROMACS had for its first decade. Its weak side is the input: it inherits a parser whose silent defaults are not gromosXX's — an `.imd` with no MULTIBATH keeps `TempBathParameters::default()` (τ=0.1, imd.rs:306-317), so `thermostat_on` is true at simulation.rs:488-493 and md.rs:362-367 while gromosXX runs no bath — and G7 as written (PLAN.md:594) would freeze that. Add presence-aware parsing with a grompp-style diagnostics echo before step 2, and per-term energy slots before `TermSpec` gains a second additive variant.

## Scorecard
**D1 — Modify.** Grouped, typed data is right (the flat 200-key mdp namespace is GROMACS's unfixable regret), but "defaults in one place" must mean *recipe* defaults over presence-aware parsing: MULTIBATH, PRESSURESCALE, PERTURBATION, DISTANCERES, POSITIONRES, ENERGYMIN become `Option<Block>` in `ImdParameters` (imd.rs:219-304 defaults them). `passthrough` must be loud: md.rs:646-668 re-parses GAMD/EDS out-of-band, so a passthrough GAMD block yields plain MD with no error.

**D2 — Agree.** `grompp`→`.tpr`→`mdrun` in miniature, on one condition: `AlgorithmSpec` is *fully resolved* (DOF, virial type, NSM, `four_pi_eps_i` as values), never today's `Option` + `unwrap_or(imd.x)` (algorithm_sequence.rs:1190-1194). Otherwise `instantiate` re-derives — the one thing `mdrun` never does.

**D3 — Agree**, plus `#[serde(deny_unknown_fields)]` on every variant so `Term("schnet", cutof=…)` is grompp's fatal "Unknown left-hand" rather than a dropped key. `gmxapi.modify_input` converged on string-key dicts anyway; the typed pyclass layer was the drift.

**D4 — Modify.** Additive providers are exactly `IForceProvider` (density fitting, QMMM-CP2K, PLUMED never replace classical forces; A8 matches). But GROMACS gives each provider its own energy term (`F_DENSITYFITTING`, `F_COM_PULL`, `F_EQM`); here `orchestrator_algorithm.rs` *overwrites* `special_total` and states "alongside position restraints isn't supported", while posres/distres live inside `Forcefield` (simulation.rs:329-410, md.rs:1140-1147). Two `TermSpec`s cannot coexist today.

**D5 — Agree.** A7 is the plan's best idea. Make the bundle something the run *writes*, not only reads (checkpoint lesson below).

**D6 — Agree.** GROMACS's renames (`nstxtcout`→`nstxout-compressed`) survived because grompp rewrote them with a note; each shim should print its `Recipe(...)` equivalent.

**D7 — Modify.** Fine as `convert-tpr`/`modify_input`, but gmxapi confined edits to parameter values because arbitrary tpr edits break invariants. Add `validate_plan(&[AlgorithmSpec]) -> Diagnostics` for the order invariants (Orchestrator after Forcefield, TemperatureCalculation before Barostat, EnergyCalculation last) and stamp `Simulation.plan_edited`.

**D8 — Agree.** Strict-xfail first is what FUTURE.md:518-521 demanded. Add to G6 a grep for `parse_file(` on the input path outside `run/` (a second IMD reader is how mdp drift began), and to step 0 one non-EM reference without MULTIBATH: of 41 reference `.in` files, 39 carry it and the 2 without are EM, so the hazard is invisible.

## What GROMACS learned the hard way that this plan should adopt
1. **`mdout.mdp` — echo effective values.** `Recipe.from_imd` returns `(recipe, diagnostics)`: blocks absent→defaulted, blocks passed through, values coerced (imd.rs:1072-1081 turns unparseable numbers into `0.0`/`0` silently, and short lines keep defaults). `md` writes `<name>.recipe.toml` beside `.tre`; `Simulation.recipe` exposes the same. The divergence table's "cause" column then comes from diffing recipes before energies.
2. **`gmx dump`.** `md @input x.imd @dump` prints recipe + plan (`AlgorithmSpec` JSON) and exits. This is the Rust-side A/B that step 0 declares impossible (PLAN.md:611), available at step 2 for free.
3. **Fatal on unknown keys, `-maxwarn` for the rest.** `deny_unknown_fields`; `build_plan` *errors* on a physics-bearing passthrough block (GAMD, EDS, LOCALELEV) unless `allow_passthrough=[...]` names it explicitly.
4. **The checkpoint embeds the tpr.** No checkpoint exists yet; define `.recipe.toml` + topology hash + plan as its future header now, so restart work never re-reads a possibly edited `.imd`.
5. **Schema versioning from day one** (`recipe_version = 1`). `tpx_version`/`tpx_generation` came late; `convert-tpr` exists to pay for it.

## What NOT to copy from GROMACS
1. **gmxapi's data-flow graph** (`mdrun(modify_input(read_tpr()))`, futures, `while_loop`, `subgraph`): right for ensemble workflows, heavyweight for `sim.step(10)`. Stay eager; REMD later is a list of recipes, not a graph. architecture.md's "don't foreclose" is enough.
2. **Sentinel-typed values** (`nstcalcenergy = -1` "auto"; `-1`/`0`/`no` meaning different things per key). `parallel: Auto|Serial|Parallel` is right; `thermostat_tau = -1.0` meaning off (simulation.rs:491) is the mdp disease — use `ensemble.thermostat: Option<…>`.
3. **Opaque binary run input.** Keep recipe and plan as text (serde TOML/JSON); never bincode the plan for speed. Every `gmx dump` is an apology for the `.tpr`.

## Risks the plan misses
1. Absent MULTIBATH → Berendsen bath on (above). A7 cannot catch it: the only bath-less references are EM, which skips the thermostat.
2. Energy-slot collision: moving posres into `Orchestrator{terms}` changes where `E_Posrest` lands relative to the `.tre` gromosXX writes; G2 compares front-ends of the *same* build and would pass while the reference comparison shifts.
3. Writer coverage is n=1 for PRESSURESCALE and DISTANCERES, n=4 for PERTURBATION; anisotropic pressure tensors, NHC chain lengths, LINCS orders reach `to_imd` untested. Step 2's gromosXX script must also run on `to_imd(Recipe.npt(...))` factory outputs, not only the round trips.
4. If `AlgorithmSpec` keeps `Option` fields, "Python edits the plan, never the instantiation" (PLAN.md:530-531) is false: an edited `Thermostat{}` and a `total_dof` re-derived from `recipe.constraints` can disagree.
5. `ntc_from_constraints` and the `nve/nvt/npt` factories on `ImdParameters` (imd.rs:932-956) are a fourth defaults site the audit at PLAN.md:469 does not list; step 3 must move them, not copy them.

## The one change you would insist on before implementation starts
Before step 2: make `ImdParameters` presence-aware (`Option` per optional block), have `from_imd` and `build_plan` return `Diagnostics` (notes / warnings / errors; `md` refuses to run on errors and on unacknowledged physics-bearing passthrough), and restate G7 as "`RunRecipe::default()` is what gromosXX does with the 13 blocks every reference carries — verified by running gromosXX on `to_imd(RunRecipe::default())`", not `from_imd(ImdParameters::default())`. Everything else in the plan gets stronger once the input side cannot lie.

### Review 3 — HOOMD-blue v3/v4 perspective

## Verdict (3 sentences)
The plan lands on the shape HOOMD reached only after its v2→v3 rewrite — declarative parameter objects, one attach point, topology-dependent validation at that point — and improves on HOOMD where it is still weak (a run reproducible from a data file, which HOOMD cannot do). Plan-as-data with no live editing is the better choice than HOOMD's attached `SyncedList` model for a project whose step order *is* the specification. What is unsettled is that "term" is not one contract in the code today (restraints are `Forcefield` fields; the orchestrator overwrites `special_total`), and an editable ordered plan has no validator — both must be fixed on paper before step 2.

## Scorecard
**D1 — Agree.** v2's `context.initialize()` global made construction equal registration and left defaults split between Python and C++; a topology-free `RunRecipe` with `from_imd(ImdParameters::default()) == default()` (PLAN.md:594) is exactly v3's fix, anchored in Rust.

**D2 — Modify.** Endorse plan-as-data over HOOMD's attach-time object creation (source of its shared-`nlist` detach/re-attach bugs), but it is three stages, not two: construction is entangled with `pairlist_algorithm.update(&topo,&conf,…)` before `Forcefield::new` (simulation.rs:527), and `init()` mutates positions under NTISHK (shake_algorithm.rs:54-73). Name `start = init + step 0` as a shared stage.

**D3 — Modify.** Equivalent to `TypeParameter`/`ParameterDict` validating at assignment — good — but: bridge with `pythonize`, not `serde_json` strings (JSON drops `inf`/NaN; errors carry byte offsets, not kwarg names); `#[serde(deny_unknown_fields)]` on every variant or `Term("schnet", cutof=…)` silently drops the typo — write that test first; an exhaustive-`match` registry yields names only, whereas HOOMD's validators carried *types* — emit `(name, type, default)` so G5 checks types.

**D4 — Modify.** This is `hoomd.md.Integrator(forces=[…])`, which works because every element is a `Force` with one accumulation contract. Here restraints are `Forcefield` fields (forcefield.rs:124-128, not `PotentialProvider`s) and `ProviderOrchestratorAlgorithm` *overwrites* `special_total`, the slot posres writes (orchestrator_algorithm.rs:23-30, forcefield.rs:1252) — `terms=[posres, schnet]` is unsupported today. See "one change".

**D5 — Agree.** HOOMD never solved this: GSD round-trips state, not operations, so no HOOMD run is reproducible from data. Note `raw_blocks` already stores modelled blocks verbatim (imd.rs:351-353); `to_imd` must choose per block "regenerate" vs "passthrough", and G3 should add a byte-level check for passthrough blocks.

**D6 — Agree.** HOOMD 3.0 broke everything with no shims and retaught its users; v4's `NVT → ConstantVolume(thermostat=…)` used one deprecation cycle and went smoothly. At 0.1.0 with two in-repo kwarg users, one cycle suffices — provided shims route through `Recipe` (G6), never a second path.

**D7 — Modify.** "No live editing" is the right lesson from `SyncedList`. But HOOMD's integrator *owns* the intra-step order (forces→methods→constraints); users physically cannot put a barostat before the pressure compute. Exposing the order means `instantiate` must validate it — `resolve_algorithm_sequence` checks nothing and the `pressure` getter documents the resulting garbage (simulation.rs:1107-1114). Also name-keyed `insert_after("Forcefield")` (algorithm_sequence.rs:700-716) is ambiguous once a kind repeats; HOOMD keys by object identity — add index/identity forms.

**D8 — Agree**, one gap: G1 (PLAN.md:571) names `prepare_system` + `build_sequence` but not `init`/step 0, which is physics. And the existing "parity" test compares only `.names` (test_gromosXX_references.py:598) — precisely FUTURE.md's trap; step 0's `array_equal` v0 must supersede it.

## What HOOMD learned the hard way that this plan should adopt
1. **Two validation moments, each with a path in the error.** Assignment-time (`TypeConversionError` naming `lj.params[('A','A')]['epsilon']`) and attach-time completeness against the state ("type B has no value"). Map: `Term()` → serde; `instantiate` → topology-dependent checks (region resolves to ≥1 atom, `elements` cover it, posres indices in range), with errors naming term and field rather than "Failed to resolve algorithm sequence".
2. **The engine owns ordering invariants even when the list is editable.** Declare them as registry data (`Orchestrator.must_follow = Forcefield`; `Barostat.requires = PressureCalculation`; `SteepestDescent.excludes = {LeapFrog*, Thermostat, TemperatureCalculation}`) and run `validate_plan` inside `instantiate`.
3. **Per-term energy namespace** (`Force.energy`; `hoomd.logging.Logger` loggables like `md.pair.LJ.energy`) rather than a fixed table. `run()` returns 12 fixed columns (simulation.rs:927-935); a term's energy lands in `special_total`, invisible from Python. Give each `TermSpec` an energy slot keyed by its registry name.
4. **`run(0)` was explicit**; the gotcha was `lj.energies` raising before it. Eager step 0 avoids the gotcha but means nothing observes frame 0 unless `outputs` is instantiated with the plan — make writers part of `instantiate`, shared with `md.rs`.
5. **Ship the migration table in the same PR as the shims**; HOOMD's v3 guide arrived after and users reconstructed it from tracebacks.

## What NOT to copy from HOOMD
1. **`hoomd.custom.Action`/`CustomForce`** — Python-side state mutation; violates zero-physics-in-Python and makes `md` parity impossible. Any future hook must be a read-only observer.
2. **Live mutation after attach** (`sim.operations.integrator.forces.append` mid-run). pyo3-gromos CONTEXT.md still lists `Simulation.add_ml_potential` — don't; it recreates the kwarg problem as a method.
3. **`Trigger`/`filter` class hierarchies on every operation.** GROMOS frequencies are IMD integers (`NSCM`, `NSNB`, `NTWE`); keep them integers so `to_imd` stays trivial — Trigger objects are one reason HOOMD operations cannot be serialized.

## Risks the plan misses
1. **`np.array_equal` (A3) cannot hold under `parallel`.** Chunking uses `current_num_threads()` and rayon's fold/reduce tree depends on work-stealing (innerloops.rs:576-607). G2 must pin `parallel: Serial`, or the reduction must be made thread-count-independent; A2's "pinned" is insufficient.
2. **Plan names vs instantiated names.** `LeapFrogIntegrator` expands one spec into two algorithms (algorithm_sequence.rs:1265-1268); `AlgorithmSpec::name()` and `Algorithm::name()` will disagree unless one spec = one algorithm. G2's `sequence=` path should assert equality.
3. **MULTIBATH truncation.** Every path reads `temp_bath[0]` only (simulation.rs:483-492). A scalar `ensemble.thermostat {T, tau}` makes `to_imd` silently drop extra baths; model a list or make `to_imd` refuse. (No current reference has >1 bath, so G3 will not catch it.)
4. **`#[serde(default)]` policy.** Needed on recipe groups for `from_dict` leniency, wrong on `TermSpec` variants where a missing `model` must fail. Decide per struct; don't derive uniformly.
5. **Thin restraint oracle.** Only `nacl_1water_distres` exercises restraints; any move of restraints between `Forcefield` and providers has almost no guard.

## The one change you would insist on before implementation starts
Define the term contract before writing `TermSpec`: every `TermSpec` instantiates to exactly one `PotentialProvider` registered in the orchestrator with its own energy slot, and `Forcefield` either stops owning restraints (ported to providers, parity-tested on `nacl_1water_distres`) or the plan states that restraints are `AlgorithmSpec::Forcefield` parameters and *not* terms. The text currently claims both (PLAN.md:518 vs forcefield.rs:124), the `special_total` overwrite (orchestrator_algorithm.rs:89) makes them mutually exclusive, and step 5's "two files touched" metric is meaningless until this is settled.

### Review 4 — ASE / i-PI / PLUMED / OpenMM-ML (pluggable force providers) perspective

## Verdict (3 sentences)
The skeleton — plain-data `RunRecipe` → `build_plan` → `instantiate`, serde-tagged terms, one builder under two oracles — is the right shape and cleaner than what ASE or OpenMM-ML converged on for the same problem. Its one serious flaw is that `TermSpec` (PLAN.md:518-521) says *which* provider to run but not *what it owns*, so `Term("schnet"|"xtb"|"mopac")` under A8 (PLAN.md:559) silently yields a non-GROMOS ML/MM (inner zone computed twice, admitted at `ml_potential.rs:145-150`) while the type reads like the QM/MM this project already reproduced to 1.3e-9 (PLAN.md:342-344). Fix ownership and the energy channel before step 2; the rest is refinement.

## Scorecard
**D1 — Modify.** Grouping is right, but `terms` mixes two ownership models: `PositionRestraints`/`DistanceRestraints` live *inside* `Forcefield` today (`forcefield.rs:1247-1276`, writing `special_total`/`distanceres_total` at 1e-8 parity) while SchNet/Xtb go through an orchestrator that *overwrites* `special_total` (`orchestrator_algorithm.rs:89`, docs 23-30). Either restraints stay `Forcefield` config (then "one variant + one arm" is false for them) or they become providers and the two clobber each other. Also: gromosXX already has a file model for QM terms (`QMMM` IMD block + `.qmm` with `QMZONE`/`QMUNIT`) and `gromos-io` reads none of it — `passthrough` (PLAN.md:516) would carry it opaquely, so `to_imd` of a QM term is not round-trippable and A7 quietly narrows to classical runs.

**D3 — Agree**, with two must-haves: `#[serde(deny_unknown_fields)]` on every variant (a typo'd `cutof=` otherwise vanishes — ASE's `Calculator(**kwargs)` accepted anything for years), and term paths resolved against the bundle directory, not CWD. Exhaustive-match registries are OpenMM-ML's `registerImplFactory` done at compile time — better.

**D4 — Disagree** (with the framing, not the deferral). "Additive" is correct for restraints, `LocalElevation` (a PLUMED-style bias) and Δ-learned models — LAMMPS `hybrid/overlay`. It is wrong for the three engines listed: the project's own training target is the *isolated* QM energy (PLAN.md:359-361). OpenMM-ML `createMixedSystem` and LAMMPS `hybrid` (not overlay) exist because this is the default trap. Deferring zone-aware `Forcefield` is fine only if the term *declares* its coupling and `instantiate` refuses the unsupported one — the "reject, never silently drop" rule already applied to `Embedding`.

**D5 — Modify.** The `[input]` shape (`ch4_water_fep/input.toml`) works for paths; add provenance — model checksum, libtorch/xtb/MOPAC version. The MOPAC-11% lesson (PLAN.md:335-337) is the archetype. `elements` must be in the bundle: `Topology` has no atomic numbers (`ml_potential.rs:10-14`); gromosXX keeps them in `QMZONE`.

**D7 — Agree.** ASE's `atoms.calc` + stateful calculators + `calculation_required` caching is the counterexample; PLUMED requires restart for input changes. Eager construction also surfaces "xtb not on PATH"/bad model at `Simulation()`, not step 1.

**D8 — Modify.** G1-G7 are good; step 5's metric (PLAN.md:673-677) is wrong. "Two files" proves wiring, not physics. `Term("xtb")` should have to reproduce `water_dimer_mechst` — which it *cannot* under A8, since that oracle was hit with `LjCrfInteraction::with_zone_partition`, a replacing scheme. Either step 5 is `Term("xtb", coupling="delta")` with the discrepancy stated, or 2.8 moves before step 5.

## What the calculator/plugin ecosystems learned the hard way that this plan should adopt
1. **Units are input data, not constants.** gromosXX's `QMUNIT` block carries per-engine length/energy/force/charge factors; PLUMED makes the host declare units (`setMDEnergyUnits`/`setMDLengthUnits`); i-PI forces atomic units on the wire; ASE's CODATA-2014 switch (`ase.units.create_units`) shifted every calculator's numbers. Here `HARTREE_TO_KJMOL = 2625.4996` (8 figures) sits beside `BOHR_TO_NM = 0.0529177210903` (12) — different vintages, invisible to the 1e-8 suite because it contains no QM run. Put `units` on engine terms with `QMUNIT` defaults. SchNet assumes nm→kJ/mol (`schnet.rs:200,266`); an eV/Å model is off by 96.5/964.9 silently — SchNetPack's own `SpkCalculator` takes `energy_unit`/`position_unit` for exactly this. `schnet.rs:200` also casts positions to f32: declare it.
2. **Ownership is part of the term.** LAMMPS `hybrid` vs `hybrid/overlay` is the additive/replace choice at input level; OpenMM-ML `interpolate=True` keeps both under λ. `zones.rs` already has the vocabulary (`PairOwner`, `owner()`, `lj_owner()`, `qm_lj`); the term only needs to name which table it claims.
3. **External-process lifecycle.** `qm_subprocess::run_subprocess` spawns per step, no timeout; `xtb.rs:292-295` uses fixed filenames in `work_dir` — two `Term("xtb")` entries, or `XtbPotential.evaluate()` from two threads, corrupt each other (ASE `FileIOCalculator` `directory` collisions; i-PI gives each client its own scratch and a `timeout`). Derive `work_dir` from the term index, add a timeout, and use the `&mut self` the trait grants for SCF restart (gromosXX keeps the xtb calculator alive; cold SCF per step also changes the drift figure at PLAN.md:349).
4. **Per-term energy channel.** ASE's `SumCalculator` lost component energies; PLUMED reports `.bias` per action; the `.tre` reports posres/distres/LE separately. `Contribution.energy` is one scalar summed at `orchestrator.rs:96` — with `Vec<TermSpec>`, G2 parity can only compare totals, and compensating bugs pass.
5. **Failure classes.** i-PI separates transient (retry/reconnect) from fatal; `ProviderError` has two variants and `run_step` (`algorithm.rs:181`) propagates a `String` and loses the frame. Write the current configuration before propagating.

## What NOT to copy
1. **ASE's atoms-hash caching (`check_state`).** One evaluation per step, providers owned by the sequence — caching is what produced stale-result bugs.
2. **i-PI's socket as default transport.** In-process `tch` and subprocess are simpler; sockets pay off only for a persistent remote engine — a future provider, not the seam.
3. **OpenMM-ML's "edit the MM system" surgery.** It mutates classical force objects to remove pairs; here that would mean Python editing `Forcefield`. The static ownership table in `zones.rs` plus 2.8 is the better design.

## Risks the plan misses
1. `special_total` exclusivity vs `forcefield.rs:1252`: a recipe with `PositionRestraints` + `SchNet` is silently wrong today and no reference system has both — make it a *required* row in G5's table.
2. Virial: `XtbInteraction` returns `Mat3::ZERO` (`xtb.rs:359`). `Term("xtb")` in an NPT recipe gives wrong pressure with no error; `instantiate` should reject `Barostat` + a provider lacking virial — an `implemented_properties`-style flag on the trait is the one ASE idea worth taking.
3. Regions are resolved once; PLAN.md:434's proximity zones would make ownership step-dependent. State now that regions are static per run.
4. G2's `np.array_equal` (A3): f32 tensors and libtorch's thread pool may not be bit-reproducible — pre-register this finding.
5. G4's "one example value" must be a *runnable* example, or `test_every_kind_has_a_parity_case` is satisfied by dummies for `charge`/`multiplicity`/`gfn`/MOPAC method.

## The one change you would insist on before implementation starts
Add `coupling` to every engine `TermSpec` — `Delta` (additive, what exists) | `Replace { qm_lj }` (claims `PairOwner::Provider` pairs; `instantiate` returns "requires zone-aware Forcefield (2.8)" until it lands) — and give `Contribution`/the orchestrator a per-term energy so `.tre` output and G2 see each term. Step 5's proof then becomes `Term("xtb", coupling="replace")` matching `water_dimer_mechst`, not a file count.

### Review 5 — Polars / pyo3 (Rust↔Python plan-as-data) perspective

## Verdict
The shape is right: `RunRecipe → build_plan → Vec<AlgorithmSpec> → instantiate` is Polars' `DslPlan → IR → Executor` split, and the audit (PLAN.md 443-483) names the real disease. Two signatures betray it: `instantiate(&[AlgorithmSpec], &RunRecipe, …)` (l.528) lets stage 2 read past the plan, and `recipe.terms.append(Term(...))` (l.657) assumes a mutability model pyo3 does not give you. Fix those, put defaults and factories in the library rather than the binding, and the rest is execution.

## Scorecard
**D1 — Modify.** Grouped plain data: yes. "Defaults in one place" is already false by construction: `ImdParameters::default()` (imd.rs 219-304) *is* a defaults table and the factories (imd.rs 958-1049) are a second. Make `RunRecipe::default()` *derived* (`from_imd(&ImdParameters::default())`), not a third table G7 checks for equality. Step 3 (l.649) moves factories into `pyo3-gromos` — the thin layer would own defaults, exactly how Polars' eager `read_*`/lazy `scan_*` defaults drifted for years (Python signatures and Rust `Default` impls both held them). Factories are `RunRecipe::nvt(...)` in the library. On A4: take the `gromos-run` fallback now — `gromos-md` carries clap/env_logger and optional `mpi`/`cudarc` (Cargo.toml 21-34), and splitting is cheap before the code exists.

**D2 — Modify.** `build_plan(&RunRecipe, &Topology)` is DSL→IR resolution against the "schema" — correct. But Polars builds executors from IR *only*; `instantiate` re-reading `&RunRecipe` (l.528) means `AlgorithmSpec` is not self-contained — the `unwrap_or(imd.x)` layering the audit condemned (l.471-473) under a new name.

**D3 — Modify.** Tagged-enum serde as the bridge: agree — it is how `PyExpr`/`PyLazyFrame` pickle. Conditions: `#[serde(deny_unknown_fields)]` on both enums plus a test that a typo'd kwarg fails (serde silently drops unknown fields; the internally-tagged path buffers through `Content` and does not compose with `flatten`, which `passthrough` will tempt you toward); `Term(kind, **params)` via `pythonize`/`depythonize` (pyo3 0.21-compatible), not a JSON-string detour; exhaustive-`match` registries — agree, add `#![deny(clippy::wildcard_enum_match_arm)]` in `run/`. The `.pyi` substring test is too weak (D8).

**D5 — Modify.** `to_imd`/`from_imd`: agree, that is the stable format. `to_dict`: Polars learned at 1.0 that a visible JSON plan becomes a format users depend on, switched the default to binary and declared plans same-version-only. Declare it now: a `version` field, `#[serde(default)]` on every group (partial dicts load; new fields don't break old dicts), `BTreeMap`/`IndexMap` instead of `raw_blocks: HashMap` (imd.rs 195) so `to_imd`/`to_dict` ordering is deterministic, and note `serde_json` writes NaN as `null` and cannot read it back.

**D6 — Modify.** One cycle at 0.1.0 is right (Polars purged its 0.x deprecations in one 1.0 break). But shims live in Python (`polars/_utils/deprecation.py`; `PyLazyFrame` never carried legacy signatures) — the plan puts them in pyo3 (l.656-658), so step 4 deletes Rust instead of a Python file. Add `filterwarnings = ["error"]` to `[tool.pytest.ini_options]` (pyproject 74-81) or `DeprecationWarning` stays silent outside `__main__`. `__init__.py:73` hardcodes `0.1.0` while the workspace is `0.0.26` — derive from package metadata.

**D7 — Modify.** Polars deliberately offers no IR editing from Python. If you do, validation is the price: ordering invariants live in comments today (simulation.rs 578-580, 619, 634); `validate_plan(&[AlgorithmSpec]) -> Result<(), RunError>` must precede `instantiate`, or the escape hatch itself violates "Rust owns invariants". Keep `PyAlgorithmSequence` a `Clone` plain-data pyclass; only `PySimulation` needs `unsendable` (l.784).

**D8 — Modify.** G1-G4 and G7's equality check: agree. G5: replace the grep with `python -m mypy.stubtest gromos` — it checks a hand-written `.pyi` against the runtime module (needs `#[pyo3(signature=…)]` everywhere so `__text_signature__` exists). G6: `just lint` is only `cargo clippy` (Justfile 24-25), and `AlgorithmSequence::new()` has six test hits in `gromos-core/src/algorithm.rs` (250-394) — exclude `#[cfg(test)]`. G7's literal grep for `0.5`/`-1.0` matches physics everywhere; scope it to `pyo3-gromos/src` and `bin/`. Steps: Step 0 first is right; Step 1's "verbatim" move is 19 `println!`/`process::exit` sites (md.rs 414-572, 987-1298) that must become `Result` — decide the error type there. Rust-side A/B is declared impossible at step 0 (l.611) and never added: put `cargo test -p gromos-md` parity in step 2's gate so G1 does not depend on maturin.

## Adopt from Polars/pyo3
1. **`polars-python` layout, recipe below it.** `pyo3-gromos` already is `polars-python`; `RunRecipe`, factories and `RunError` go in the library with `serde` non-optional — Polars gates serde on `polars-plan`, then py-polars must always enable it; the gate buys nothing.
2. **Error newtype.** Orphan rules forbid `impl From<gromos_run::RunError> for PyErr` inside `pyo3-gromos`; Polars' answer is `PyPolarsErr(PolarsError)` + one `From` + `create_exception!` (`polars.exceptions.ComputeError`). Today: 50 ad-hoc `PyErr::new::<Py*Error>` sites, `Algorithm::apply -> Result<(), String>`, and A14 "asserts the error text". Define `RunError {MissingFeature{term, feature}, UnknownKind, UnknownParameter, InvalidPlan, Io}` and `gromos.exceptions`.
3. **Pickle via serde.** `__getstate__/__setstate__` on `Recipe` and the plan (`PyLazyFrame` does exactly this); `multiprocessing` drivers need it.
4. **Feature flags.** A14 is right (Polars cfg-gates `Expr` variants, making serialized plans build-dependent). Add `gromos.build_info()` and `available: bool` per registry entry so Python can tell `schnet` is data-only before `instantiate`; CI builds `--no-default-features` and `--all-features` because exhaustive matches differ under cfg.
5. **Immutable builders.** If `.terms` is a `#[getter]` returning `Vec`, `recipe.terms.append(...)` mutates a temporary. Polars avoids nested mutation by construction — every op returns a new object: `Recipe.with_term(...)`, `.with_control(dt=...)`.

## Do NOT copy
1. cfg-gated enum variants in `Expr`/`DslPlan` — keep data unconditional (A14).
2. `Arena<IR>`, optimizer, lazy `collect()` — the GROMOS step order is fixed by the reference; `Vec<AlgorithmSpec>` and a stateful `Simulation` are right.
3. The pure-Python mirror layer (`polars/lazyframe/frame.py` restating every signature) — with ~60 recipe fields that layer *is* the drift; keep pyclasses direct, `Term(kind, **params)` generic, Python only for deprecation and sugar.

## Risks missed
1. Recipe expressiveness < `ImdParameters` for *modelled* blocks: MULTIBATH is `Vec<TempBathParameters>` (imd.rs 63, 200-207) yet every builder reads `temp_bath[0]`; `PressureParameters` has 3×3 tensors (213-214). G3 fails on such systems — and a multi-bath `.imd` already runs different physics than gromosXX (a divergence-table row). Step 2 needs a field-by-field table: modelled / passthrough / rejected, for all ~90 fields.
2. `passthrough` is a physics hole: a block gromos-rs ignores but gromosXX honours passes A7 (gromosXX reproduces *its own* energies) while the recipe runs something else. Passthrough must be an allowlist of inert blocks; unknown blocks are errors.
3. `AlgorithmSpec` is never stated to derive serde/`PartialEq` or carry a version, yet `Algorithm(kind, **params)` (l.650) requires it.
4. `pyo3-gromos/ml` forwards only `gromos-forces/ml` (Cargo.toml 30-32); it must also forward `gromos-md/ml` or `SchNet` instantiates in one crate and errors in the other.
5. `.local/polars`, the "reference" both CONTEXT.md files cite, is not in the checkout — the pattern is unverifiable in review.

## The one change
`instantiate(&[AlgorithmSpec], &Topology, &Configuration) -> Result<AlgorithmSequence, RunError>` — no `&RunRecipe`. Everything stage 2 needs (resolved DOF, cutoffs, box kind, restraint paths, `four_pi_eps_i`, the `parallel` decision) is written into the spec by `build_plan`, so an edited plan is exactly what runs; `RunRecipe` feeds `build_plan` and a `RunControl {dt, nstlim, t0}` only. Polars' correctness story rests on "executors see only IR"; without it G2's third front-end tests nothing and step 5's two-file count quietly becomes four.

