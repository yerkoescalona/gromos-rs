# gromos-rs — Progress Dashboard
> Auto-generated 2026-06-16 from PLAN.md audit

State legend:
- ✅ **DONE** — implemented, wired, reference-tested
- 🔶 **PARTIAL** — code exists, but unwired / disconnected / incomplete
- 🔧 **STUB** — code exists with wrong/placeholder logic (needs rewrite)
- ❌ **TODO** — not yet started

---

## Reference Test Suite
```
34 / 34  ████████████████████████████████████████  100%  ALL PASS
```

---

## Priority 1 — MD Engine Physics (gromosXX-faithful)

| #   | Item                                              | State    | Notes |
|-----|---------------------------------------------------|----------|-------|
| 1.1 | NTIVEL velocity gen + unit-conversion audit       | ✅ DONE  | bit-for-bit, ref-tested |
| 1.2 | SETTLE, LINCS, COM rotation removal               | ✅ DONE  | all ref-tested |
| 1.3 | Nosé-Hoover thermostat (single + chain NHC)       | ✅ DONE  | ref-tested |
| 1.4 | Triclinic / truncated-octahedron boundary         | ✅ DONE  | ref-tested (`aladip_trunc_oct`) |
| 1.5 | O(N) cell-list pairlist                           | 🔶 PARTIAL | built + set-equality-validated; **not wired** into `md` binary |
| 1.6 | Restraints (distance, dihedral, angle, J-val, order-param, dist-field) | ❌ TODO | 0/6 — no code yet |
| 1.7 | Virtual atoms (`addvirt_top`)                     | ❌ TODO  | no code yet |
| 1.8 | Advanced sampling (EDS, GaMD, FEP/TI, REMD)      | 🔶 PARTIAL | **code exists** in integrators/forces (eds.rs, gamd.rs, fep.rs, remd.rs); **untested, unvalidated** |

```
P1 overall   ██████████░░░░░░░░░░░░░░░░░░░░░░  14 / 21   67%
             (counting 1.8 as ½ credit — code exists but zero ref tests)
```

---

## Priority 2 — Analysis Foundations (no-duplication layer)

| #    | Item                                                        | State       | Notes |
|------|-------------------------------------------------------------|-------------|-------|
| 2.1a | Queryable per-atom metadata (molecule map, residue/name incl. solvent) | 🔶 PARTIAL | basic `AtomSelection` exists; `parse_solvent` ignores spec; `parse_molecule` rejects non-first |
| 2.1b | `AtomSpecifier` gromos++ facade                             | ❌ TODO     | not started |
| 2.2a | Rotational fit (Kabsch/SVD)                                 | 🔧 STUB     | `simple_rotation_fit` in `rmsd.rs:161` only re-centers — **no rotation** |
| 2.2b | Statistics + block-averaging error (`ee()`)                 | ❌ TODO     | not started |
| 2.2c | PBC gathering / molecule unwrapping (unified primitive)     | ❌ TODO     | not started |
| 2.2d | Single-point energy entry point (`ener` → `gromos-forces`)  | 🔧 STUB     | `ener.rs` **hardcodes LJ σ/ε** instead of calling engine |
| 2.3a | `rmsd` — real rotational fit                                | 🔧 STUB     | runs, but re-centering only; blocked on 2.2a |
| 2.3b | `ext_ti_ana` — real dH/dλ parsing                          | 🔧 STUB     | falls back to synthetic data; blocked on 2.1+2.2 |
| 2.3c | `nhoparam` — NMR N-H order parameters S²                   | 🔧 STUB     | **wrong algorithm** — currently computes Nosé-Hoover params |
| 2.3b | Free-energy estimators (BAR, ext_ti_merge, reweight, Widom, dg_ener) | 🔶 PARTIAL | code scaffolds exist; no reference tests |
| 2.4  | Other stubs: `visco`, `frameout`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy` | 🔧 STUB | code exists but placeholder/wrong logic |

```
P2 overall   ████░░░░░░░░░░░░░░░░░░░░░░░░░░░░   code exists throughout
             0 items fully done (ref-tested); ~8 partial/stub; 2 not started
```

---

## Priority 3 — py-gromos API & Education

| Phase | Item                                                         | State       | Notes |
|-------|--------------------------------------------------------------|-------------|-------|
| 3.1 ✅ | Simulation API, AlgorithmSequence API, 62 Python ref tests, `.pyi` stubs | ✅ DONE | |
| 3.1 ✅ | Topology/Configuration/Simulation wrappers + numpy interop   | ✅ DONE     | |
| 3.1   | Expose ForceField eval (single-point energy/force)           | ❌ TODO     | |
| 3.1   | Expose SHAKE / constraint info                               | ❌ TODO     | |
| 3.1   | Expose energy decomposition (bonded, LJ, CRF, kinetic, pressure) | ❌ TODO | |
| 3.1   | Study Polars pyo3 patterns for API design                    | ❌ TODO     | |
| 3.2   | Method chaining: `sim.run().energies().plot()`               | ❌ TODO     | |
| 3.2   | Energy timeseries as DataFrame (Polars/pandas)               | ❌ TODO     | |
| 3.2   | `md_runners.py` simplify; `analysis.py` Python-expose        | 🔶 PARTIAL  | files exist, need update |
| 3.2   | Rich `__repr__` / `_repr_html_` for Jupyter                  | ❌ TODO     | |
| 3.3   | Rewrite notebooks (01/02/03) on new API                      | 🔶 PARTIAL  | old notebooks exist, need rewrite |
| 3.3   | Rewrite 17 example scripts on new API                        | 🔶 PARTIAL  | scripts exist, need update |
| 3.3   | Fix `test_basic.py`, `test_advanced_features.py`             | 🔧 STUB     | tests exist but broken |

```
P3 overall   ████░░░░░░░░░░░░░░░░░░░░░░░░░░░░   2 items fully done (Simulation API + wrappers)
             rest: partial / needs update / broken
```

---

## Priority 4 — Code Quality & Consistency

| Item                                                        | State       | Notes |
|-------------------------------------------------------------|-------------|-------|
| Clippy (~390 warnings: forces 89, integrators 77, io 31, core 15) | 🔶 PARTIAL | warnings known, not addressed |
| Replace bare `unwrap()` with `.expect()` / `?`             | 🔶 PARTIAL  | present in non-test code |
| Add missing `#[test]` (SHAKE — currently 0, improper dihedral) | 🔶 PARTIAL | gaps known |
| Split large files (nonbonded ~1500, bonded ~1300, topology ~1200 LOC) | 🔶 PARTIAL | oversized, not split |
| Unify CLI error types; audit `pub` visibility               | ❌ TODO     | |
| Benchmarking infra (baseline + CI regression)               | 🔶 PARTIAL  | `nonbonded_bench`, `math_bench` etc. exist; no saved baseline |

```
P4 overall   ██░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░   infrastructure exists, not addressed
```

---

## Overall Picture

```
                     fully done   partial/stub   not started
Priority 1 (engine)     62%          19%            19%
Priority 2 (analysis)    0%          82%            18%   ← code is there, needs fixing
Priority 3 (py-gromos)  18%          36%            46%   ← base layer done
Priority 4 (quality)     0%          83%            17%   ← known gaps, not addressed
```

**The real situation:** The codebase is much further along than a raw 22% suggests.
Most analysis tools, py-gromos bindings, and the advanced sampling code exist as working
(if imperfect) implementations. What's missing is **reference tests and correct algorithms**
in specific places — not blank files.

---

## Nearest Next Steps (unblocked)

1. **Wire `CellListPairlistAlgorithm` into `md.rs`** — closes P1.5, one heuristic remaining
2. **Distance restraints** (P1.6 first item) — unblocked, port `distance_restraint_interaction.cc`
3. **Benchmarking baseline** — `cargo bench --save-baseline v0.1` before any perf work
4. **Fix `AtomSelection` gaps** (P2.1a) — small, high-leverage, unblocks all of P2
5. **Kabsch/SVD rotational fit** (P2.2a) — fixes `rmsd`, unblocks `nhoparam` rewrite

---

## Parked / Blocked

| Item               | Blocked by                          |
|--------------------|-------------------------------------|
| `ext_ti_ana`       | P2.1 selection + P2.2 fit + stats   |
| `nhoparam`         | P2.2a fit (needs full algo rewrite) |
| P2.3b free-energy estimators | P2.3b `ext_ti_ana` + P2.2 stats |
| Tutorials t01–t06 end-to-end | compute-limited; deferred to last |
