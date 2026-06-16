# gromos-rs — Progress Dashboard
> Auto-generated 2026-06-16 from PLAN.md audit

---

## Reference Test Suite
```
34 / 34  ████████████████████████████████████████  100%  ALL PASS
```

---

## Priority 1 — MD Engine Physics (gromosXX-faithful)

| #   | Item                                         | Status | Progress                          |
|-----|----------------------------------------------|--------|-----------------------------------|
| 1.1 | Reproducibility (NTIVEL velocity gen + unit audit) | ✅ DONE | `████████████` 2/2 |
| 1.2 | Constraints (SETTLE, LINCS, COM rotation)    | ✅ DONE | `████████████` 3/3 |
| 1.3 | Nosé-Hoover thermostat (single + chain NHC)  | ✅ DONE | `████████████` 1/1 |
| 1.4 | Triclinic boundary (truncated octahedron)    | ✅ DONE | `████████████` 4/4 |
| 1.5 | O(N) cell-list pairlist                      | 🔶 PARTIAL | `█████████░░░` 3/4 — built+validated, **not wired** into `md` binary |
| 1.6 | Restraints (distance, dihedral, angle, J-val, order-param, dist-field) | ❌ TODO | `░░░░░░░░░░░░` 0/6 |
| 1.7 | Virtual atoms                                | ❌ TODO | `░░░░░░░░░░░░` 0/1 |
| 1.8 | Advanced sampling (EDS, GaMD, FEP/TI, REMD)  | ❌ TODO | `░░░░░░░░░░░░` 0/4 |

```
P1 overall   ████████░░░░░░░░░░░░░░░░░░░░░░░░  13 / 21   62%
```

---

## Priority 2 — Analysis Foundations (no-duplication layer)

| #    | Item                                                      | Status   | Progress |
|------|-----------------------------------------------------------|----------|----------|
| 2.1  | Atom selection (per-atom metadata + AtomSpecifier facade) | ❌ TODO  | `░░░░░░░░░░░░` 0/3 |
| 2.2  | Shared infra (Kabsch fit, stats/ee, PBC gather, single-point energy) | ❌ TODO | `░░░░░░░░░░░░` 0/4 |
| 2.3  | Fix blocked programs (rmsd rotfit, ext_ti_ana, nhoparam)  | ❌ TODO  | `░░░░░░░░░░░░` 0/3 |
| 2.3b | Free-energy estimators (BAR, TI merge, reweight, Widom, dg_ener) | ❌ TODO | `░░░░░░░░░░░░` 0/5 |
| 2.4  | Clean up known stubs (visco, frameout, amber2gromos, ener, sasa_hasel, dssp, solute_entropy) | ❌ TODO | `░░░░░░░░░░░░` 0/7 |

```
P2 overall   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░   0 / 22    0%
```

---

## Priority 3 — py-gromos API & Education

| Phase | Item                                                  | Status  | Progress |
|-------|-------------------------------------------------------|---------|----------|
| 3.1   | Rust bindings: ForceField eval, SHAKE, energy decomp, Polars study | ❌ TODO | `░░░░░░░░░░░░` 0/4 |
| 3.2   | Python API: method chaining, DataFrame, md_runners, repr | ❌ TODO | `░░░░░░░░░░░░` 0/4 |
| 3.3   | Notebooks & examples (17 scripts, 3 notebooks, tests) | ❌ TODO | `░░░░░░░░░░░░` 0/3 |

```
P3 overall   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░   0 / 11    0%
```

---

## Priority 4 — Code Quality & Consistency

| Item                                              | Status  |
|---------------------------------------------------|---------|
| Clippy (~390 warnings across workspace)           | ❌ TODO |
| Replace bare `unwrap()` with `.expect()` / `?`   | ❌ TODO |
| Missing `#[test]` (SHAKE, improper dihedral)      | ❌ TODO |
| Split large files (nonbonded, bonded, topology)   | ❌ TODO |
| Unify CLI error types; audit `pub` visibility     | ❌ TODO |
| Benchmarking infra (baseline + CI regression)     | ❌ TODO |

```
P4 overall   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░   0 / 6     0%
```

---

## Overall Plan Progress

```
Priority 1   ████████░░░░░░░░░░░░░░░░░░░░░░░░   13 / 21   62%
Priority 2   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░    0 / 22    0%
Priority 3   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░    0 / 11    0%
Priority 4   ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░    0 /  6    0%
─────────────────────────────────────────────────────────────
TOTAL        ████░░░░░░░░░░░░░░░░░░░░░░░░░░░░   13 / 60   22%
```

---

## Nearest Next Steps (unblocked, highest leverage)

1. **Wire `CellListPairlistAlgorithm` into `md.rs`** (P1.5 finish) — already implemented + validated
2. **Distance restraints** (P1.6 first item) — unblocked, start with `distance_restraint_interaction.cc`
3. **Benchmarking baseline** (P4) — `cargo bench --save-baseline` before any more perf work
4. **Atom selection solidification** (P2.1) — unblocks all of P2

---

## Blocked / Parked

| Item           | Blocked by              |
|----------------|-------------------------|
| `ext_ti_ana`   | P2.1 selection + P2.2 fit + stats |
| `nhoparam`     | P2.1 + P2.2 (needs full algorithm rewrite) |
| Tutorials t01–t06 end-to-end | compute-limited; deferred to last |
| EDS/GaMD/REMD (P1.8) | optional; big, untested code exists |
