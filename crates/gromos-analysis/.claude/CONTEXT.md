# gromos-analysis — stage contract

## Job
L4 facade. The gromos++-style analysis programs (rmsd, rdf, ene_ana, hbond, ...) as a thin
layer over the shared engine core. The facade must contain zero physics.

## Inputs (consumes from)
`gromos-core` + `gromos-io` today. Will add `gromos-forces` for single-point energy (P2.2 — the no-duplication change).

## Outputs
Analysis binaries producing results comparable to gromos++ programs.

## Status — known issues
- `rmsd`: works but `simple_rotation_fit` only re-centres, no rotational fit → P2.2 (replace with Kabsch/SVD). Ref: `.local/gromosPlsPls/gromos++/src/fit/RotationalFit.cc`
- `ener`: hardcodes LJ σ/ε instead of calling gromos-forces → P2.2 (single-point energy entry point)
- `nhoparam`: **wrong algorithm entirely** — computes Nosé-Hoover thermostat params, not NMR N-H order parameters S² → P2.3 (rewrite to gromos++ algorithm). Ref: `.local/gromosPlsPls/gromos++/programs/nhoparam.cc`
- `ext_ti_ana`: falls back to synthetic data → blocked on P2.1 + P2.2
- Stubs needing real implementations (P2.4): `visco`, `frameout`, `amber2gromos`, `sasa_hasel`, `dssp`, `solute_entropy`
- **No reference tests yet** — target of P2 roadmap

## Key files
```
src/bin/structural/rmsd.rs   — simple_rotation_fit (P2.2: replace with Kabsch)
```

## Crate-specific rules
- **Zero physics re-implementation.** If you're about to write LJ / CRF / bonded code here, stop — it belongs in gromos-forces. This is exactly how gromos++ accumulated 79k LOC.
- **The trap:** adding analysis logic that works but forks a physics term from the engine is the failure mode the whole architecture exists to prevent.
- **AtomSpecifier:** build as a facade over gromos-core primitives once P2.1 lands; do not build a parallel selection system.
- **Rotational fit** (Kabsch/SVD) belongs in gromos-core or a shared utility, consumed by both `rmsd` here and any future engine use.
