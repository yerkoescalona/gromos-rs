# gromos-rs — workspace guide

gromos-rs ports **gromosXX/md++** (MD engine) + **gromosPlsPls/gromos++** (analysis + tools) into
one Rust workspace with a single shared core. The architectural invariants live in `architecture.md`.
The porting discipline lives in `porting.md`. The long-horizon bets live in `FUTURE.md`. The live
roadmap and reference test matrix live in `PLAN.md`.

## When working in a crate, read its stage contract

| Crate | Layer | Stage contract |
|-------|-------|----------------|
| gromos-core | L0 data core | `crates/gromos-core/.claude/CONTEXT.md` |
| gromos-io | L0 I/O service | `crates/gromos-io/.claude/CONTEXT.md` |
| gromos-forces | L1 pure compute | `crates/gromos-forces/.claude/CONTEXT.md` |
| gromos-integrators | L2 steppers | `crates/gromos-integrators/.claude/CONTEXT.md` |
| gromos-md | L3 orchestration | `crates/gromos-md/.claude/CONTEXT.md` |
| gromos-tools | L3 system builder | `crates/gromos-tools/.claude/CONTEXT.md` |
| gromos-analysis | L4 analysis facade | `crates/gromos-analysis/.claude/CONTEXT.md` |
| pyo3-gromos | L4 Python bindings | `crates/pyo3-gromos/.claude/CONTEXT.md` |
| py-gromos | L4+ Python package | `py-gromos/.claude/CONTEXT.md` |

Each contract covers: the crate's job, what it consumes/produces, current status, key files,
and crate-specific rules that extend the workspace-wide guides.

## Prime directives

1. **The reference suite is the oracle.** A change that alters wired, tested output is a regression
   until proven a deliberate, documented decision. See `porting.md`.
2. **One core, no second codebase.** One `(Topology, State)` in gromos-core; the MD engine and the
   analysis facade both *call* it. The moment analysis re-implements a physics term, gromos++ has been
   recreated. See `architecture.md`.
3. **Shape force/state/neighbor work as composable providers over shared data** — so QM/MM/ML and
   GPU become additive, not forks. See `architecture.md`.
4. **Don't port the bugs, don't diverge by accident.** Tests before wiring; derivation before
   transcription; grep before porting; document every deliberate quirk. See `porting.md`.

## Reference sources (read-only)

| Source | Path |
|--------|------|
| gromosXX / md++ | `.local/gromosXX/md++/src` |
| gromos++ | `.local/gromosPlsPls/gromos++/src` |
| GROMOS theory book | `.local/doc/gromos_book` |
| Tutorials | `.local/gromos_tutorial_livecoms/tutorial_files` |
| Force fields | `.local/gromosXX/forcefields` |

## How to verify

`cargo build --release --bin md`, then `cargo test -p gromos-md --test test_gromosXX_references`.
See `crates/gromos-md/.claude/CONTEXT.md` for how to add a reference test.

## Global decisions

- **f64 everywhere** (not f32)
- **gromosXX `@` CLI convention:** `@topo @conf @input @fin @trc @tre @trf @trv @verb ...`
- **Tolerances:** force_abs=1e-6, energy_rel=1e-8, position_abs=1e-9
- **CLI arg parsing:** clap `#[derive(Parser)]` with `gromos_args()` pre-processor (`@key → --key`, `@f argfile` expansion). No custom arg parsers.
- **Doc style:** Rust → KaTeX + `[^label]` footnotes; Python → NumPy docstrings + `.. math::`
- **On commit:** update CHANGELOG.md and Cargo.toml version.
