# gromosXX's own regression oracle (`md++/src/check/`)

Input files of gromosXX's `check` suite (`md++/src/check/data/`, GPL, same origin as every other
reference here). `check/aladip.t.cc`, `aladip_unperturbed.t.cc` and `aladip_special.t.cc` hard-code
per-term energies for the alanine dipeptide in 20 SPC waters — `aladip.topo`, `aladip.pttopo` and
`aladip.conf` in the parent directory are byte-identical to the suite's copies. The values are
evaluated by `crates/gromos-md/tests/gromos_check_suite.rs`, term by term, independent of the `md`
binary's trajectories: a second source of truth for the bonded, perturbed, soft-core and
reaction-field terms (PLAN.md "Cross-cutting — reference tests").

`aladip_distres.in` is `aladip_special.in` with the ANGLERES, DIHEDRALRES, XRAYRES, LOCALELEV and
ORDERPARAMRES blocks removed — the restraint kinds gromos-rs does not model yet (PLAN.md 1.6);
the distance-restraint values are evaluated with it.
