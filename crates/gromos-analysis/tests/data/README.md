# gromos-analysis reference outputs (gromos++, `.local/gromosPlsPls`, 2026-08-30)

`shared/` is `crates/gromos-md/tests/gromosXX_references/shared/`, `aladip.trc` the trajectory
`gromosXX_references/aladip_solvated/expected/trajectory.trc`. gromos++'s gather/Boltzmann warnings
are stripped from the stored outputs.

| dir | gromos++ command |
|-----|------------------|
| tser | `tser @topo shared/aladip.topo @pbc r @prop "d%1:1,2 a%1:1,2,3 t%1:1,2,3,4 tp%1:2,3,4,5%60" @traj aladip.trc @dist 5` |
| mdf | `mdf @topo shared/aladip.topo @pbc r @centre 1:1 1:6 @with 1:3-12 @traj aladip.trc` (writes `MIN_1:1.dat`, `MIN_1:6.dat`) |
| dg_ener | `dg_ener @temp 300 @stateA eA.dat @stateB eB.dat` and `… @col 3 3` (synthetic series, `eA`/`eB` from a seeded generator) |
| dfmult | `dfmult @temp 300 @stateR eR.dat @endstates eX1.dat eX2.dat` |
| matrix_overlap | `matrix_overlap @m1 m1.dat @m2 m2.dat @dim 3` |

`tests/gromospp_references.rs` runs ours with the same arguments and compares every number
(1e-6 relative; `-nan` matched literally) and every label.
