# Approved activation recovery — same eight logical restart tests

Root explicitly approved this technical correction after the eight seed-only
attempts completed as SAD calculations. Preserve job1210179 and all its costs,
outputs and unavailable restart fields. The original plans remain immutable.
The broader user authorization stays pinned from the shared round PLAN.

The initial-guess manual specifies a same-basename GBW as the AutoStart trigger;
the native manual specifies xtbw precedence once restarting. The initial no-GBW
preparation omitted that documented trigger and produced SAD in all eight calls.
The recovery adds each seed's **matching, already archived GBW**, alongside its
xtbw, at `endpoint.runtime.gbw` and `endpoint.runtime.xtbw`. No extra input keyword
is introduced. Both source files come from the exact same pinned prior calculation.
[AutoStart documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/initialguess.html#autostart-feature),
[native restart documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/modelchemistries/semiempirical.html#native-gfn-xtb-and-gfn2-xtb).

Exactly eight recovery attempts: old/new destination ×vacuum/ALPB ×self/opposite
seed, with the same source XYZ/state/parameters/MaxIter500/300K/native mixer as
before. Copy immutable pre-run versions of both seed files and retain post-run
GBW/xtbw and any AutoStart GES artifact. Verify source file identity, absence of
other restart files and actual output. Only XTBRESTART plus native mixer and all
existing electronic-state/parameter/charge/receipt checks qualifies a native
restart; orbital or SAD fallback remains unavailable. No ordinary-SCF, new state,
solution selection, standalone backend or additional logical comparison.

Use a fresh finite manifest and the same CPU-only gpu allocation:64CPU/128GiB,
eight8-rank calls. Count all16actual attempts across original and recovery. No
further call follows automatically if native restart remains unconfirmed.
