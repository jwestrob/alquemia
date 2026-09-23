# Completed parallel pilots — 2026-09-22

Jacob approved the three finite experiments with “I approve.” All jobs and
collections are terminal. Do not relaunch them. [REPORT.md](REPORT.md) and
[RESULT.json](RESULT.json) are the final integrated result; PLAN.md/INPUTS.json
preserve exact scope and sources. Production and the larger goal remain unchanged.

## Actual findings

- Alternative starts: four of16 fit the existing physical domain, producing eight
  searches, all returning the old native basins. The four difficult folds admit
  neither template. No added discrimination.
- Solvent probes:27admitted points, all eight scores available, unchanged7correct/
  1wrong under transferred adaptive bands. No improved selected-source spread.
- Basin grid: all512MACE calls succeed;1022/1024GFN2 converge. Two1H4I La vacuum
  failures remain missing. Two of four curves pass frozen integration gates;
  width contributions are−0.107/−0.084modelkcal/mol. Q9Z4J7 fails quadrature/domain
  checks. The expanded pool retains old adaptive minima on all three available
  sources. No demonstrated classifier benefit or routine entropy correction.
- Q88JH5 exact repeats reproduce a4.805kcal La solvent-transfer difference between
  geometries only6.46e-7Å apart. Large finite-step solvent jumps also occur on
  some basin curves. A unique numerical cause/remedy has not been established;
  passing quadrature is not validation of the electronic surface.

Final actual common collection for the basin is
`workspaces/local_basin_breadth_20260922/scoring_v1/after_solvent_3_1210125.json`.
It was reused directly, with no duplicate collector. Integral result is
`workspaces/local_basin_breadth_20260922/RESULT_v2.json`; a partialv1 remains
archived after a JSON boolean serialization failure. Only algebra was replayed.
Root comparison is `workspaces/nikasha_parallel_pilots_20260922/COMPARISON_v1.json`.
The common-pool minimum and the integral's own65-node minimum are distinct;
no width contribution is added to an unrelated adaptive winner.

## Execution and ownership

Second_shell owns the completed starts/numerical-repeat report; Khoury owns the
completed solvent probes/native-control note; water_basins owns the completed
integral and figures. Root integrates. Existing scientific inputs/recipes and
frozen acceptance criteria were retained; failures were not rescued.

Finished jobs:1210098/1210103/1210104/1210115 (starts and exact repeats),
1210100/1210101 (solvent probes),1210107 (basin MACE),1210122–1210125 (basin GFN2).
Four original basin allocations1210108–1210111 were canceled only after pending/
zero-runtime checks and replaced with identical manifests directly on gpu;
inherited whole-node exclusivity had prevented CPU sharing. No chemistry duplicated.
Two replacement jobs have FAILED scheduler status from the explicit SCF failures;
their successful cells and final collection remain available. No active watcher.

Total861MACE/1172GFN2 attempts (1170converged), zeroDFT;
265472allocated core-seconds and380requested GPU-seconds. These include allocation
and collection overhead; local setup/tests/reporting add unmetered work.
No production latency or new calibrated accuracy estimate follows from this pilot.

## Remaining proposals, not queued

Close the tested bounded-start/probe rules and omit routine width integration.
Native.xtbw self/cross-restart (eight proposed calls), standalone-xTB fallback,
collective scaffold response and PQQ-redox work remain proposed-only. No new
backend, chemical state, production rescore, default change or remote push occurred.
Recover their exact plans and current user agreement before execution. The earlier
adaptive improvement remains recorded separately in nikasha_next_phase_20260922;
this phase does not establish another accuracy gain or close the larger goal.
