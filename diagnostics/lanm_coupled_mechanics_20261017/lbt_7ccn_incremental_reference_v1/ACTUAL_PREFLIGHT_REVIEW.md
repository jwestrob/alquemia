# Independent actual preflight review — 15319

**Verdict: pass; no blocker to root binding the declared one-cell run.** Actual Slurm preflight reports `preflight_passed` under the reviewed source and configuration. The new reference is unattempted and unavailable; no incremental reference comparison is yet possible.

## Exact pins

| Artifact | SHA-256 |
| --- | --- |
| Actual preflight `RESULT.json` | `e74495a811be5066ad97515cb9a57d5dd9cf285c67282450741948421f428003` |
| `REFERENCE.py` | `bf75f502a2e3c53e51329cfd59b9ff34f661d85acc0ed24483f51868bbf194cd` |
| `CONFIG_PREFLIGHT.json` | `d6a94e65299cecf46b74b58e1faaa381d7c1bfb5607524156b780eddb3f087ff` |
| Single-cell `MANIFEST.json` | `49cb08da8e829d1021e421715abf1a77086ea2c6c014bd8fdbf5aaab3579dbf5` |
| New native input template | `b94c22fc46ff47bd8c678f2e49634948fa99195638f6f69823d0b99f075b102d` |
| New native XYZ | `ff56bcb82c9a77ce719c52bf5a4417c57721edf0f8fa41510b488735d5456c18` |
| Incremental `DIRECTIONS.json` | `cd9fe02ce05325080314cabecd45e25c38482220724d47765681c5e05e1adff0` |

Actual workspace: `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_incremental_reference_v1/preflight/execution/`. Producer job15319 ran on `node-64-768g-1` and reports 4.250164747 seconds application elapsed; scheduler accounting remains authoritative.

## Evidence checked

- All 470 unique explicit file/hash pairs traversed in the result, manifest and cached solvent receipts verify. Scientific binding equals the frozen configuration; future preflight-result, execution-decision and endpoint-result fields remain null in that preflight configuration.
- The declared comparison denominator remains `old_La,new_La`; the manifest contains exactly one analytic-gradient task, `new_La`, with one worker and 12 MPI ranks. It binds exact POINT `eb96df29`, its positions and physical IDs, 256 atoms, charge −2, La physical/effective singlet, target255 and ECP46/1030 explicit electrons. No caps or external charges are declared.
- The staged input template is byte-identical to actual old14880's input. The new output and gradient files are absent. The preflight selects no execution tasks and records zero native runner invocations, DFT, MACE, ALPB, GFN, MM, optimizer steps, new geometry or other-metal calls.
- Saved terminal/search/checkpoint and physical-guard records join exactly to completed continuation15254. Source rebasing remains false. The reported native model-work replay is exactly the stored continuation work: −1.705988079 MACE, −2.520989587 ALPB and −4.226977666 composite kcal/mol. No values were recalculated during this review.
- Both cached solvents are admitted with zero new calls: old receipt `69646f73` and new receipt `37e90528`. The executed preflight passed the reviewed raw energy/gradient, coordinate, charge, cache and native-state checks; this review inspected their saved evidence without numerical replay.
- The old native row is exactly equal to actual14880's admitted row, including raw receipt `f5bb82f7`, full gradient and residual observations. The old force report has all 256 physical IDs; the sole direction is the exact incremental chord. `new_La` remains `unavailable/not_attempted`, `full_two_available` is false and comparison is unavailable.

## Retained limitations

Old14880's native SCF admission remains valid under the unchanged protocol, but its 18-cycle stopping evidence does not pass every printed residual criterion. Only orbital gradient passes; energy, maximum/RMS density, orbital rotation and DIIS flags remain false. `ConvCheckMode=Total+1el-Energy`, `ConvForced=0`, the small-gap warning and `all_printed_residuals_pass=false` are preserved. These are not rewritten as an all-criteria convergence or stable-state claim.

The preflight qualifies inputs and reuse for one fixed new reference. It establishes neither the new endpoint's native admission nor numerical/physical force accuracy, electronic stability or reference agreement. Physical-state and force-consistency qualifications remain false. The 1 kcal/mol resolution remains descriptive; no optimizer extension, retry, extra endpoint or parameter change follows from this pass.

Review actions were saved JSON/source inspection, literal metadata equality, artifact dimensions and hashing only. No scientific imports, numerical replay, molecular calls, submissions or production changes were performed.
