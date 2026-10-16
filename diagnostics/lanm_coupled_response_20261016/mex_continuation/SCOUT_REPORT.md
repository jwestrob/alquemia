# Mex La continuation scout 12727

**The declared 12-step continuation passed.** It lowered La own energy by
3.335800 kcal/mol beyond the accepted pilot, for total work −3.588288 kcal/mol
from the original source. All 12 requested own composite evaluations were
accepted; the search stopped at its accepted-step limit.

| Incremental component | Work, kcal/mol |
|---|---:|
| Local MACE vacuum | −1.253377 |
| Embedded-minus-isolated GFN2 | +0.213806 |
| Existing MM ledger | −2.416383 |
| Fixed FF-charge global ALPB | +0.120153 |
| Total | −3.335800 |

All original-source physical guards passed. Maximum source displacement was
0.181529 Å, target displacement 0.102301 Å and maximum total torsion 0.077454 rad.
The fixed complement stayed exact; no worsened severe contacts, chirality or
peptide failures were recorded. The seven-coordinate metric retained rank 7.

The terminal represented-span gradient remains substantial: maximum atom
14.056838 and Cartesian RMS 0.194325 kcal/mol/Å. This is a bounded, nonstationary
response at EF3 with three fixed Nd spectators. Own La descent alone does not
establish the paired accommodation contrast or protein affinity. The unchanged
Dy continuation and shared common5 comparison have separate result records.

## Calls and cost

12 MACE, 24 local GFN2, 12 MM configurations/48 group queries and 12 ALPB calls;
one native model load. No new DFT or whole-protein SCC. Application wall time
119.218240 s, main-process CPU 28.415112 s, observed GFN2 child CPU 113.624913 s.
Application allocation intervals were 1430.618883 CPU-s and 119.218240 GPU-s.
Root's scheduler receipt records 121 s on 12 CPUs and one H200; these intervals
are distinct from measured process utilization. See SCOUT_ACCOUNTING.json.

## Exact saved evidence

- Result: `workspaces/lanm_coupled_response_20261016/mex_continuation_scout_v1/execution/RESULT.json`, SHA256 `3fac70c7bf1c3d5b86340ad8d041e188356c080beea762db74f9ec3f79a218ad`.
- Search: same execution directory, `La_search/SEARCH.json`, SHA256 `c969f90a69933455c1c1046b5b11b2de406e8bca3b02e0ad2972b65aa86aeff7`.
- Source/plan: `MEX_CONTINUATION.py` `c571f20e…`, `MEX_CONTINUATION_PLAN.md` `88dc7df6…`.
- Root admission and scheduler fields: `mex_continuation_scout_v1/ROOT_REVIEW.json`.

This report transcribes admitted saved results. It performs no new numerical
analysis, changes no bounds and authorizes no additional continuation.
