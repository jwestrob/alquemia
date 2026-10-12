# Corrected Nikasha parameter/context work — 12 October 2026

## Result and current decision

The live Linux native-ORCA parameter correction is reproduced. Sixty actual
controls are terminal:52 admitted, eight original/default-roundtrip density
failures preserved. Corrected La and14 corrected transition-metal aquo cells
match the archived independent controls. Zn/water population parity remains
unqualified. Read runtime/REPORT.md and the independent runtime_review report.

The strict208 corrected calibration remains unsupported (gap−4.556602kcal/mol).
The complete prior2796-point correction was imported, not rerun; see
recovered_strict208/IMPORT.json. Corrected origin-only203/3/2 and frozen geometry
208/208 remain separately scoped results over correlated folds from25proteins.
No threshold tuning or production promotion. Current operating guide and PLM SOP
now identify the earlier native energetic performance claims as historical.

## Completed context experiment

All40 composite endpoints are now available. The independent review agrees with
all energies, selections, interactions and costs. Read context_comparison/REPORT.md
and RESULT.json. Root acknowledged the actual context-factorial-solvent-v2 wake.
No factorial workers remain; do not rerun completed80input pairs.

Two-anchor ordering (La reference c5 minus Ca reference q9) in original/expanded
contexts, kcal/mol:

- Origin: +8.867782 / +12.175484.
- Original three proposals: −4.556602 / −1.212842.
- Expanded three proposals: +0.437421 / +3.472171.
- Common five: −3.714195 / +1.265264.

Expanded context alone does not restore the original proposal pool's ordering.
q9's original La proposal still gives~12.1kcal/mol La-favoring motion response
in both contexts. Context adds a mostly static offset; geometry selection also
matters. A positive two-anchor ordering is not a full canonical calibration.
No production promotion or claim that accommodation is repaired.

MACE:24 reused energies+16 new cross-context+four exact CPU/GPU replays. Actual
20calls used113.884s elapsed/590.696 processCPU-s/4.25GiB peakRSS, zeroGPU.
GFN2:80primary+80independent+four causal diagnostic calls,2490.296 measuredCPU-s.
All40 solvent components pass v2; every original absolute-energy comparison fails
v1 and stays recorded. Density-independent nonsc dispersion cancels exactly
within same-coordinate water−vacuum terms; no fitted offset or changed parameters.
Reports: context_MACE/, context_GFN2/, context_comparison_review/.

## Canonical calibration passed; structural transfer is running

All25 canonical pools are available (14Ca/11La), with all150 solvent components
passing the frozen component/population gates. The unchanged extrema rule yields
an accommodated class gap of3.472171047kcal/mol (Ca edgeq9z4j7, La edgec5axv8).
Origin-only also separates, gap10.297678912 (La edgeA0A3F2YLY8). These are
calibration results on consumed training inputs, not independent predictive
accuracy or evidence that accommodation is better. Original strict208 calibration
and the two-anchor common-five sensitivity remain separately reported.

Actual canonical cost:552 fresh calls (276xTB+276tblite),10085.975583CPU-s,
744.688672s molecular-stage elapsed at14workers, peak single-processRSS127052KiB.
24 input pairs/48 backend calls were reused and are not counted as new compute.
Maximum component residual1.29869250e-7Eh; population residual1.09345701e-7e.
No task exceptions or failed canonical components. The original absolute-energy
parity gate remains failed; no raw absolute totals are made equivalent.
Root inspected raw receipts and acknowledged the actual completion wake.
Independent raw reconstruction also passes: see `expanded_canonical_review/REPORT.md`
and RESULT.json (commit24401b2), including all300 input records and preserved
absolute-energy failures.

The already-frozen transfer stage now runs on14 one-thread numerical workers:
worker1490237 (startticks41917155), completion watcher1490238. Exact commands,
input/collection hashes and pre-execution declaration are in
`workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/TRANSFER_LAUNCH.json`.
Frozen manifest SHA256:
`4f842a2d92364e1559dd5d266ae86e6b682a8ed81181a9e58f0a17cbd59ccfe1`.
Frozen canonical collection SHA256:
`929dedab7ba46c7088de7a4a695a8f19fc5af60066ecd3599645f687bd4e6896`.
Do not duplicate the completed canonical stage or active transfer worker.

Question: do those unchanged canonical-only bands retain reference fidelity on
structural transfer? The1320 further input pairs cover3 crystals,6 probes and
101 additional pools; the complete declared100 triples use104source/context
pairs from98 structures, with six original preparation exclusions. Both pooled
and origin-only results retain their own fixed calibration and availability.
All810 MACE energies across the full135pools are reused; no optimization,
new structures, state changes, threshold tuning on test outcomes or production
promotion. Charges remain source-specific (Ca/La−4/−3,−3/−2,−2/−1).

Coverage is primarily La-conditioned plus three Ca-conditioned probes, not the
full original two-conditioning panel. Unknown PLM labels are excluded. Curie
inventories additional real Ca-conditioned preparation/proposal reuse, read-only;
no extension is running or implied. The expanded own-three method stays distinct
from the two-anchor common-five sensitivity. A larger training gap alone is not
an accuracy benefit.

Next precise action on completion: read
`workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/FINAL_COLLECTION.json`,
all unavailable cells and calibrated paired/triple comparisons, then acknowledge
`transfer_completion_event.json`. Its watcher queues the real root thread.
Do not use a stale stage collection to monitor an interrupted/resumed attempt.
If a worker is externally killed, inspect its owned START.json process groups
before reallocating threads: underlying children use separate process groups.

Root owns14 numerical workers. Independent canonical review is finished; Curie's Ca-conditioned inventory
owns<=1 read-only thread. Total<=16. All previous
runtime/factorial/canonical workers are terminal. Driver/PLAN remain frozen.

The q9 diagnosis is complete in `q9_response_diagnosis/REPORT.md` (ed2cf8f).
CN7/denticity stay unchanged; two included hydrogen-bond geometries change in
opposite directions; both saved La proposals reach0.8Å displacement boundary.
Omitted packing atoms have no new heavy contact within the fixed3.5Å bin.
The archived73-atom DFT origin has different charge/context and does not validate
the145/192-atom donor motions. Zero new molecular calls in that diagnosis.

Exact read-only completed comparison reproduction (choose a fresh output path):

```bash
python scripts/compare_context_factorial.py \
  --mace diagnostics/orca_parameter_audit_20261012/context_MACE/RESULT.json \
  --solvent workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json \
  --output /tmp/nikasha_context_factorial_recheck_20261012
```

## Ownership and safeguards

Jacob cleared the pause and explicitly authorizes<=16 TOTAL numerical login
threads. Recorded verbatim atop AGENTS.md. Root's current worker owns14;
all remaining light preparation/analysis shares at most2. Earlier runtime,
rootMACE and solvent allocations are released. Canonical review is complete; Ca-conditioned inventory uses<=1 thread. Independent guard/runtime/factorial/expanded-launch reviews
are complete. New guarded dispatch is committed/pushed asdd758db; its molecular branch
remains explicitly unexecuted. Root owns scoring/reporting, current-facing docs and integration.

Parameter guard commitcd97a66 implements exact task-local full-file imports,
shell-aware audits, runtime/profile/state admission and pre/postflight hashes.
Historical original/default controls remain diagnostic-only. The guard's Ca/La
registry uses actual current receipts; it does not globally patch old production
executors or qualify new-metal mixed profiles automatically. The new guarded
dispatch wrapper is separately tested; inspect its report before use.

Two headless Claude attempts produced zero science: first OAuth refresh failed;
one technical retry refused the explicitly authorized login scope. Their receipts
remain under workspaces/claude_parameter_audit_20261012/. Codex completed runtime
work without changing authentication or asking Jacob to repeat authorization.

No Slurm science jobs were submitted. The30minute readiness monitor remains
activePID669874; latest probe119 still fails account/partition admission. On recovery
inspect actual live jobs; do not duplicate login work or automatically restart
old LanM campaigns. Historical outputs, existing production and concurrent files
are preserved. Scoped commits/pushes remain authorized.
