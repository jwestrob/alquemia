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

## Active corrected expanded-reference evaluation

Root launched the reviewed canonical stage on14 one-thread numerical workers.
Worker1475139 (startticks41828186), completion watcher1475140. Exact commands,
manifest hash and pre-execution declaration are in
`workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/CANONICAL_LAUNCH.json`.
Frozen manifest SHA256:
`4f842a2d92364e1559dd5d266ae86e6b682a8ed81181a9e58f0a17cbd59ccfe1`.
The watcher writes `canonical_completion_event.json` and queues the real root
thread. No Slurm dependency or new optimization. Do not duplicate the worker.

Question: does the existing expanded4.3Å own-three-candidate method retain
reference fidelity after the parameter correction? Use all25 canonical references
(14Ca/11La), existing actual geometries/MACE and the unchanged calibration rule.
The canonical stage needs276 fresh xTB and276 fresh tblite inputs;24 input pairs
are reused from the completed factorial. Independent preparation/launch review
passed. Source-specific Ca/La charges span−4/−3,−3/−2,−2/−1.

Transfer is prepared but not executing. It requires a valid canonical operational
calibration under the published extrema rule (gap>0.02kcal/mol). No bands are
adjusted using transfer outcomes. Its1320 further input pairs cover3 crystals,
6 probes and101 additional pools; the complete declared100 triples use104
source/context pairs from98 structures, with six original preparation exclusions.
Coverage is primarily La-conditioned plus three Ca-conditioned probes, not the
full original two-conditioning panel. Unknown PLM labels are excluded.

All810 MACE energies are reused; neither optimization nor production changes.
The expanded own-three method remains distinct from the two-anchor common-five
sensitivity. Positive two-anchor ordering does not substitute for calibration.
Planning extrapolation for the canonical stage is~8440CPU-s (~603s ideal at14),
not a new measured runtime. Real process receipts determine actual cost.

Next precise action on completion: read
`workspaces/orca_parameter_audit_20261012/expanded_candidate_reevaluation_v3/CANONICAL_COLLECTION.json`
and its raw failures/calibration. If admitted, root executes the already-frozen
transfer stage with the same runner and a fresh transfer completion watcher.
Do not use a stale stage collection to monitor an interrupted/resumed attempt.
If a worker is externally killed, inspect its owned START.json process groups
before reallocating threads: underlying children run in separate process groups.

Curie finishes the preparation report/tests and own scoped commit; driver/PLAN
are frozen while the worker runs. Heisenberg finishes the<=1-thread geometric
q9 diagnosis, exact archived DFT inventory and vault note. No fresh DFT calls.
Root owns execution, shared checkpoint, calibrated comparisons and promotion.

Structural finding so far: CN7/denticity stay unchanged. The saved motion swaps
strengthened/weakened geometries of two included hydrogen bonds and reaches the
0.8Å displacement boundary. Omitted packing atoms have no new heavy contact
within the declared3.5Å bin; proximity alone does not establish their energy.
The exact archived q9 DFT origin uses73 atoms and different charge/context,
so it does not validate the145/192-atom saved donor motions. See the final
`q9_response_diagnosis/` report when available.

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
rootMACE and solvent allocations are released. Heisenberg uses<=1 coordinate
analysis thread. Independent guard/runtime/factorial/expanded-launch reviews
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
