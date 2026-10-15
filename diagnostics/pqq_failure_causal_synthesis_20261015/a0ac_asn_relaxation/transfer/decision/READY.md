# Saved transfer decision summary

Source preparation only. No aggregation, numerical test, submission or molecular call has run. The exact completed reporter files are pinned in [PINNED_INPUTS.json](PINNED_INPUTS.json), SHA256 `a8e548187ee88daeaf4d77702d5e0e45348fe51f02f1c402c037bf98ea50a6fc`.

From the repository root, root submits the single stdlib command on one normal Slurm CPU:

```bash
sbatch diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/decision/RUN.sbatch \
  "$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/decision/SUMMARY.py" \
  "$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/decision/PINNED_INPUTS.json" \
  a8e548187ee88daeaf4d77702d5e0e45348fe51f02f1c402c037bf98ea50a6fc \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_decision_v1/execution"
```

The explicit script path avoids the Slurm spool-directory problem. The output directory must be new. Python `/usr/bin/python3` resolves to `/usr/bin/python3.12`, SHA256 `e50d468e8b0adfb05733f5b87b3cff34829c4a8c1aea50c865aa8bdfe4bb150f`; no third-party dependencies. Script SHA256 `8f40e005bb5e59f1e219455af6dde4e92310b57960d9b01c4aef67089634e4a5`; wrapper `359b70fa7b6002fc954c4dbdb863e049b50499aba99caae779a16060b7acc35b`.

## What it reports

- Unique-context outcome counts for all three pools under operational and mathematical policies, separately from the recorded cohort membership counts. All 257 contexts / 260 memberships remain. Context aliases must agree in identity, scores, selections, decisions and every exported component/transition before deduplication.
- Exact identities and saved values for changed decisions, scores and selections across common3→4, common4→6 and common3→6; newly correct contexts, regressions and lost coverage remain explicit. `marginal_relaxation` separates fixed-flip improvements from subsequent relaxation and compares the newly correct context sets without assuming they differ.
- Every remaining common6 inconclusive, wrong or unavailable context, with all membership IDs. Fixed canonical bands are copied and checked against the recorded original bands; no decision is recalculated.
- Saved per-protein/conditioning/cohort ranges and source identities from `SOURCE_SPREADS.tsv`. Range differences are descriptive and only computed when available membership is identical. They are not pooled across overlapping cohorts.
- Existing 100-triple outcome counts, changed medians/decisions and all unavailable triples, retaining original definitions. Triples and related folds are correlated observations.
- All 18 superseded records separately, preserving their failures without inflating the 257-current-context denominator. Exact exported stage costs and resource links are copied without a new campaign total or relabeling allocation time as measured utilization.

Outputs: `RESULT.json`, exact `INPUTS.json`, and `EXECUTION.json` with Slurm placement, implementation hash, duration and failure traceback if needed. The source report remains authoritative for scientific admission. This summary does not recompute endpoint scores, select geometries, fit bands, optimize, or make new molecular/classification calls. Partial scientific coverage is reportable when the saved reporter itself completed.

Heisenberg independently cleared the source/actual-schema review at the hashes above; no concrete blocker was found. No imports, tests or numerical execution were performed in review. Runtime aggregation remains unrun; root owns submission and result interpretation.
