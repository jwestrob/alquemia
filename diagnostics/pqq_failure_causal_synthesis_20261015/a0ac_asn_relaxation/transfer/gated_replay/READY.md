# Retrospective gated replay — executed as job6045

Completed15October2026 and verified after the reboot on16October. Read
[REPORT.md](REPORT.md) and [RESULT.json](RESULT.json). The command below is the
executed command, retained for provenance; its output directory now exists.
A deliberate new replay needs a new output directory. No rerun is needed.

## Original pre-execution instructions (historical)

The [frozen plan](PLAN.md) declares the post-inspection rule. Gauss independently cleared the source/schema review; no imports, numerical tests, aggregation or submissions occurred in preparation/review.

Root submits one normal Slurm CPU from the repository root:

```bash
sbatch diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/gated_replay/RUN.sbatch \
  "$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/gated_replay/REPLAY.py" \
  "$PWD/diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/transfer/gated_replay/PINNED_INPUTS.json" \
  789f2537c5670d62681e7d45fd12712581994282e0707fb52ab0caeddbc79322 \
  "$PWD/workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_gated_replay_v1/execution"
```

No workspace exists from this task. The output directory must be new. The explicit script path avoids Slurm spool-directory resolution. Only pinned stdlib Python `/usr/bin/python3.12` is used (invoked by `/usr/bin/python3`), SHA `e50d468e8b0adfb05733f5b87b3cff34829c4a8c1aea50c865aa8bdfe4bb150f`.

## Exact pins

| Artifact | SHA256 |
|---|---|
| `REPLAY.py` | `f178da838df8ab358b471d773f90bdf9bc60f996043960b557f9be77de756f34` |
| `PLAN.md` | `505488bbc8b33613ac4216c635c39a08582a6b6e71e5a8c743a843d9b54870ad` |
| `RUN.sbatch` | `9bde8f8c287c8c534cf4e542fe0ba9b961c3cc573e2bc18798db14d8c8255d87` |
| `PINNED_INPUTS.json` | `789f2537c5670d62681e7d45fd12712581994282e0707fb52ab0caeddbc79322` |
| Completed decision `RESULT.json` | `f70378b775db2844a18f1cbfb5a974b3b3eb6b930aab2b0ea65a0e6831ba61eb` |
| Completed decision `INPUTS.json` | `390c6b34e4a66d7eea8c50252228232c36b9b09f34b6559e9bdf91c3419ab0e6` |
| Completed decision `EXECUTION.json` | `aa4200a0df8136cf94d2f6dda373873c165676ab8b47154d553fa07bab96044c` |

The decision files live in `workspaces/pqq_failure_causal_synthesis_20261015/asn_transfer_decision_v1/execution/`. Their exact reporter files and executed `decision/SUMMARY.py` SHA `8f40e005bb5e59f1e219455af6dde4e92310b57960d9b01c4aef67089634e4a5` are reused. Every nested input pin is checked by the unchanged saved-summary replay, which must exactly reproduce the completed decision result. `CANDIDATES.tsv` is separately pinned for explicit failed endpoint evidence. There is no scoring, classification or scientific-kernel import.

## Outputs and interpretation

`RESULT.json` retains all 257 routed contexts and 260 memberships, alias identities, both policy comparisons, both selected geometries, raw saved scores, any score/selection/decision difference, regression, unavailable screen or required result, and all untriggered differences. It copies the unchanged canonical bands, 18 superseded records and saved failures. `INPUTS.json` and `EXECUTION.json` retain exact provenance, placement, duration and any failure traceback.

The trigger is **only** the existing common4 operational selection of `amide_flip` by either metal. Both reporting policies follow that same route: common6 if triggered, common3 otherwise. Unknown screens and missing required pools remain unavailable. New class decisions are not computed; the routed saved pool's existing decision is copied.

Two search slots per triggered or untriggered context quantify replay eligibility or avoidance across the declared population, including contexts whose searches were historically reused. They are not measured runtime/cost savings or counts of new jobs. All triggered cases still require both metal searches and their cross-scores in any future implementation. The screen itself also has a cost. An unfavorable fixed flip can relax usefully; agreement on this consumed archive cannot establish unseen-basin recall or justify automatic production promotion.
