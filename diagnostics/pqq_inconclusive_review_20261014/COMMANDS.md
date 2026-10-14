# Saved-artifact inconclusive review

Run from the repository root. All numerical analysis uses Slurm. These are
completed jobs, not a continuation queue; inspect their outputs before any replay.
None invokes MACE, xTB SCC, DFT, folding or optimization.

| Job | Scope | Actual receipt/result directory under `workspaces/pqq_inconclusive_review_20261014/` |
|---|---|---|
| 4037 | Score algebra, margins, all48 peers | `scores_v1/` |
| 4038 | Saved donor/PQQ geometry, confidence, optimizer flags | `geometry_v1/` |
| 4041 | Corrected contact descriptor and compact peer table | `geometry_contacts_v2/` |

Each directory retains `SUBMISSION.json` with the exact executed command,
input hash, completion-wake command and PID; `completion_event.json` contains
actual scheduler accounting. Numerical results are in `execution/RESULT.json`.
The two geometry directories also retain `execution/EXECUTION.json` with
measured process time and peak memory. Root's compact joined evidence is
[EVIDENCE.json](EVIDENCE.json).

Read the actual selected-score table without repeating calculations:

```bash
cat workspaces/pqq_inconclusive_review_20261014/scores_v1/execution/SOURCES.tsv
cat workspaces/pqq_inconclusive_review_20261014/geometry_contacts_v2/execution/PEERS.tsv
```

The frozen finite population is
`workspaces/pqq_inconclusive_review_20261014/REVIEW_INPUT.json`, SHA256
`b3ffa84989a40d02fa2371afd6d6ba26d831474bb340d2c8b9a1d00247da0a97`.
It contains nine inconclusive targets and every available same-protein peer
from the declared cohorts, deduplicated by source/context identity. It preserves
cohort aliases, the five existing inconclusive triples and old unavailable inputs.

## Scientific interpretation limits

Distances to class bands are model-score distances in kcal/mol, not error bars
or affinity estimates. Source conditioning is distinct from the protein label.
Origin scores on pooled bands are mechanistic comparisons; independently
calibrated origin decisions are retained beside them. Peer medians do not define
a new voting rule. The model-change decomposition and origin-to-accommodation
decomposition answer different questions.

Job4038's `closest_mapped_nonbonded_heavy_pairs` and `nonbonded_scope` fields
are invalid because the kinematic graph omits fixed PQQ internal bonds. Use4041's
separate inter-residue descriptor. The original output is preserved; no molecule,
energy, score, donor distance or source confidence changed in this correction.

Fresh replay requires a new output directory: the scripts deliberately refuse
to overwrite completed execution directories. Metadata-only preparation may
run on the login host; scientific imports and execution are Slurm-guarded. All
runtime source/helper/map pins must still match; do not bypass mismatches.
