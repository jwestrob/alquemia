# Existing chemistry resolves all eight abstentions with a canonical-only midpoint

**A threshold fixed from the designated canonical25 correctly classifies all257
examined source contexts.** No score, geometry, energy method, source selection or
label changed. The eight previously inconclusive calls resolve correctly. This
supports changing the decision policy as the immediate development direction;
more molecular adjustments are not necessary merely to resolve these eight calls.

## Actual comparison

| Decision rule | Correct | Wrong | Inconclusive | Unavailable |
|---|---:|---:|---:|---:|
| Existing canonical bands | 249 | 0 | 8 | 0 |
| Canonical midpoint, primary | 257 | 0 | 0 | 0 |
| Calibration-sensitivity band, secondary | 255 | 0 | 2 | 0 |

Both operational and mathematical selection give these outcomes. All25 canonical
contexts remain correct. Outside those contexts, the midpoint gives232correct;
the sensitivity band gives230correct/2inconclusive. These structures share proteins
and development history with calibration; this is not a232-protein blind test.

All94 available pre-existing three-source summaries now classify correctly under
both new rules, compared with89correct/5inconclusive previously. The six original
preparation exclusions remain unavailable. The100 overlapping triples use the
same members and median-score rule, not a new vote or favorable subset.

## Why the result changed

The old rule requires a score to reach the canonical Ca or La class range before
assigning a class. Scores between those ranges remain inconclusive. The primary
new rule instead places one boundary halfway across that canonical gap:

`T = Ca_max + (La_min - Ca_max)/2 = -405457.73794054205 kcal/mol`

Ca is below T, La above T, and an exact tie is inconclusive. This raw contrast is
the existing model score, not an affinity or binding free energy. Its large offset
is unchanged. Thresholds were computed and written before transfer outcomes were
loaded. The canonical extrema are PedE Q88JH5 and C5AXV8; evaluation extrema did
not enter fitting, and no threshold was adjusted after this comparison.

The secondary rule spans all midpoints obtained by removing one calibration
biological group at a time, including the full midpoint. Its bounds are
`[-405458.3097447878, -405454.31015960994] kcal/mol`. This is a calibration-leverage
diagnostic, not a statistical confidence or prediction interval. It retains two
abstentions:

- ExaA Q9Z4J7, La-conditioned sample3; expected Ca.
- MMOL1770, Ca-conditioned sample0; expected La.

Thus all current class calls can be correct under the primary rule while these
two calls remain sensitive to which calibration proteins are included. The
midpoint result does not establish that source dependence or uncertainty vanished.

Separately, leaving each calibration protein out and predicting that protein
with its own omitted-group threshold gives25/25 correct. Full-envelope predictions
were not presented as held-out. Historically consumed model-development examples
still do not become a new blind validation set through this calculation.

## Interpretation and recommendation

Proceed with the canonical midpoint as the named candidate decision rule and
retain the calibration-sensitivity result as an explicit warning on borderline
calls. This experiment shows useful signal was withheld by the old bands on the
examined benchmark. It does not add new physical accuracy, establish broad-family
transfer, supply probabilities, or justify fitting thresholds to future PLM outputs.
Do not add more physics merely to force these eight archived cases across old
bands. Production defaults and released calibrations have not been changed.

## Execution, checks and reproducibility

Slurm11426 completed normally in3seconds on1CPU (3allocatedCPU-seconds;
2.296996226s application). Zero MACE/GFN2/DFT calls, no GPU allocation, and no
memory/time override. Exact archived-band decisions, alias joins, canonical scores,
triple medians and exclusions reproduced. Four real-calibration algebra checks
passed. Independent review verified source/input/result pins and the interpretation.

Protocol: `pqq_calibration_only_decision_v1`. Immutable attempt:
`workspaces/pqq_decision_policy_20261016/comparison_v1/`.
Input SHA256: `6ccd51a4a2ac48183241ca2d0c914b98d8676752f4b48a3f886071daf3895a11`.
Full result SHA256: `adf8cce4e41a8e5fae30a05e826d1e745b9d489cc8ac688174018fc7fd8e6684`.

Read [RESULT.json](RESULT.json), [SUMMARY.tsv](SUMMARY.tsv),
[every changed call](TRANSITIONS.tsv), and
[held-out canonical records](HELD_OUT_CANONICAL.tsv). Full per-context and triple
exports, frozen thresholds, code snapshot and receipts remain in that workspace.
The [plan](PLAN.md) and [implementation](COMPARE.py) specify every rule.

The completed job should not be resubmitted. To reproduce this saved-data analysis
in a separate output directory, run from the repository root (one Slurm CPU):

```bash
sbatch --parsable \
  --output=workspaces/pqq_decision_policy_20261016/comparison_v1/replay-%j.out \
  --error=workspaces/pqq_decision_policy_20261016/comparison_v1/replay-%j.err \
  workspaces/pqq_decision_policy_20261016/comparison_v1/RUN.sbatch \
  workspaces/pqq_decision_policy_20261016/comparison_v1/COMPARE.py \
  workspaces/pqq_decision_policy_20261016/comparison_v1/INPUTS.json \
  6ccd51a4a2ac48183241ca2d0c914b98d8676752f4b48a3f886071daf3895a11 \
  workspaces/pqq_decision_policy_20261016/comparison_v1/replay_execution
```

This optional reproduction refuses an existing output directory. Root's original
submission and completion watcher are recorded in SUBMISSION.json. No reproduction
or additional experiment is needed to obtain the reported result.

## User decision — 16 October 2026

Jacob selected the canonical midpoint with the sensitivity flags retained.
See [the preserved decision](DECISION.md) and [exact policy](SELECTED_POLICY.json).
The full scalar context/triple exports and frozen calibration are now committed
alongside this report. PLM manuscript inclusion is deferred; research returns to
La/Dy profiling. Historical production commands and results remain unchanged.
