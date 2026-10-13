# Ca reference coverage recovery — 12 October 2026

Jacob authorized redo of the eighteen unavailable Ca-conditioned cases. This is
separate recovery evidence; the completed primary125 collection and frozen25
canonical calibration remain unchanged. No production promotion.

## Actual findings and execution

- Job2866 completed the real15-source PQQ preparation audit. Every source has
  intact PQQ; all selected Ca atoms lie25.0–25.2 Å away at another six-O site.
  Normalization replay and heavy coordinates are exact. No metadata repair.
- Job2864 completed the MMOL_2048 Ca0 disulfide diagnosis. Cys127–128 S–S is
  intact (2.118 Å); existing Cys128 HB2 is1.377 Å fromCB and loses its parent
  under the source-graph cutoff. Repair2870 passed: only that H moved0.287 Å to the ff19SB1.09 Å
  C-H length; source heavy atoms/core/inventory are unchanged.199-atom envelope
  and frozen state pass; two untouched controls replay. Native scoring preparation
  follows; no repaired score exists yet.
- MxaF P12293 Ca1/Ca4 recovery job2869 runs one warm7CPU worker onnode-64-768g-5.
  Exact maps/selector preflight passed. Reuse eight native cells; restart only
  two interrupted La searches from original q0, plus required cross cells.
  Watcher3490926 is armed; manifest/submission/watch under mxaf_v2 workspace.
  Both native startup energy/force controls passed exactly within tolerance.
  Final result does not yet exist. Solvent recovery is prepared next and reuses
  all compatible archived components, with the original frozen bands.

## Targeted structure generation

ExaA crystal1FLG and its primary paper explicitly show a second N-terminal Ca
stability site besides PQQ Ca (https://www.rcsb.org/structure/1FLG;
https://doi.org/10.1006/jmbi.2000.3603). This independently motivates a two-Ca
source hypothesis. Agent Heisenberg prepares exact mapping/input/runtime for
one finite AF3 batch each for ExaA, PedE and ExaF. PedE/ExaF occupancy is a
homology-based hypothesis, not experimentally established by this campaign.
Preserve original protein/MSAs/templates/PQQ and seed1; add one Ca, five samples
per protein,10recycles. Keep all outputs. Select the unique PQQ-associated metal
by the existing structural gate, never by predicted class. No fold launched yet.
Current GPU QoS caps one GPU and12GiB host memory/job; use admitted resources
without changing QoS or artificial memory/time request. Runtime scout precedes
expansion; missing inference remains unavailable.

## Ownership / restart

Root owns jobs/scoring/integration; Curie owns MxaF recovery+solvent adapter,
Gauss disulfide preparation, Heisenberg refold preparation. All science uses
Slurm. No readiness probes; individual completion wakes stay armed. Existing
LanM2665/2666 and collector2669 remain preserved. Stale report allocations2832/
2833 on mismatchednode-48-256g-18 are unrelated and must not be copied/rerun.
Read actual SUBMISSION.json, WATCH.json and terminal artifacts before any action.

Next: read workspaces/ca_coverage_recovery_20261012/mxaf_v2/FINAL_COLLECTION.json
when present; launch exact two-row solvent recovery only after native admission.
Submit separately pinned H-repair and fold-input preparation jobs when ready.
All original failures, costs and outcomes remain visible alongside recovery.
