# Midpoint prepared; common-seven scoring submitted

Preparation4136 completed: both real ExaA hosts and all eight new midpoint cells
pass unchanged physical gates. No donor atom, chemistry or membership changed.
Preparation took3scheduler seconds on1CPU; no energy calls.

Root submitted job4137 on seven shared CPUs for the finite MACE/qualified frozen-
charge solvent sequence. Config SHA256:
`f0edcf081a019a5b5bffb5278050fb1ce6052361837064527a339d213b2c7147`.
Completion watcherPID1131772 is armed for event `exaa-midpoint-common7-v1` in the
current root thread. Exact submission/rearm argv is in
`workspaces/exaa_metal_midpoint_20261014/scoring_v1/SUBMISSION.json`.
Inspect actual state before acting; do not repeat completed preparation or old20
cells. Expected result: `scoring_v1/execution/FINAL_COLLECTION.json` in that workspace.

The declared comparison is historical3/previous5/new7 common geometry pools,
with all20 old cells reused. The previous full transfer gave exactly zero selected
score benefit; it remains a negative result. New scoring is not yet interpreted.
This midpoint uses a force-motivated fixed fraction0.5 and both saved donor
arrangements. No optimizer, DFT, thresholds or production change.

After collection, report actual per-metal and paired-contrast changes and costs.
No automatic finer fractions, bisection or joint search follow this experiment.
Root owns execution/integration; three disjoint preparation/scoring/review agents
have completed implementation. See PLAN.md and EXECUTION.json for scope and pins.
