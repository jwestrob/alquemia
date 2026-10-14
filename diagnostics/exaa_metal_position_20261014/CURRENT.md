# Reciprocal placement supported in vacuum; solvent4132 running

Preparation4129 and MACE4130 completed8/8 new cells plus8 exact reuses.
Control delta_D +3.373337 and target −9.328603 kcal/mol pass the frozen reciprocal
rule. Moving only the ion changes the saved donor response in both predicted
directions. Read REPORT.md/RESULT.json; no optimum/classifier claim.

Solvent4132 runs eight shifted geometry/state cells on7 shared standardCPUs,
node-64-768g-5. Actual internal preflight passed. Exact8 originals were verified;
all16 MACE values are reused. Maximum32 SCC plus16 native ALPB component calls.
Uses qualified standaloneGFN2 component-v2 and native frozen-charge recipe;
historical absolute parity failures stay failed. No newMACE/DFT/optimization.

Completion watchPID1093803 is armed, event exaa-reciprocal-metal-solvent-v1.
Exact submission/rearm: workspaces/exaa_metal_position_20261014/solvent_v1/
SUBMISSION.json. ConfigSHA0f162cb433686011436a4c9f04a71ebded1e6c3e8de112383032e7ec6933d708.
Read execution/FINAL_COLLECTION.json and EXECUTION.json when terminal; keep
unavailable cells. Apply SOLVENT_PLAN.md's frozen primary reciprocal rule, and
report its separate two-vs-four common geometry minima without class bands.

Root owns execution/integration. Independent inputaudit and staticreview passed.
Preparation/scoring/review agents are complete; no duplicate jobs or automatic
optimization continuation. Production/defaults/calibration remain unchanged.
Earlier MACE4130 used118s/7CPUs (826allocatedCPU-s), eight model forwards; all
previous jobs terminal. Delayed old callbacks require acknowledgement only.
