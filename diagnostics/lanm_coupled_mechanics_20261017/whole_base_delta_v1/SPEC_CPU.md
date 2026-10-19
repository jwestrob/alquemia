# Development CPU continuation: four calls, exact scout reuse

The original root design SPEC, eight-state inventory, preflight/scout source and
their actual results remain immutable. This separate adapter implements only
the next CPU allocation. It cannot load MACE, run solvent/MM/DFT, change geometry,
or begin the 7CCN challenge.

Require actual passed zero-call preflight14799 and admitted one-call scout14800,
including native wall time ≤600seconds. Rerun the saved-only preflight admission
and exact eight-state signatures using frozen BASE_DELTA.py. Re-admit the scout's
raw outputs with the original native criteria and compare its exact task,
physical IDs/XYZ and recorded cost. Reuse that result without evaluation.

New evaluations, sequential and fixed before execution:

1. 7CCO after-La, whole256 vacuum GFN2 analytic gradient.
2. 7CCO before-Dy, whole256 vacuum GFN2 analytic gradient.
3. 7CCO after-Dy, whole256 vacuum GFN2 analytic gradient.
4. 7CCO after-Dy, existing baseline172 capped region, isolated GFN2 gradient.

The fourth request uses the exact archived after-La map/geometry and changes
only the already-declared target element to Dy, effective UHF0/charge−2.
Physical multiplicity6 remains separate metadata. There is no external field
in either whole or regional requests. Full750/regional538 valence counts and
physical spin parity are supplied by the passed common preflight.

Use the pinned qualified corrected profile, complete parameter file, executable,
libraries/environment and original native admission. GFN2/300K/accuracy0.0001/
500iterations/no restart remain unchanged. Each fresh native directory stages
the complete parameter bytes and native outputs. No damping variant, fallback,
repeated default failure or retry. Stop the loop on the first native/admission
failure; preserve all eight rows and every attempted raw record. Missing later
calls remain explicitly not requested. Molecular counters are incremented
before the native call, including failed calls.

Successful completion requires exactly four accepted new calls: three whole
and one regional. Including the unchanged scout, cumulative calls are five
(four whole, one regional); MACE remains zero. Only after this actual admission
may a separately reviewed one-GPU stage perform the single missing regional
MACE calculation. The candidate-wide ten-call ceiling and seven cached regional
pairs are unchanged. The future challenge is gated by the original both-metal
1kcal/mol development triage; no automatic scheduling is implemented here.

All execution uses Slurm, one CPU, native threads1, no GPU and scheduler default
time/memory. Every invocation uses a fresh output directory. `CPU.py` imports
scientific helpers only after its Slurm guard; source checks on login are AST
and metadata only. The actual scripts/configs/results are hash-linked. Parent
SCC evidence, native wall/CPU/RSS receipts and allocation identity remain
available; numerical derivative/state/accuracy qualification remains false.

The result exposes ordered eight `rows` with whole-GFN receipt pins and the one
new regional-GFN receipt, plus the saved-only eight-state manifest. Separate
GPU/assembly source supplies the candidate/base scalar and force records for
the frozen saved comparator. No survivor work or accuracy conclusion is formed
by this CPU stage.
