# Q9 two-Ca preparation passed; preference remains unscored

All five real Q9 two-Ca sources are admissible under the unchanged local
PQQ chemistry and canonical context policy. The full source retains both native
Ca and all 4,849 heavy atoms with **zero coordinate change**, including after
protonation. Q9 sample3 correctly uses chain D's PQQ Ca and retains chain B's
structural Ca. The nearest retained context heavy atom is 23.79–23.99 Å from
the spectator, comfortably outside the frozen 4.3 Å exclusion gate.

Each resulting paired context has 192 atoms and the exact original canonical
state/protein signature (Ca charge −3, La charge −2). The actual one-Ca P12293
sample0 regression reproduces both core/context XYZ and expanded mapping/state
with zero displacement and identical XYZ hashes. This verifies that the narrow
source-inventory extension preserves the previous single-Ca behavior.

These are source preparation results. They do not yet establish better metal
classification, nor do they include the distant Ca in the local Hamiltonian.
The original unavailable one-Ca results remain unchanged; the new sources are
separate two-Ca occupancy hypotheses ready for matched native/solvent scoring.

Job2969 completed normally in 28 scheduler seconds on six CPUs (168 allocated
CPU-seconds); measured preparation wall time was 26.54 s and process CPU time
81.25 s. Batch MaxRSS was 738,908 KiB. No model forwards or electronic endpoint
calculations ran. Exact inputs, states, rows and process/accounting receipts are
in Q9_RESULT.json and spectator_preparation_q9_v2/execution/FINAL_PREPARATION.json.

The remaining ten C5/Q88 sources are now bound to the **same code/settings**,
reusing this passed regression receipt. Root submitted preparation2985. Their
actual results are pending; no additional preparation variant was introduced.
