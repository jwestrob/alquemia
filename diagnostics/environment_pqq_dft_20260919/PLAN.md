# Native DFT check of charge-associated PQQ context shifts

Approved by root under Jacob's explicitly delegated next-round accuracy work:
“Proceed with frozen 8 DFT endpoints on all four charge−1 PQQ contexts and archived
core/neutral controls.” This is consumed development evidence, not blind testing.

## Question and frozen scope

The full compact native OMOL panel retains 25/25 PQQ calls, but the gap contracts
79.030608→2.315459 model-kcal. Every neutral expansion shifts its contrast downward;
all four contexts adding formal charge−1 shift upward. Real added anionic groups
can physically favor La; OMOL's global categorical charge conditioning can also
alter the descriptor. Existing outputs establish association, not cause.

Use all four charge−1 calibration contexts, with no favorable exclusions:
A0ACD6B9F2, A8R3S4, Q88JH5, Q9Z4J7. Run both Ca/La endpoints for each: eight new
native r2SCAN-3c/CPCM(Water)/DefGrid3 single points. Preserve each archived native
SP recipe and the exact complete-context XYZ previously used by OMOL; no nuclear
motion, new water, protonation, composition, cutoff, cap or reference change.
Reuse eight actual archived original-core DFT endpoints, and both original-core
and context outputs for neutral-expanded 1H4I/4MAE. Reuse corresponding native
OMOL core/context outputs. No new MACE calls, gradients, optimizations or training.

Record each endpoint's original/context energy and their difference, each Ca−La
contrast and context-minus-core contrast. Parse archived CPCM dielectric terms
and report their differences as diagnostics already included in the quantum
energy, never as extra terms or a unique decomposition of the environmental
change. Charge, composition, cavity and interactions all change together.
The endpoint energy change compares different compositions and is not itself
an interaction energy or an affinity. Raw contrasts and paired class gaps are
reportable; no old bands or universal zero are inherited. No threshold is fit.

## Execution

One ordinary64CPU/256GiB allocation, four16-rank tasks concurrently, using the
existing pinned native-ORCA runner and MPI/thread policy. Eight fixed tasks;
no arbitrary elapsed/CPU cap. Previous smaller context jobs took minutes and
one larger202-atom PQQ La endpoint took37minutes; this pilot has139–161atom inputs.
Count actual allocated resources, failed attempts and reuse. No production change.

## Interpretation and next decision

If native DFT reproduces the large relative charge-associated shifts, diagnose
physical context/representation dependence rather than automatically deleting
charge features. If it does not, the discrepancy motivates a carefully separated
model question; it still does not uniquely prove global conditioning is causal.
PQQ accuracy and structural robustness matter, not lower endpoint energy alone.
Preserve every output. Do not select an altered charge, context or threshold to
recover a favorable gap. Follow-on hypotheses require a new recorded experiment.
