# Saved LBT partition forces against the whole reference

## Question and fixed denominator

At the same 7CCO/cartesian_01 coordinates, does changing the existing baseline
partition to the existing expanded partition change the terminal force mismatch
against whole-peptide PBE0? Inspect exactly four archived hybrid cells:
baseline/expanded × La/Dy. The reference is the admitted before-La/before-Dy pair
from full collector 13424. Both metals use the same physical coordinates.

The baseline includes full residues3–14; expanded includes2–15. The left cap
therefore changes from C2–CA2 to C1–CA1. Phe1 CA is an ordinary MM atom in baseline
and the omitted cap anchor in expanded. The N1–CA1 bond stays wholly MM; the
N1–CA1–C1 angle stays retained, classified cross-unrepresented in expanded.
This check reports existing forces, not a new partition or revised ledger.

Inputs are fixed in READY.json before execution: collector
`9229f7bebb533dd44a024bb80eca10eded2b3a64cc3972b245d65d82b5ecda97`,
recovered104-cell pool `d98aafb629185e785a3bd291802e46c374cfeb35b46b52486147c51303f1d592`,
the four exact POINTs, immutable mapping helper, original frozen directions,
physical maps, state and component lineage. No favorable subset or survivor
comparison replaces the four-cell denominator.

## Execution and admission

One Slurm CPU; no GPU, executable subprocess, model load, molecular evaluation,
new geometry, finite difference, optimizer, fit or repeated52-term MM audit.
Use the original MACE environment with NumPy2.4.6 and single-thread BLAS/OpenMP.
Only NumPy and the existing partition helper's pure saved-map `scatter` are
used; its `initialize`, geometry construction and molecular routines are not
called. A process audit tripwire prevents subprocesses. Fresh output directory;
every failed slot and traceback is retained, with overall incomplete status.

Require full collector admission, original four-cell recovered pool admission,
exact POINT/state/source/coordinate/physical-ID linkage, full charge−2,
physical MACE multiplicities La1/Dy6, GFN effective singlets, and exact common
fixed-FF charge-vector/solvent receipts. Reference frozen-f state distinctions
and native SCF caveats remain unchanged. The archived reference collection is
the admission authority; this is not a new DFT parser or native qualification.

Replay each hybrid full256 derivative from pinned MACE, embedded/isolated GFN,
MM and ALPB component receipts through the original archived cap Jacobians.
Reuse PREPARE_V2.scatter, including both cap anchors. Preserve
HybridEvaluatorV2.restore's **exact** component order and array-equality closure;
there is no newly relaxed tolerance. Native energies are checked only for saved
lineage/closure, not compared between partitions or methods.

Preserve each producer's floating-point order: hybrid uses `raw*HA/BOHR`, while
REFERENCE.admit_alpb uses `raw*(HA/BOHR)`. Replay each exactly against its archive;
require identical raw ALPB gradient and receipt and report the conversion-order
difference separately. Reproduce both original baseline force tables and their
projections exactly. The before-Dy baseline is at cartesian_01, not the different
own-Dy start used by the earlier four-start MM audit.

## Outputs and interpretation

Retain all256 force/error rows per cell, reference and hybrid vector RMS,
maximum atom norm/index, mean atom norm and vector sum. Errors are
`hybrid force − reference force`, in kcal/mol/Å. Report both already-frozen
full physical directions: full_deformation_chord and initial_hybrid_Phe1_CA_force.
Phe1 CA component projections use each original direction's atom3 vector,
without local renormalization; they are contributions to those projections.
Retain cap/region classification, all component vectors and native state.

Expanded-after remains explicitly unavailable. The result can distinguish a
persistent MM force imbalance from partition-dependent compensation. It does
not identify a unique bad parameter, establish gradient/physical-state accuracy,
rank metal preferences, authorize optimization, or reopen the closed hybrid
response branch. Original reference vacuum-dianion, terminal, frozen-f and SCF
limitations remain attached. No automatic numerical follow-up follows a failure.

Source/AST/metadata review alone occurs on login. Root owns independent review,
actual Slurm execution, cost receipts and interpretation.
