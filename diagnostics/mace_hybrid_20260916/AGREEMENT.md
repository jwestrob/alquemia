# Approved MACE capability and partition pilot

Date: 2026-09-16. Status: **approved; implementation/preparation in progress**.

Jacob's authorization, following the resource estimate and proposed 12-call
pilot: “I approve the experiment; if it shows we can't fit it in, we can reserve
multiple GPUs (even if we don't use them) to effectively reserve more RAM on the
H200 node in order to accommodate.” He also specified exactly one eighth of
the node's host RAM for each GPU share, and did not restrict the work to eight
CPU cores. Earlier instructions removed project CPU/wall-time stopping budgets.

## Scope

Execute stage A of PLAN.md: one consumed 1H4I catalytic chain A, its archived
protonation/PQQ(3-)/dry state and identical source coordinates; both Ca and La;
archived qm33/qm36 cores (47/54 atoms). The common full physical system has
9,141 atoms. Its formal charge is protein -8 + PQQ -3 + metal (+2/+3):
Ca -9, La -8, both singlets. Verify these statements against pinned artifacts
before admission. Caps occur only in the archived cores, never in the full system.

Use pretrained MACE-POLAR-1 medium, float64, isolated nonperiodic realspace
electrostatics, no external field, no solvent, no geometry changes and no
training. Pin software, checkpoint and actual charge/spin conventions. Reuse
the four converged native vacuum r2SCAN-3c endpoints; **zero new DFT**.

Twelve energy/force evaluations: four core states; two full states; repeat,
translate and rotate the full pair (six more). Translation is [10,-7,3] A;
rotation is 37 degrees about normalized axis [1,2,3], about the original metal
position. Transform forces back before comparing. Repeats must execute rather
than hit the primary cache. Original numerical results are already consumed;
this is method development, not prospective validation.

Outputs: endpoint energies/forces, direct and subtractive hybrid Ca-minus-La
contrasts, partition shift, repeat/rigid-transform checks, charge/state accounting,
software/input pins and actual wall/CPU/GPU/memory receipts, including failures.
No aquo reference, calibrated classification or solution-phase prediction.

## Frozen interpretation checks

Before observing MACE outputs, use absolute endpoint and paired-contrast
repeat/rigid-transform tolerance 0.01 kcal/mol; maximum transformed Cartesian
force discrepancy 0.001 eV/A; learned total-charge tolerance 1e-5 e. These are
numerical diagnostic tolerances, substantially smaller than the existing
multi-kcal chemical signals. The **2 kcal/mol** partition diagnostic is inherited
from the prior consumed boundary experiment for comparison, not fitted to this
pilot. Report exact errors and failures; passing is not an accuracy claim.

## Execution and recovery

One H200 share starts with 28 CPUs (224/8) and 2063701/8 MiB of configured
Slurm host memory; retain the exact request and any mandatory scheduler
rounding in receipts. Execute endpoints sequentially on one GPU. A memory-only
recovery may reserve additional GPU/CPU/host-memory shares proportionally, as
explicitly authorized. GPU memory is separate and does not pool automatically.

Technical retries, checkpointed differentiation or host offload may preserve
the exact same model/system/electrostatic expression. Record the implementation
and all failed attempts; verify any rewritten numerical kernel against the
original on these real core fixtures. No changed cutoff, float precision,
periodicity, model weights, atom inventory or smaller protein may be substituted
to obtain a successful result. Surface such a scientific change for agreement.

Stage B, new solvent models, training, relaxation, additional biological cases,
new DFT and default promotion remain outside this approval. Baseline preserved.
