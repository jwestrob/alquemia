# PModel scout — out of memory after successful initialization

21 October UTC: **16611 OUT_OF_MEMORY**,962s on24CPUs/96000MiB (23,088CPU-s).
PModel and its native La atomic density completed; requested beta514/515 rotation
printed. Main LEANSCF warned3100MB/rank available vs4650.8MB estimated, then
Slurm killed a rank at first Fock formation. No molecular energy or forces.
Read ROOT_TERMINAL_REVIEW.json. Preserve original outputs and failed admission.

A separate resource-only recovery is now running as16680:12MPI ranks in a
24CPU allocation, same96000MiB, one numerical thread/rank, MaxCore6200 and an
unchanged scientific template. Preparation16667 passed; actual runtime passed.
Read ../lbt_omol_memory_recovery_v1/CURRENT.md. This failed branch is closed;
no unchanged retry or duplicate submission.

## Historical startup checkpoint

# PModel recovery — central scout running

21 October UTC: **16611 RUNNING**,24MPI,96,000MiB,MaxCore3100MB/rank,
standard_expand/default memory and time,noGPU. Completion66459 and health
watcher66460 are active and recorded in the workspace/run directory. One requested central
EnGrad; native atomic-density initialization is allowed within this endpoint.

Run submission archiveefe8b8ee preserves source, release, runtime and watcher
evidence; preparation archive1e51bb64 preserves the actual zero-call result.

Root startup7c423afe and independent startupb597fd8e confirm exact runtime
input delta and native SCNL/LaECP46; actual guess/rotation/SCF still pending.

Preparation16609 PASSED on1CPU,7CPU-s,zero engine calls. RESULT95e1852b,
manifestd0ed827a,input evidence48ee5022,archive1e51bb64. Independent static
review76932bda and actual813440eb; root actual7eceb21a. Source27ac4e80,
SPEC183a7f10,run configb654d846 and root release6df38e02 are fixed.

The complete template changes only PAtom→PModel; XYZ and external basis are
byte-identical to16602. All other electronic settings, rotation and physical
state stay fixed. Actual SCNL is expected from the unchanged EnGrad request;
review native method/guess/rotation/basis/ECP/spin/residuals at completion.
Changing the initial guess can change the converged electronic solution.
This is a derived-recipe variant, not exact historical training-label replay.

Failed16602 remains preserved (initial-guess La orbital-data failure; no final
energy or forces;12,936CPU-s; archive4fee749e). No new comparison exists yet.
No duplicate endpoint, optimizer, displaced pair or automatic retry. Root owns
all launches; terminal success still needs actual method/state review.
