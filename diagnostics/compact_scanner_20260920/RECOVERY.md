# MPI task-layout recovery — 2026-09-20

Job1203299 completed its wrapper but returned 0/4 complete scores. All8 MACE
scalar endpoints succeeded and matched archived energies exactly. All16 ORCA
attempts stopped in startup before SCF: OpenMPI saw one Slurm task slot under
`--ntasks=1 --cpus-per-task=32` and refused eight ranks. This is an allocation
layout failure, not electronic nonconvergence. Every failure remains preserved.

Recovery uses the same H200 node,32CPUs,1GPU,200000MiB, unchanged24-call science.
Request32tasks×1CPU, matching the existing working CPU runner convention. Run
each MACE pair as one explicit `srun --exact --ntasks=1 --cpus-per-task=32
--gres=gpu:1` step; it then has the real32CPU task allocation required by its
thread setting. The CPU stage retains the32-task allocation for four8-rank
ORCA processes. No oversubscription flag, forced MPI neutralization or shared
permissions/environment changes.

A fresh integrated attempt is necessary to measure complete score-path timing.
The successful8 MACE calls are therefore repeated for this explicitly recorded
recovery; no cached timing claim. Total across attempts is16MACE calls and32GFN2
launches, of which the first16 GFN2 never reached SCF. Initial failed allocation
cost43wall/GPU-seconds and1376core-seconds is added to final development cost.
The scientific manifest, geometry, model, solvent settings, published bands and
repeat tolerances remain unchanged. This is technical recovery under the
approved pilot, not another model experiment.
