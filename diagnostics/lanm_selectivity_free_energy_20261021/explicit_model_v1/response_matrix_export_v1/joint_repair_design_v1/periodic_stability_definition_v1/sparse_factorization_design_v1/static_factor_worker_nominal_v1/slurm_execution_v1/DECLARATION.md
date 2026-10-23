# Scheduling-only static-factor variant

Derive only from writable_indices_recovery_v1 RUN5af8de7f. Preserve the exact
matrix, candidate, PME/nominal shift, backend, mode qualification, writable CSC
indices, one symbolic/numeric attempt and maximum ten solves. Preserve 32 GiB
RLIMIT_AS, 900 s supervisor wall limit, 10800 CPU-s limit, OpenBLAS12/OpenMP1.

Change execution from biotite login to one real Slurm allocation: banfield,
standard, standard_expand, one node, one task, twelve CPUs per task, default
memory/time, no requeue, no exclusive-node request. The wrapper requests no
memory or time override. Current default4000MiB/CPU gives48000MiB total.

Before array access, verify the running job directly with scontrol, joining
SLURM_JOB_ID, actual BatchHost, user ID, account, partition, QoS, task/CPU/node
counts, AllocTRES and no-requeue/shared-node policy. Require the exact48000MiB
default grant, also explicitly greater than the32GiB process cap. Environment
runtime pins/imports and actual BLAS/OpenMP thread checks remain inherited and
occur before loading full matrix arrays.

Replace the login-only MemAvailable>=64GiB heuristic with the actual scheduler
grant check. Host free memory is recorded descriptively and is not called an
allocation. Slurm cgroups and the unchanged worker address-space cap govern
this shared-node execution. Record real job ID, compute node, CPU/memory grant,
raw scontrol readback, process CPU and wall costs. Scheduler allocated costs
remain pending sacct rather than being inferred from parent/batch duplicates.

Source preparation and finite stdlib controls only: synthetic allocation records
check missing/incorrect job/host/user/resources, requeue and scheduler policy;
AST/source checks retain the complete numerical worker after removing only its
new allocation metadata field. No scheduler submission, full matrix read,
symbolic/numeric factorization, solve, native call or package installation.
Preserve the original failed attempt and writable-index recovery controls.
