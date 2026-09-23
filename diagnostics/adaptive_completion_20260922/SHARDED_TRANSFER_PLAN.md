# Primary225 execution: four fixed independent task shards

Root approved full transfer after the actual original30 common pool completed
and its canonical-only reference became available (25/25, gap2.54669018555
model kcal/mol). The frozen reference is
`workspaces/adaptive_completion_20260922/original30_pool_v1/REFERENCE.json`.
This is an execution-layout change only; no scientific policy or candidate changes.

Preserve the unlaunched `primary225_v1` manifest. Prepare
`primary225_v2_sharded` with the updated immutable implementation and identical
410 source tasks,225 declared cases and450 endpoint denominator. Assign task
index modulo4 to four disjoint searches:103,103,102,102 endpoints. Each task
keeps exactly one original-q0 start, identical selected four angular coordinates,
optimizer/unit scale/bounds, native MACE energy and analytic forces. There are
no additional candidates, folds, chemical states or DFT/GFN2 calls.

Each shard has a separate lock and exact selected-index/task receipt. Its warm
calculator uses job-specific logs; task-specific proposal/evaluation directories
remain disjoint. Four separate Slurm jobs request one H200,32 CPUs and200000 MiB
each. Existing scheduler policy controls concurrency; no job priority or another
allocation is changed. The aggregate active maximum is4 GPUs,128 CPUs and800000
MiB. This trades parallel resources for elapsed completion time, not fewer
scientific evaluations. Actual total cost is the sum of the four allocation
CPU/GPU seconds; concurrent elapsed time is not substituted for that total.

Collection remains225/450 including all20 inherited unavailable cases. No
calibration refit uses transfer outcomes. Preserve all failed trials/searches
as explicit results. Root owns subsequent five-geometry common-pool scoring
and comparisons. Existing baseline/production and old experiments are unchanged.
