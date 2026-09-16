# Archived: global GFN2 feasibility pilot

Jacob stopped this route after reviewing its progress: “Oof. OK, then this is really not going to be feasible. File all this away; we can still use the baseline to start doing our benchmarking, right?” This supersedes the queued/running status in historical REPORT.md and RESULT.json. No retry remains authorized or scheduled by this workflow.

The 9,141-atom native ORCA GFN2/ALPB system exported metal parameters but produced **zero converged endpoint energies**. On the first attempt, 172 MPI ranks per endpoint replicated an excessive memory footprint. The 8-rank retry still reached 1,540,838,408 KiB batch peak RSS; La was killed with signal 9 during startup, and Ca was cancelled at Jacob's direction. No SCF throughput or viable production rate was established. This backend/configuration is unsuitable for affordable routine scoring. This does not rule out every global representation.

| Job | Outcome | Allocated CPUs | Top-level elapsed s | Allocated core-s |
|---|---|---:|---:|---:|
|1198999|cancelled during memory pressure; cleanup failed|344|866|297904|
|1199003|cancelled before execution|1|0|0|
|1199004|La startup failure; remaining work cancelled|112|841|94192|

Total recorded top-level allocation: **392,096 core-s**. Batch cleanup overlaps these records and must not be added as independent compute. Residual processes after Slurm lost containment on node-344-8t-1 are not fully accounted for; this is not an exact all-inclusive cost. The node remained drained for administrative cleanup; no node restart, permission change or undrain was attempted.

Inputs, retries, execution receipts, collection_1199004.json, terminal records and the saved boundary component audit remain under `workspaces/global_representation_20260915/`. Scratch was retained because failed cleanup could leave processes referencing it. No scientific output has been rewritten. The APBS challenger remains excluded following its approximately 10 kcal/mol partition discontinuity and other failed numerical checks. Neither challenger contributes to baseline benchmarking.

Recommendation: retain the baseline and benchmark its separately versioned peptide-amide repair. No global numerical score, predictive validation, or successful physical check is claimed.
