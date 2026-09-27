# Strict SCF test complete: attempted remedy rejected

Job1219450 completed four outputs in746s on344CPUs:256624allocatedCPU-s,
zeroGPU. Collector1219451 completed2s on1CPU; actual completion wake acknowledged.
Frozen parser rejects all4 endpoints' printed density/rotation summaries. All
runs switched to TRAH, so those summaries cannot establish final TRAH density
convergence. Keep the failed collection and report this interpretation limit.

Direct diagnostic parsing (not admitted new scores) shows rigid shifts remain
0.011142489kcal/mol(Ca),0.009014492(La). The stricter stopping changes endpoint
energies by at most0.000055kcal/mol. It did not address the dominant discrepancy;
close this remedy without more SCF/iteration increases. Detailed actual
solver evidence is in scf_check/FAILURE_DIAGNOSIS.md.

Next independent work is the already-prepared expanded Thr159 partition test,
explicitly a coarse representation diagnostic under PARTITION_DIAGNOSTIC_PLAN.md,
not passed numerical qualification, full hybrid relaxation, or affinity scoring.
The exact ML checkpoint remains unavailable. Production untouched.
