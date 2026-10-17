# Coupled La/Dy mechanics: current checkpoint

17 October2026UTC. The activated PLAN remains active and unmet. Root owns all launches and integration. PQQ production stays unchanged.

## What we learned

1. **Accessible motion matters.** Six broad-coordinate searches pass physical admission but none converges. The conditional Mex/LBT source discrepancy falls24.120991→6.209455kcal/mol. This compares unchanged models with more motion; scaffold Hamiltonians still differ. See broad_response/REPORT.md and BROAD_RESPONSE_CHECKPOINT.json.
2. **Matching scaffold treatment alone makes a modest difference.** The complete104-cell hybrid matrix retains the same four selected geometries. Source discrepancies are5.034374kcal/mol for the predeclared baseline region and5.125352 for the expanded region. Larger component changes compensate; small net region dependence does not validate each term. See lbt_hybrid_pool_recovery/result/REPORT.md and POOL_RECOVERY_CHECKPOINT.json.
3. **The new model produces different restoring forces.** Its12-step7CCO/La scout passes physical admission and lowers energy15.734267kcal/mol, mostly MM. Maximum force falls51.063006→15.450115kcal/mol/A; it remains nonconverged. See lbt_hybrid_response/scout_result/REPORT.md.

LBT is a peptide control, not a second natural LanM.7CCN is Lu-conditioned. Mex has an EF3 target, three Nd spectators and167waters. No binding free energy, robust preference classifier or experimental transfer is established.

## Owned active work

| Jobs | Work | Resources / gate |
| --- | --- | --- |
|13017–13020|Baseline hybrid full-coordinate response:7CCO La/Dy,7CCN La/Dy|Each oneH200/12CPUs; normal one-GPU queue;80steps/200evaluations per state.13017 includes the12 successful scoutsteps.|
|12977/12978|Independent embedded PBE0 external-charge La± derivative scout|12MPI each onnode64-28. Electronic residuals decrease. Actual force/cost gate precedes six remaining directions.|

Every job has an actual completion watcher; DFT also has manifest-aware health watchers. Exact commands, configurations and watcher PIDs are in workspaceSUBMISSION.json. Read HYBRID_RESPONSE_CHECKPOINT.json and mex_dft/README.md. All current molecular computation uses Slurm.

**Next:** inspect all four response terminals. If all physical/backend admissions pass, bind their exact results and cross-score the shared pool with at most four opposite-metal calls. The reviewed hybrid_response_comparison/COMPARE.py then separates fixed-coordinate model change, reselection and new response while keeping Mex fixed. FIGURE.py renders only a complete result. Saved-only force audit13031 failed its upstream exact raw derivative replay on node64/1CPU before producing rows. Source/config/native receipts are unchanged. Versioned13032 repeats only saved arithmetic with the original qualified node/12CPU environment; exact gates remain intact. It is queued normally with a completion watcher. Host/thread arithmetic is a hypothesis until measured, not a proven roundoff diagnosis. See hybrid_response_force_audit/failure_v1/.

First terminal13017 is admitted after80steps, nonconverged (max3.317400/RMS0.901487kcal/mol/A). Work from its hybrid start is−28.244859kcal/mol. It includes12scoutsteps and68new continuation evaluations; continuation allocation315s/12CPU/oneH200. 13018 also completes80steps/admitted/nonconverged (max1.788255/RMS0.536287), work−28.397116.13019 is unavailable after2acceptedsteps: its third trial embeddedGFN oscillates to500iterations; isolatedGFN passes. Do not pool partial survivors.13020 continues. A finite four-call same-Hamiltonian embedded solver diagnostic is being prepared; no restart or global fallback is launched.

## Numerical qualification and failures

Hybrid total/component derivatives pass the declared boundary and field directions at both steps and both regions (12986/12989). Numerical force consistency is distinct from electronic accuracy. The original12975 exact-map failure is preserved;12984 demonstrated only derived-Jacobian roundoff, with all coordinates/state exact. Versioned686eba84 permits bounded derived-coefficient roundoff only.

Original pool12992 preserves66 successful cells, one genuine isolated-GFN SCC failure and37 untouched cells. Four-call13009 qualifies one isolated-component numerical recovery without changing geometry/state/temperature/parameters/iterations/gates. Saved13012 validates its explicit certificate;13013 completes only the37 untouched cells with the original solver. All104 cells now pass. The original failure is immutable; this is not a global solver fallback or a unique-ground-state claim.

The primary baseline electronic region (complete residues3–14) was frozen before full-pool outcomes. Expanded2–15 remains sensitivity only. Existing actual raw data, states, caps, MM ledger, native parameter lineage and source-local geometry pools are retained.

## Limits and restart

The prior broad searches and current scout are not minima. No automatic optimizer-budget extension, new source selection, clamped-core restart, parameter fitting or production promotion follows from success. A new failed required cell leaves its dependent result unavailable.

Mex outside-QM motion is distributed: the largest3-residue window contributes only7.42% of peptide squared displacement. One local expanded-DFT region cannot validate global MM lowering; that branch is not launched. See mex_motion/FEASIBILITY.md.

Standard scheduling currently permits12CPU/job and120CPU/user; GPU partitions have separate limits. Use actual allocations and current policy, no memory/time/priority overrides. Existing unrelated jobs remain untouched. Read actual receipts before acting; delayed events never request repeated chemistry. Root owns status/submission. Agents own only disjoint source/review tasks.
