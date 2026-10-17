# Coupled La/Dy mechanics: current checkpoint

17 October2026UTC. The activated PLAN remains active and unmet. Root owns all launches and integration. PQQ production stays unchanged.

## What we learned

1. **Accessible motion matters.** Six broad-coordinate searches pass physical admission but none converges. The conditional Mex/LBT source discrepancy falls24.120991→6.209455kcal/mol. This compares unchanged models with more motion; scaffold Hamiltonians still differ. See broad_response/REPORT.md and BROAD_RESPONSE_CHECKPOINT.json.
2. **Matching scaffold treatment alone makes a modest difference.** The complete104-cell hybrid matrix retains the same four selected geometries. Source discrepancies are5.034374kcal/mol for the predeclared baseline region and5.125352 for the expanded region. Larger component changes compensate; small net region dependence does not validate each term. See lbt_hybrid_pool_recovery/result/REPORT.md and POOL_RECOVERY_CHECKPOINT.json.
3. **The new model changes terminal bonded mechanics.** Saved force audit13032 passes all four states. The largest force change is Phe1 CA, outside the electronic region, with hybrid maximum forces47–53 kcal/mol/Å at the old selected geometries. Saved MM groups make this a bonded-led change; no obvious duplicate interaction or charge-shift error was found. This does not establish which model is accurate. See hybrid_response_force_audit/result_v2/.
4. **Original response jobs are terminal.** Three finish80 accepted steps and remain nonconverged;7CCN/La fails embedded GFN after2 steps. The three admitted works are−28.245/−28.397/−31.301 kcal/mol. There is no complete paired response pool. See lbt_hybrid_response/original_four_state_result/.
5. **The finite numerical recovery failed at a later geometry.** Four-call13036 passes its three frozen agreement gates, and preflight13047 exactly reconstructs the original failed trial. Continuation13048 resumes correctly but fails embedded GFN2 again after10 accepted steps/11 evaluations. The new failed trial reaches500 SCC iterations with large residuals; isolated GFN2 converges. The frozen second-failure rule closes this branch. No four-state pool or partial-survivor substitute is available. See lbt_hybrid_response_recovery/failure_result/REPORT.md.

LBT is a peptide control, not a second natural LanM.7CCN is Lu-conditioned. Mex has an EF3 target, three Nd spectators and167waters. No binding free energy, robust preference classifier or experimental transfer is established.

## Owned active work

| Jobs / task | Work | Resources / gate |
| --- | --- | --- |
|12977/12978|Embedded PBE0 external-charge La± derivative scout|12 MPI each on node64-28; completion and health wakes armed. Actual derivative/cost result gates six remaining directions.|
|13051|Whole256/q−2 PBE0 before-La force/cost scout|Actual preflight13049 passes all four states,1030 explicit electrons. Only one new reference cell,12MPI. Completion and alert-only health wakes armed.|

Read ACTIVE_CHECKPOINT.json for exact active submissions and watcher records. All four original response jobs13017–13020 and recovery13048 are terminal. Recovery costs151s on12CPUs/oneH200; partial values remain diagnostic only. Solver qualification13033/13036 and saved force audit13032 are terminal. No paired pool or figure is submitted. Root owns submissions; agents own adapters/reviews. All numerical science uses Slurm.

**Next:** collect the one whole-peptide DFT force scout to test the large terminal response without inheriting MM stiffness. This is an explicit conditional vacuum/frozen-f test, not a new source or state search. Collect the Mex derivative scout before its remaining checks or response references. The Mex eight-cell response inventory is reviewed; a source-only adapter is being prepared with execution blocked until actual derivative/admission/cost gates pass. No further DFT cells or SCC repairs are launched.

The saved audit's first CPU-host attempt13031 failed an upstream exact derivative replay before producing rows. Unchanged source/config/gates pass on the original qualified host/thread setup as13032. Cross-host exact-arithmetic portability remains unresolved; use the qualified replay environment, and preserve the first failure. See hybrid_response_force_audit/failure_v1/ and result_v2/.

## Numerical qualification and failures

Hybrid total/component derivatives pass the declared boundary and field directions at both steps and both regions (12986/12989). Numerical force consistency is distinct from electronic accuracy. The original12975 exact-map failure is preserved;12984 demonstrated only derived-Jacobian roundoff, with all coordinates/state exact. Versioned686eba84 permits bounded derived-coefficient roundoff only.

Original pool12992 preserves66 successful cells, one genuine isolated-GFN SCC failure and37 untouched cells. Four-call13009 qualifies one isolated-component numerical recovery without changing geometry/state/temperature/parameters/iterations/gates. Saved13012 validates its explicit certificate;13013 completes only the37 untouched cells with the original solver. All104 cells now pass. The original failure is immutable; this is not a global solver fallback or a unique-ground-state claim.

The primary baseline electronic region (complete residues3–14) was frozen before full-pool outcomes. Expanded2–15 remains sensitivity only. Existing actual raw data, states, caps, MM ledger, native parameter lineage and source-local geometry pools are retained.

## Limits and restart

The prior broad searches and current scout are not minima. No automatic optimizer-budget extension, new source selection, clamped-core restart, parameter fitting or production promotion follows from success. A new failed required cell leaves its dependent result unavailable.

Mex outside-QM motion is distributed: the largest3-residue window contributes only7.42% of peptide squared displacement. One local expanded-DFT region cannot validate global MM lowering; that branch is not launched. See mex_motion/FEASIBILITY.md.

Standard scheduling currently permits12CPU/job and120CPU/user; GPU partitions have separate limits. Use actual allocations and current policy, no memory/time/priority overrides. Existing unrelated jobs remain untouched. Read actual receipts before acting; delayed events never request repeated chemistry. Root owns status/submission. Agents own only disjoint source/review tasks.
