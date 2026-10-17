# Coupled La/Dy mechanics: current checkpoint

17 October2026UTC. The activated PLAN remains active and unmet. Root owns all launches and integration. PQQ production stays unchanged.

## Latest: whole-peptide DFT challenges the extreme hybrid terminal force

Whole256 before-La13051 and original saved collector13283 are complete/admitted.
Along the frozen initial hybrid Phe1CA force direction, PBE0+sameALPB gives
3.341027 versus hybrid51.063006 and whole-MACE0.798547kcal/mol/A. Global force
difference RMS is10.141916 for hybrid versus3.511987 for whole-MACE. The large
hybrid terminal response is poorly supported by this conditional reference;
MACE remains quantitatively different and is not validated by this comparison.
Read [the actual scout report](lbt_whole_reference/scout_result_v1/REPORT.md).

Native18-cycle termination/full256 gradient and exact collector replay pass.
Four printed residual flags and all electronic-state/force-accuracy limitations
remain explicit. Cost284892CPU-s plus3CPU-s for collection. Root approved the
original remaining three fixed cells after this measured-cost/force review,
with no repeated scout, retries or optimization. Exact decision and config are
in workspaces/lanm_coupled_mechanics_20261017/lbt_whole_reference_continuation_v1/.
Post-scout source/config review, actual Slurm preflight13288 and all seven
negative controls13289 pass (5CPU-s total). The original remaining cells now
run as13290afterLa/13291beforeDy/13292afterDy,12MPI each, with completion and
actual singleton-manifest health wakes. Native output creation is verified;
all three printed molecular state/ECP/MaxIter settings and actual12MPI/MaxCore3100
match the frozen inputs. Independent startup review passes; terminal admission
remains pending. No waiting collector is allocated.
Read [the execution checkpoint](lbt_whole_reference_continuation/activation_v1/REPORT.md).
Before-Dy13291's first-build silence alert was reviewed against advancing CPU
accounting and the completed scout's1949.7-second first iteration. Continue
unchanged; convergence remains unestablished. The checkpoint preserves the
actual alert, scheduler observations and independent review.
After-Dy13292's matching first-build alert was separately reviewed with its own
output and increasing CPU accounting; continue unchanged, with no recovery.
After-La13290 has the same reviewed first-build status; all three continuation
jobs remain active without a justified cancellation or retry.
All three initial silence alerts have now resolved: first SCF iterations are
printed (after-La1967.8s, before-Dy1948.6s, after-Dy1969.1s). Exact rows and live
accounting are in the execution checkpoint's `first_iteration_review_v2/`.
Mex central13281 has reached iteration6. All four remain active without terminal
admission; no new mechanics or preference result is established by this progress.

## Latest: Mex mismatch localized to final integration; original gate still fails

Saved diagnosis13277 completed with exact original replay and zero molecular
calls. Printed final COSX contributes -0.401606 of the -0.414201 kcal/mol/A
derivative mismatch. Native point-charge gradient mapping reproduces exactly.
Subtracting the printed correction leaves -0.012595; subtracting the full
post-SCF shift leaves -0.011124. Both exceed the old0.005409 tolerance.
The post-SCF shift includes an explicit remainder beyond printed rounding;
neither diagnostic energy definition is promoted. Read
[the saved diagnosis](mex_dft/gridx3_diagnosis_v1/REPORT.md).

Cost2CPU-s, no new molecular calls. The original remaining6 qualification and8
response cells remain stopped. A separate consistent-grid qualification is
specified under mex_dft/consistent_grid_v1/EXECUTION_SPEC_V1.md: fresh central
gradient first, then at most two displaced energies after actual admission/cost
review. Only SCF grid2 is brought to the original final grid3 construction;
original derivative tolerance and state/SCF settings remain. Root accepted this bounded design. Actual saved-only preflight13280 passes all
source, shared-component and three legacy-parser controls with zero molecular
calls (10CPU-s). Failed13279 remains preserved (10CPU-s); versioned adapter1600a0a5
corrects only a native text delimiter. Root and independent review admit the
prepared inputs. Fresh central13281 is submitted on12MPI, config28abed13,
completion watcher2759506 and alert-only health armer2759507. Actual native startup
passes full grid/state/runtime checks: grid2 equals original/final grid3; XC/grid1/
grid3 stay unchanged. Charge0/multiplicity1/NEL832, MaxIter125 and12MPI/MaxCore3100
match. Terminal gradient admission and measured cost still gate both displacements.
No plus/minus or original-six/response-eight cells are submitted. Exact receipts:
mex_consistent_grid_central_v1/SUBMISSION.json; see the consistent-grid report.
Whole13051 and its original collector13283 have completed as recorded above.

## Mex derivative scout failure

Both Mex endpoints12977/12978 terminate normally and pass native admission.
Original saved collector13271 fails PBE0/total directional consistency: embedded
PBE0 analytic -3.673754 versus finite difference -4.087955 kcal/mol/A;
residual -0.414201 exceeds 0.005409. MM and fixed-charge ALPB pass.
Read [the failed scout report](mex_dft/scout_pair_failure_v1/REPORT.md).
The original stop rule closes the remaining six qualification points and
blocks the dependent eight-cell response panel. No solver/step retry is launched.

Independent source review finds no mapping/sign/unit defect. The saved
COSX/gridx3 decomposition is complete under mex_dft/gridx3_diagnosis_v1/;
it makes no molecular calls and cannot override the original failed gate.
Whole-peptide reference13051 subsequently completed; its result and continuation
are recorded above. The goal remains unmet.

## What we learned

1. **Accessible motion matters.** Six broad-coordinate searches pass physical admission but none converges. The conditional Mex/LBT source discrepancy falls24.120991→6.209455kcal/mol. This compares unchanged models with more motion; scaffold Hamiltonians still differ. See broad_response/REPORT.md and BROAD_RESPONSE_CHECKPOINT.json.
2. **Matching scaffold treatment alone makes a modest difference.** The complete104-cell hybrid matrix retains the same four selected geometries. Source discrepancies are5.034374kcal/mol for the predeclared baseline region and5.125352 for the expanded region. Larger component changes compensate; small net region dependence does not validate each term. See lbt_hybrid_pool_recovery/result/REPORT.md and POOL_RECOVERY_CHECKPOINT.json.
3. **The new model changes terminal bonded mechanics.** Saved force audit13032 passes all four states. The largest force change is Phe1 CA, outside the electronic region, with hybrid maximum forces47–53 kcal/mol/Å at the old selected geometries. Saved MM groups make this a bonded-led change; no obvious duplicate interaction or charge-shift error was found. This does not establish which model is accurate. See hybrid_response_force_audit/result_v2/.
4. **Original response jobs are terminal.** Three finish80 accepted steps and remain nonconverged;7CCN/La fails embedded GFN after2 steps. The three admitted works are−28.245/−28.397/−31.301 kcal/mol. There is no complete paired response pool. See lbt_hybrid_response/original_four_state_result/.
5. **The finite numerical recovery failed at a later geometry.** Four-call13036 passes its three frozen agreement gates, and preflight13047 exactly reconstructs the original failed trial. Continuation13048 resumes correctly but fails embedded GFN2 again after10 accepted steps/11 evaluations. The new failed trial reaches500 SCC iterations with large residuals; isolated GFN2 converges. The frozen second-failure rule closes this branch. No four-state pool or partial-survivor substitute is available. See lbt_hybrid_response_recovery/failure_result/REPORT.md.

6. **The terminal force is specifically bond-and-angle dominated.** MM-only audit13122 passes all four original starts and all208 local terms, with exact native parameter/term checks and fullforce closure. The largest term in every case is N–CA stretching, followed by N–CA–C bending. Along the total MM force, bonds contribute24.96–30.62 and angles16.70–22.52kcal/mol/A; torsions/nonbonded are smaller. This identifies the modeled cost, not the accurate model. No copying/counting defect found. Audit+saved report13130 cost10allocatedCPU-s, zero other molecular calls. Read lbt_mm_term_audit/result/REPORT.md and ALL_TERMS.tsv.

7. **The two models disagree on the cost of exactly the same motion.** Whole-model completion 13140–13142 passes all native and saved replay gates. The 13017 La deformation lowers hybrid energy by 28.244859 kcal/mol but raises whole-MACE plus fixed-FF ALPB by 11.290076. Shared solvent cannot explain this disagreement; differing electronic regions and MM mechanics prevent unique attribution to stiffness. Whole-model Dy-minus-La deformation work is -0.207676 kcal/mol, below the existing 1 kcal reporting resolution. Read [the report](lbt_whole_model_completion/REPORT.md). Four energies and three force sets are available; hybrid after-Dy remains missing. No accuracy or preference claim follows.

LBT is a peptide control, not a second natural LanM.7CCN is Lu-conditioned. Mex has an EF3 target, three Nd spectators and167waters. No binding free energy, robust preference classifier or experimental transfer is established.

## Owned work

| Jobs / task | Work | Resources / gate |
| --- | --- | --- |
|12977/12978 complete;13271 failed derivative gate|Embedded PBE0 external-charge La± derivative scout|Both endpoints available; PBE0/total derivatives fail, MM/ALPB pass. Remaining six stopped.|
|13051/13283 complete|Whole256/q−2 PBE0 before-La scout and saved collector|Native admission and exact saved replay pass; one new reference cell.284895allocatedCPU-s total.|
|13281|Mex central La gradient with consistent SCF/final COSX grids|Actual zero-call preflight13280 and native startup checks pass.12MPI; terminal gradient admission and measured cost gate the two displaced calculations. Completion and alert-only health wakes armed.|
|13290/13291/13292 running|Three remaining whole-LBT reference cells|Original source/settings unchanged. Actual13051/13283 admitted; preflight13288 and seven controls13289 pass.12MPI per singleton, completion/health wakes armed.|
|Stopped at prerequisite|Mex eight-cell DFT response adapter|Prepared source/tests remain preserved. Original scout13271 fails derivative qualification; no response cells submitted.|
|13122/13130 complete|Four-start MM term attribution and saved report|4/4 states,208 terms, all closures pass. N–CA stretching and N–CA–C bending dominate.10allocatedCPU-s total. See lbt_mm_term_audit/result/REPORT.md.|
|13140/13141/13142 complete|Whole-LBT fixed-deformation model completion|Native admission and saved replay pass. Two MACE plus one ALPB calls; four energies/three force sets. 140 CPU-s/11 GPU-s. DFT gates unchanged.|

Read ACTIVE_CHECKPOINT.json for exact active submissions and watcher records. All four original response jobs13017–13020 and recovery13048 are terminal. Recovery costs151s on12CPUs/oneH200; partial values remain diagnostic only. Solver qualification13033/13036 and saved force audit13032 are terminal. No paired pool or figure is submitted. Root owns submissions; agents own adapters/reviews. All numerical science uses Slurm.

Historical single-endpoint checkpoint, before plus completion: Mex minus 12978 completed on 17 October at 09:38 UTC: 32 SCF cycles, normal native termination, adapter status complete/DFT available. Its final gridx3 exchange step took 1,920.673 s; the long silence resolved without intervention. Scheduler cost: 22,132 s × 12 CPUs = 265,584 allocated CPU-s. Native energy-based convergence does not satisfy every printed density/gradient criterion; all flags and physical-state limitations remain. Read [the single-endpoint checkpoint](mex_dft/scout_minus_result_v1/REPORT.md). No paired derivative or mechanics accuracy result is established. Plus 12977 and whole-peptide 13051 remain active; no dependent reference is launched.

**Next:** collect the whole-LBT four-cell comparison after actual13290–13292 terminal results. Preserve all missing/failed slots; no repeated scout or retries. For Mex13281, review the actual terminal central gradient and measured cost before either displaced calculation. Preflight13280 and central submission are already complete; do not repeat them. No displaced Mex or original response cells are submitted.

The saved [startup timing review](mex_dft/consistent_grid_v1/startup_review/PHASE_TIMING_REVIEW.json)
records that the prior Mex first SCF iterations took about20minutes after startup
and guess. Output silence during the new first build alone does not establish a
stall. Both jobs were RUNNING with advancing CPU accounting at12:00UTC; whole13051
was calculating the COSX gradient. These are activity observations, not convergence
or force-accuracy evidence. Existing health/completion wakes remain active.

The saved audit's first CPU-host attempt13031 failed an upstream exact derivative replay before producing rows. Unchanged source/config/gates pass on the original qualified host/thread setup as13032. Cross-host exact-arithmetic portability remains unresolved; use the qualified replay environment, and preserve the first failure. See hybrid_response_force_audit/failure_v1/ and result_v2/.

## Numerical qualification and failures

Hybrid total/component derivatives pass the declared boundary and field directions at both steps and both regions (12986/12989). Numerical force consistency is distinct from electronic accuracy. The original12975 exact-map failure is preserved;12984 demonstrated only derived-Jacobian roundoff, with all coordinates/state exact. Versioned686eba84 permits bounded derived-coefficient roundoff only.

Original pool12992 preserves66 successful cells, one genuine isolated-GFN SCC failure and37 untouched cells. Four-call13009 qualifies one isolated-component numerical recovery without changing geometry/state/temperature/parameters/iterations/gates. Saved13012 validates its explicit certificate;13013 completes only the37 untouched cells with the original solver. All104 cells now pass. The original failure is immutable; this is not a global solver fallback or a unique-ground-state claim.

The primary baseline electronic region (complete residues3–14) was frozen before full-pool outcomes. Expanded2–15 remains sensitivity only. Existing actual raw data, states, caps, MM ledger, native parameter lineage and source-local geometry pools are retained.

## Limits and restart

Experimental LBT terminal chemistry has already had a complete bounded primary-source review: [paper, full supplement and deposition audit](../lanm_coupled_response_20261016/LBT_TERMINAL_CHEMISTRY.md). Those records do not establish the synthesis terminal groups and retain a charge-declaration discrepancy. The free NH3+/COO− model remains conditional. Repeating those same sources cannot resolve the missing experimental specification; this uncertainty does not change the ongoing fixed-state force comparison or authorize a new terminal-state campaign.

The prior broad searches and current scout are not minima. No automatic optimizer-budget extension, new source selection, clamped-core restart, parameter fitting or production promotion follows from success. A new failed required cell leaves its dependent result unavailable.

Mex outside-QM motion is distributed: the largest3-residue window contributes only7.42% of peptide squared displacement. One local expanded-DFT region cannot validate global MM lowering; that branch is not launched. See mex_motion/FEASIBILITY.md.

Standard scheduling currently permits12CPU/job and120CPU/user; GPU partitions have separate limits. Use actual allocations and current policy, no memory/time/priority overrides. Existing unrelated jobs remain untouched. Read actual receipts before acting; delayed events never request repeated chemistry. Root owns status/submission. Agents own only disjoint source/review tasks.

Saved preparation13059 failed an upstream exact terminal physical/force replay before molecular calls. Its Python binary matched, but NumPy dependency metadata differed (1.26.4 versus original2.4.6). Recovery13060 exposed a missing Gemmi dependency before replay.13063 keeps the original environment and supplies only the exact existing Gemmi package/native libraries in an owned overlay. Neither original environment is modified. The specific terminal numerical field remains unresolved until actual evidence. See mex_dft_response/preparation_failure_v1/.

Mex response adapter delivery: mex_dft_response/READINESS_REPORT.md. All six preparation/diagnosis/preflight allocations are terminal and preserved; no molecular calls were made. The original full physical/source/force checks and native acceptance rules remain unchanged. Mex12977/12978 are complete and scout collection13271 failed its derivative gate. The separate consistent-grid central scout13281 and whole-LBT remaining13290–13292 are active; the original response panel remains stopped.

## Fixed-deformation model completion — 17 October2026UTC

A bounded cache audit found missing whole-MACE energies at the exact 13017 La terminal for both metals, plus its Dy fixed-FF ALPB. Root froze lbt_whole_model_completion/SPEC.md before execution; independent source review passed. Preflight 13140 passed with zero molecular calls. Run 13141 completed exactly those three evaluations on 12 CPUs/one H200. Saved collector 13142 replayed actual native outputs, input states, full force arrays and the comparison. Both before points and after-La solvent were reused. Energy coverage is 4/4; whole-MACE force coverage is 3/4 because archived before-Dy was energy-only. The missing hybrid after-Dy point stays unavailable.

The same La displacement has hybrid work -28.244859 versus whole-model work +11.290076 kcal/mol. This establishes Hamiltonian dependence, not which model is accurate or a unique MM defect. Total cost is 140 allocated CPU-seconds and 11 GPU-seconds. Full results, source and actual commands are preserved under lbt_whole_model_completion/; exact native pins are in ACTIVE_CHECKPOINT.json. Independent integration review agrees.

The ongoing before-La DFT scout tests the initial forces and cost. It alone cannot validate the finite deformation work; the existing admission/cost gates still control further reference cells. No closed-response recovery, new geometry or DFT continuation was launched by this adjunct. All three completion watchers are recorded; delayed callbacks require acknowledgement, never repeated chemistry.

## Historical source-preparation checkpoint — superseded by activation above

The original runner permits only before-La. A separate [continuation draft](lbt_whole_reference_continuation/CONTINUATION_DRAFT.py) now implements the three already-declared remaining cells as singleton runs and a saved four-cell collector. Read [SPEC_DRAFT.md](lbt_whole_reference_continuation/SPEC_DRAFT.md) and STATIC_DRAFT_REVIEW.json. Author AST inspection and independent static source/provenance review pass; no scientific imports, numerical tests, runnable configs, continuation decision or submissions were made.

The draft replays original13049 unchanged, then separately admits the completed13142 after-Dy solvent. It forbids repeating before-La, retains unavailable comparator slots and uses the original method/state/admission rules. Actual13051 completion and its original saved collection, root force/cost review, a bound continuation decision and real-artifact Slurm checks are still required before execution. No physics result or execution readiness is implied by source review. Root owns those gates and all future launches; source preparation is complete.

The separate saved-artifact validation harness is also prepared and independently
reviewed: `VALIDATE_DRAFT.py`, interface note and `STATIC_VALIDATION_REVIEW.json`
in the same continuation directory. It requires a genuine positive replay before
seven rejection checks; supplemental geometry/metal checks explicitly exercise
the unchanged native admission helper. AST/source review only: no scientific
imports, executed tests, bound configuration, decision or submission. Run the
actual checks through Slurm only after the original scout/collector and cost
review. Native continuation execution and the later collector need their own
actual validation; this source review does not establish either.
