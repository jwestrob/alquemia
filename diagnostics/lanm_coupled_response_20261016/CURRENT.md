# Latest: range-only first scout submitted

12544 runs7CCO/La continuation with oneGPU/12CPUs and completion watcher1944984.
Frozen adapterc61a12b1/config7ac515db under lbt3_range_scout_v1. Root reviewed
warm-start/guards/own-energy loop; independent review passed. Internal real-artifact
preflight checks all four old starts before model loading. Reuse exact saved
energy/forces and original-source references; .40Å source envelope is the only
search change. No continuation result yet. Remainingthree follow on execution
and physical admission, not energy sign. Root owns submissions; Harvey owns
extended commonpool, Leibniz has separate read-only Mex transfer design.

# Latest: range-only continuation declared; implementation underway

Read LBT3_RANGE_EXTENSION_PLAN.md. All four prior searches and complete paired
pool are terminal/reused; saved diagnostic12542 identifies the source .20Å
Asp7 oxygen displacement boundary in every state. The new version changes only
that original-source envelope to .40Å, retains target .20Å/torsion .20rad and
all physical/state guards, and starts from each exact saved terminal point.
Leibniz owns new continuation source; Harvey owns review/union-pool source.
Root owns submission. No new job or extended-range result yet. First7CCO/La
scout, then otherthree on execution/physical admission irrespective of sign.
PQQ and original failures remain intact; goal remains active and unmet.

# Latest: shared-pool response pilot complete; scientific extension being selected

12541 completes4new cross-metal cells/12totalmatrixcells. Source-local paired
spread72.966739→57.895436kcal/mol; globalmenu La7CCO-derived/Dy7CCN-derived persists.
Read lbt3_accommodation/REPORT.md and POOL_RESULT.json. No robust preference or
binding-free-energy claim. 12542 saved diagnostic complete: every search blocked
atAsp7O sourceextent.20A, no damagedbonds/newseverecontacts. Significant omitted
donor/backbone loads and distinct internal carboxylate forces remain.
All owned molecular jobs terminal. Root selects next finite physical-response
extension; no newbound/search/reference authorized by delayedcompletionalone.
PQQ untouched. Preserve allcompleted states; no repeat of12533/12538/12541.

# Latest: all four accommodation searches complete and physically admitted

12538 completed three unchanged searches in225s/12CPUs/oneH200; all final
physical gates pass. Together with12533, all four own-metal candidates exist.
Every final trial set is blocked by the original0.20A source envelope; forces
remain. Read lbt3_accommodation/SCOUT_REPORT.md and REMAINING_REPORT.md.
No repeated searches. All scientific jobs terminal while the narrow four-cell
opposite-metal common-pool adapter finishes; exact cell/geometry reuse required.
Root owns execution. Separate saved-force/limiting-atom diagnostics will guide
any extension; no expanded search or changed bound has been launched.

# Latest: remaining three accommodation searches submitted

12538 runs7CCO/Dy,7CCN/La,7CCN/Dy sequentially with one warm GPU/12 CPUs,
unchanged search physics and fresh per-state history/geometry/cache. Completion
watcher1922524; implementationba41b0d/configbf84b806 in
workspaces/lanm_coupled_response_20261016/lbt3_accommodation_remaining_v1/.
Root and independent reviews passed actual state/counter scope before submission.
Scout12533 is complete and must be reused. Harvey owns the separate common-pool
cross-score adapter; Leibniz reviews it. No continuation outcome is implied yet.

# Latest: first accommodation scout complete; unchanged continuation admitted

12533 completed7CCO/La in107s/12CPUs/oneH200. Nine accepted movements lower the
whole composite energy2.880329kcal/mol; final physical gates pass. All final
trial multipliers exceed the0.20A source envelope; force remains, so this is a
finite boundary-limited candidate, not a minimum. Read lbt3_accommodation/
SCOUT_REPORT.md and exact workspace ROOT_REVIEW. NineMACE/eighteenGFN calls.
Leibniz stages unchanged otherthree; Harvey stages commonpool cross-scoring.
Root owns submission. All earlier molecular calls remain reused; no duplicate
scout or change to PQQ. Protein-level goal remains unmet and active.

# Latest: first intact-peptide accommodation scout submitted

Job12533 runs7CCO/La under the frozen LBT3_ACCOMMODATION_PLAN.md after the
actual derivative12521 gate passed. One GPU/12 CPUs; completion watcher1917544.
Exact inputs/configb9577bae and frozen implementationee1db81b are under
workspaces/lanm_coupled_response_20261016/lbt3_accommodation_scout_v1/.
Three independent static reviews passed; no optimization result is implied yet.
Root owns execution and review. Leibniz stages unchanged other-three adapter;
Harvey stages common-pool cross-scoring. Scout success means execution/physical
admission, irrespective of energy improvement or preference sign. No duplicate
submissions, changed physical gates or PQQ changes. Goal remains active.

# Latest: numerical force gate passed; accommodation adapter in preparation

FD12519 completed8/8 MACE and12520 completed16/16 GFN2 energies. Comparison12521
passes both steps for both metals and their paired derivative, residuals about
1e−5kcal/mol/angstrom. See lbt3_derivative/REPORT.md. No repeated calls needed.
All current owned jobs are terminal; actual completion wakes remain recorded.

The fixed LBT3_ACCOMMODATION_PLAN.md now passes its numerical prerequisite.
Leibniz implements the narrow composite-response adapter; Harvey independently
reviews state/contact/metric behavior. Root next submits7CCO/La scout, then the
other three under identical settings on execution/physical admission. This is
not yet an optimization result or a completed LanM profiler. PQQ stays intact.

# Latest execution: finite derivative matrix submitted

Whole-composite force profile12512 passed both sources and all eight geometric
finite-check inputs. MACE12519 (one GPU/12 CPUs) and standaloneGFN2 12520
(four CPUs) now execute the declared 8+16 energy cells. Completion wakes are
recorded in lbt3_fd_gpu_v1 and lbt3_fd_cpu_v1. Root owns collection/submission;
no derivative result is implied yet. Reuse all completed origins and forces.

LBT3_ACCOMMODATION_PLAN.md fixes the small own-metal response pilot before any
optimization: seven common coordinates, original physical guards, whole-peptide
composite energy/forces, twelve accepted steps maximum, shared cross-scored pool.
Leibniz implements the narrow adapter; Harvey independently reviews. Actual
search submission requires the finite derivative gate. No new DFT or PQQ change.

# Latest: whole-peptide solvent and forces complete

12500 completed all four native MACE forces with exact prior-energy replay.
12509 completed both 64-mode physical coordinate catalogues. 12510 completed
all eight standalone solvent gradients. See lbt3_solvent/REPORT.md and RESULT.json.
La selects7CCO and Dy7CCN with and without solvent; source dependence remains.
No Mex composite or biological accuracy claim follows from these endpoints.

LBT3_RESPONSE_DESIGN.md was frozen before composite-force rankings. Leibniz owns
force/profile and finite-displacement geometry staging; Harvey owns the separate
finite-energy adapter. Root owns submissions. Next is one small two-step actual
derivative check, then a bounded own-metal response search if it passes. Saved-force profile12512 completed with both source rankings and eight paired
finite-displacement inputs admitted. No optimizer or finite-displacement energy
job submitted yet. All12500/12509/12510/12512 are terminal; preserve queued completion notices and do not duplicate chemistry.

# Current checkpoint — 16 October 2026

## Latest actual execution

Saved LBT3 comparison12484 is COMPLETE after a narrow legacy-symlink pin parser
recovery; failed12478 remains preserved. See REPORT.md and lbt3_static/. Whole
scaffold barely changes fixed-source contrast but changes the compatible common
pool selection from La-source for both metals to La-source for La/heavy-source
for Dy. This is a dry conditional mechanism, not a validated protein preference.

Dy qualification12488 is COMPLETE8/8, with all four absolute energy/charge and
two solvent-component parity gates passed. New full-profile SHA635e111b is
preserved under dy_parameter_qualification/. No production registry was promoted.
The earlier not-executed paragraph below is historical.

Next declared LBT3_SOLVENT_PLAN.md is in execution: four nativeMACE forces12500
submitted, watcher1866007, exact receipt in lbt3_forces_v1/SUBMISSION.json.
Harvey owns eight standaloneGFN2 gradient adapter; Leibniz owns source-graph
coordinate catalogue. Root owns submissions and comparison. No optimizer yet.

Jacob's monomer-first steering is in effect. Mex-LanM plus the consumed LBT3
peptide supplies a small opposite-preference control; it is not two natural
LanMs or completion of the protein-level goal. See MONOMER_TARGETS.md.

LBT3 preparation12469 completed both sources and all eight task geometries:
whole256 atoms/charge−2 and mapped48-atom core/charge−1, La7CCO and Lu7CCN
model1. Declared terminal chemistry remains conditional; all deposited
coordinates were preserved. Actual preparation included two bounded missing-H
helper calls. Scheduler3s/1CPU; see the immutable preparation receipt.

Native MACE12471 completed all eight fixed evaluations,26s on12CPUs/oneH200,
312allocatedCPU-s/26allocatedGPU-s. The warm native model used physical La
multiplicity1/Dy6 and exact pinned inputs. This is an executed vacuum diagnostic,
not a solvent-corrected preference result. A saved-data comparison is being
prepared; no optimization or repeated MACE call is authorized by its receipt.

The current La-corrected nativeGFN2 profile still misassigns three Dy shell
arrays. DY_PARAMETER_READINESS.md records the exact schema evidence. A new
bounded real-core La/Dy runtime qualification is being staged separately;
the existing La-only registry does not admit Dy. No new full-profile import,
Dy parity result or solvent calculation has yet been executed in this branch.

Root owns execution/integration. Leibniz owns the saved MACE comparator; Harvey
and child own the narrow Dy staging/review. All three completed molecular/
preparation jobs12460/12469/12471 are terminal with completion wakes recorded.
Older no-submission/possible-panel paragraphs below are historical.

## Earlier recovery context

The active goal is a useful, source-robust La/Dy LanM preference profiler with
protein mechanical response and explicit occupancy/hydration. It is not yet met.
PQQ midpoint policy and scalar evidence are preserved/pushed in bb27323; the
saved Hans force diagnosis is preserved/pushed in46d91aa. No repeated chemistry.

The last discussion-only goal turn was no scientific progress. Execution resumes
with the finite PLAN.md preparation. Leibniz owns PREPARE.py; root owns review,
Slurm submission and integration. No new owned job has been submitted yet.
Existing unrelated DPANN allocations are left alone.

Jacob suggests measured monomeric proteins as the cleaner first benchmark.
Harvey is recovering a bounded primary-evidence monomer panel, preserving
reserved SpyCI-LAMBS outcomes. EXPERIMENTAL_TARGET.md already records the
consumed Hans/Mex target and the distinction between apparent protein affinity,
assembly and fixed EF3 energies. No dimer campaign is a prerequisite.

Preparation12460 is COMPLETE: four admitted geometries and eight actual full
classical cells,13seconds/1CPU, zero QM. See PREPARATION_REPORT.md and
workspaces/lanm_coupled_response_20261016/preparation_v1/execution/RECEIPT.json.
The prior no-submission sentence is historical; no owned job remains running.
The completion watcher queued its event; root reviewed actual result separately.

Next: select a real monomer comparison from actual evidence and stage its finite
response calculation. Mex plus consumed LBT3 is a possible opposite-direction
control; LBT3 is a peptide, not another natural LanM, and assay conditions differ.
Harvey owns the evidence note; Leibniz owns LBT3 source-preparation readiness.
No automatic eight-QM-cell Hans batch or dimer campaign.
