# Restart: Dy electronic representation diagnosis, 28 September

Read FROZEN_F_ASSESSMENT.md and EXPLICIT_F_DIAGNOSIS.md in this directory.
The saved-output audit and legacy frozen-f review are complete. Recommendation:
one separately named frozen-f Dy analytic-force qualification on the exact consumed
50-atom Hans EF3 fixture, before environment or protein relaxation expansion.
No new molecular job is active in this session. PQQ jobs belong to another session.

Physical Dy(III) is 4f9, S=5/2; ECP55 has a spin-free restricted valence model.
For this50atom charge−1 fixture,263 all-electron count minus55 core =208 explicit
electrons,104alpha/104beta. This is not the physical sextet solved as a singlet.
Pin actual lcecp-1-TZVP ECP/orbital/AutoAux assets from FROZEN_F_INVENTORY.json.
Use a new research adapter; existing STATES Dy sextet/ECP28 must remain intact.
Do not silently reuse native-3c offsets or import old SP outputs as forces.

Next implement preparation/state-aware parsing and a finite manifest for one
analytic-gradient origin plus source-derived donor displacement checks. Freeze
force tolerances, atom mappings and small displacement before energies. First
execute the origin scout, inspect actual native ECP/electron/gradient evidence,
then run dependent directional checks only if the origin is admissible. No old
campaign restart, whole-protein optimization, large matrix or reserved labels.

Publication/checkpoint status: origin/main successfully updated through776a13b.
Local diagnosis commits8be9f59 and3ebaacf exist; subsequent pushes failed twice
because github.com DNS resolution temporarily failed. Retry normal authenticated
`git push origin HEAD:main`; no force push or history rewrite. The local raw
workspaces remain on shared storage, not in GitHub. Vault updated separately.

## Update: actual running scout and successful remote backup

28 September03:06PDT: worker1220294 confirmed running40CPUs, collector1220295
pending afterany. Frozen-f molecular SCF advancing; no admitted endpoint yet.
All changes throughb3eb8c1 successfully pushed to origin/main; previous DNS
failure is resolved. No remote raw-workspace backup is implied.

The dependent force check is implemented, three real-fixture checks pass.
It explicitly refuses preparation until the origin energy/gradient is complete.
After origin collection, from repository root:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/metal_environment_frozen_f_response.py prepare --origin-manifest workspaces/metal_environment_response_20260926/dy_frozen_f_scout_v1/manifest.json --plan diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/FROZEN_F_SCOUT_PLAN.md --ranks 10 --workers 4 --output workspaces/metal_environment_response_20260926/dy_frozen_f_direction_v1
```

Inspect current free slots before selecting execution layout;40total here is a
recorded runnable layout, not a requirement to reserve a larger idle node. Use
existing batch/collector templates with the new response script name and matching
manifest CPU slots. Arm both actual terminal and live-health watches as in scout
SUBMISSION.json; use new job IDs, receipt paths and event names. Never rerun the
origin. No dependent calculations have yet been submitted.

Embedding agent independently recovers consumed-source transfer design/occupancy;
root owns all molecular execution. Source/label audit does not block collecting
the active scout. Overall goal remains La/Dy discrimination, not merely forces.

## Next matched exchange implementation — pushed fc7ecd3

`metal_environment_lady_compact.py` is ready and refuses preparation unless the
actual four-direction collection reports its declared force gate passed. Six
real-fixture tests pass, including all three actual paired source coordinates and
electron counts. No new exchange manifest or molecular jobs have been submitted.
Current task remains1220300/collector1220301, confirmed running at03:18PDT.

After successful FINAL_COLLECTION.json, prepare (from repository root):

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/metal_environment_lady_compact.py prepare --origin-manifest workspaces/metal_environment_response_20260926/dy_frozen_f_scout_v1/manifest.json --force-gate workspaces/metal_environment_response_20260926/dy_frozen_f_direction_v1/FINAL_COLLECTION.json --la-basis legacy/qmmm_lc/ecp_lib/orca_La.ecp_basis --la-aux diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/lcecp1_tzvp_La_autoauxj_orca611.inc --plan diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/COMPACT_EXCHANGE_PLAN.md --ranks 8 --workers 5 --output workspaces/metal_environment_response_20260926/lady_compact_exchange_v1
```

This is a40CPU finite five-new-cell layout; assess actual available slots before
submission. Reuse of the real Hans8DQ2_Dy origin is mandatory. No native3c La
substitution or importing old model scores. Report both Hans-vs-Mex contrasts;
unknown source/state populations and missing environmental/scaffold physics stay
explicit. Use the shared existing batch templates with this script name, preserving
separate workspace/wake IDs. If force gate fails, do not bypass it.

## Embedded execution bridge prepared, not executed

Root added metal_environment_lady_embedded.py. It collects actual state/ECP/gradient
plus required point-charge gradient and physical donor torsion projection. Execution
requires the passed local force gate and an executor pinned in the prepared snapshot.
The agent's existing preparation-only Hans manifest predates this bridge and cannot
be executed by it; prepare a NEW versioned workspace after the gate using the same
frozen source/assets, rather than editing its immutable snapshot. The old manifest
still passes read-only dry-run. No embedded molecular cells have been submitted.
The shared parser now reads charge from declared state, enabling the real Mexcharge0
region; compact charge−1 tests remain unchanged. Actual embedded analytic integration
is unrun and must not be claimed from preparation tests.

Current active work remains directional1220300/collector1220301. Latest observation
at03:27PDT:13macroSCF iterations on plus, residual improving; no failure alerts.

## Two scientific matrices running — 28 September03:51PDT

Compact comparison1220308 remains RUNNING40sharedCPUs, collector1220309. Five
new La/Dy origins plus exact reused Dy origin; report both Hans-source exchanges.

First full-region frozen-f Hans scout1220312 confirmed RUNNING on node-344-8t-1;
collector1220313 afterany. Four86MPI endpoints assign all344slots; actual runtime
receipts confirm full8256990MiB scheduler RAM policy,mem0,25%operationalheadroom.
Workspace lady_frozen_embedded_hans_v2; manifest708ba36d326f1435f1fe6d2a09e455fce0d52b2a6083bab9e1dc3d7e31627925.
Exact terminal/live-health commands/PIDs in SUBMISSION.json. No other full-region
source launched. Read electronic_diagnosis_20260928/EMBEDDED_EXECUTION_PLAN.md.

This is the fixed195atom source and frozen+2degreeAsp85 perturbation, same finite
field for La/Dy andA/B, targetcharge−1/806effectiveelectrons, physicalLa1/Dy6 with
restrictedvalence1. Parser requires actual customECP, native gradient andpcgrad,
and projects the physical torsion with existing cap Jacobians. Electronic component
only; no complete solvent/MMmechanics or affinity claim. Monitor SCFprogress and
partialcompletion; do not repeat unmonitored largeDy expenditure. PQQ unchanged.
On terminal wake collect actual per-metal work and paired response, with missing
cells explicit. A converged local origin alone does not qualify embedded forces.



## Matched C-H repair comparison launched — 28 September04:23PDT

Worker1220316, collector1220317, six new analytic endpoints,48CPU standard exclusive,
mem0,6×8MPI; no old energy reuse. Workspace lady_compact_CHrepair_v1. New protocol
nikasha_LaDy_compact_CHrepair_v1, exact reviewed repaired-input manifest compact_CH_repair_v1.
Terminal and futurev5 live-health watchers armed; commands/PIDs in SUBMISSION.json.
Read electronic_diagnosis_20260928/COMPACT_CH_REPAIR_PLAN.md and COMPACT_CH_TRANSFER_REPORT.md.
Only14/14/10 real carbon-boundH atoms changed per source, identically for both metals.
Old compact protocol collection still reproduces exact published completed values.
Original caps/exchangeableH/waters/heavies remain exact; this is a targeted repair.

All3 fullsource C-H geometry candidates separately admitted from saved forces/coords;
original serialization errors, missing optimizer flags and CH2-label volume flags
remain visible. No unrecorded repeat. Repair itself32CPU-s; electroniccostpending.
Embedded1220312 still running original normalized-H195atom Hans; collector1220313.
No new embeddedsource, production change or classifier-success claim. Root owns
execution/results. Other agents' preparation/audits completed. Inspect actual jobs
and receipts after restart; do not repeat either matrix or its successful origins.



## Repair comparison closed; coupled-response work continues — 28 September04:58PDT

Repaired compact1220316/collector1220317 completed all6 endpoints. D(Hans8DQ2,Mex)
=+11.090948176 and D(Hans8FNR,Mex)=-18.612513541kcal/mol. Changes-.239938169/-.144330681;
source reversal persists. Read electronic_diagnosis_20260928/COMPACT_CH_COMPARISON.md
and force results. Large absolute H-repair works cancel between metals; do not call
this a discrimination gain. Close compact C-H rescoring, retain improved preparation.
Actual electronicworker+collector85,106CPU-s; separate geometryprep32CPU-s.

Original embedded1220312 remains running on oldnormalizedH;1220313 collector.
La SOSCF residuals improve; Dy remains oscillatory. Saved partialoutput evidence
under embedded_SOSCF_observation_v1 is NOT accepted molecular energies.
One diagnosed nextscout1220323 runs repairedDy_A with explicitTRAH,112MPI on
node-112-1500g-1,mem0; collector1220324. Workspace lady_repaired_trah_scout_v1;
manifest aaf27076a6e0976263a28a69f9b9c1e1151d90ce123623819b37fb8b8735049d.
Newsolver AND repairedH/field differ fromold; do not claim uniquecause or reuseold
energies. Read EMBEDDED_TRAH_SCOUT_PLAN.md; exactwatches in SUBMISSION.json.

Repaired full-region source contracts for all3sources are ready under
lanm_ef3_CboundH_repaired_v2; other repaired electroniccells unsubmitted.
Selected classical scaffold ledger is tested separately; coarse andhalfstep
MM-LJ finite-difference gates fail but4.0018errorratio and Richardson~3.18e-5kcal/rad
indicate finite-step truncation. No fullhybrid force qualification claimed.
Embedding_preparation owns finalboundedreport; no additionalQM/optimization.
Root ownsallscoring/execution. PQQ and sharedproduction writer unchanged.

Email update accepted bylocalrelay04:36PDT; recipient delivery not independently
verified. Vault/report checkpoints current. Goal remains active: neither a solver
capability nor repairedgeometry is a working La/Dy discriminator.


## Embedded Hans comparison complete — 28 September05:30PDT

All four original-H frozen-f embedded endpoints1220312 and collector1220313 are complete. La work−0.04353110, Dy+0.06425619, differential+0.10778728kcal/mol for the prescribed +2degree Asp85 motion. Native gradients/point-charge gradients collected; no affinity or full force qualification. Actual cost1,828,018allocatedCPU-s,0GPU. Read electronic_diagnosis_20260928/EMBEDDED_HANS_RESULT.md. Do not rerun this matrix.

Repaired Dy TRAH scout1220323 remains live, collector1220324 pending. At last check42min, macro2 with improving electronic gradient; no accepted repaired energy yet. Finite physical force assembly now works on archived native La A/B plus matching classical terms (1220328,1CPU-s). Agent finite_field_derivatives will apply it to the completed frozen-f endpoints, reusing exact La classical data and computing only missing Dy classical terms. No new quantum stage or optimization launched. Source reversal still unresolved; goal open; PQQ unchanged.


## Repaired origins running; Glu91 matrix prepared — 28 September06:05PDT

Dy TRAH1220323/collector1220324 remains live, no accepted repaired Dy endpoint yet. Missing repaired La origin now runs independently as1220332/collector1220333,112MPI on standard node-112-1500g-2,mem0; source-matched old La orbitals initialize a NEW calculation, never fill its energy. This explicitly revises the earlier sequential scheduling gate; see REPAIRED_LA_SEEDED_PLAN.md. All displaced Glu91 electronic cells still wait for the repaired origins/forces. Both jobs have terminal and health watchers; exact commands/PIDs in each SUBMISSION.json.

La workspace lady_repaired_la_seeded_v1; the unseeded lady_repaired_la_origin_prepared_v1 remains unsubmitted. Glu91 source pair glu91_motion_v1 and four-cell lady_glu91_response_prepared_v1 pass preparation/dry-run, no submissions. Mode was selected from actual original-H differential forces: Glu91 normalized load~2.75 versus~.066 for tested N83 hinge. Do not generalize to all scaffold modes or biological affinity. Field/caps unchanged for the five moved real donor atoms; inventory/covalent checks pass.

Classical quarter-step passes unchanged thresholds after diagnosed truncation; original coarse/half failures and final stdout-only process failure remain preserved. Force assembly now includes all four original frozen-f endpoints and exact classical partners. No LanM accuracy improvement or production change yet. Orbital integration22real-artifact tests pass; actual seeded execution/collection still needs qualification. Read force_assembly_20260928/FROZEN_F_PROJECTION_REPORT.md and coupled_scaffold_20260928/QUARTER_STEP_REPORT.md.


## Direct repaired Hans/Mex comparison launched — 28 September06:13PDT

All six repaired origins now run: Hans8DQ2 Dy1220323/collector1220324 and seededLa1220332/collector1220333; Hans8FNR pair1220334/collector1220335; Mex8FNS pair1220336/collector1220337. Two new pairs use separate idle224CPU high-memory nodes,2×112MPI each,mem0, normal priority. Their terminal and health watchers are armed; exactcommands/PIDs in SUBMISSION.json. No Glu91 displaced endpoints launched.

Read electronic_diagnosis_20260928/REPAIRED_EMBEDDED_EXCHANGE_PLAN.md for the predeclared direct representation comparison and numerical-policy differences. Electronic-only conditional EF3 exchanges are not whole-protein affinity labels. Both Hans sources remain separate with source waters/spectators preserved. The joinable6cell collector and exact runnable command are in repaired_exchange_20260928/REPORT.md; three real-archive tests pass. First table0/6, no values invented. Do not rerun completed compact/original-H tests or consume reserved labels. PQQ unchanged; goal still open.


## FNR numerical divergence and explicit recovery — 28 September06:34PDT

Both Hans8FNR default/PModel solves diverged together (hundreds ofHartree energy rises, density changes>1000); no energy accepted. Stopped owned1220334 after1244s; collector1220335 completed2s. Cost278,658allocatedCPU-s. Inputs/outputs/primary collection remain unchanged and unavailable. Snapshots and checks are in electronic_diagnosis_20260928/FNR_DEFAULT_DIVERGENCE.json; no obvious atom/field coincidence or below-threshold overlap eigenvalue found. Unique cause unresolved.

Same geometry/state/Hamiltonian, explicitTRAH two-origin recovery1220342/collector1220343 submitted on released224CPU node,2×112MPI,mem0. New workspace lady_repaired_origins_Hans8FNR_trah_v1; watchers armed via SUBMISSION.json. No failed wavefunction reuse or iteration/threshold increase. Other jobs1220323(DQ2Dy),1220332(DQ2La),1220336(Mex pair) remain unchanged. Use separately named COMPARISON_RECOVERY_v1.json for recovery; primary configs still point to original failedFNR. Do not hide fallback or rerun the primary pair. Glu91 remains unsubmitted. Read FNR_TRAH_RECOVERY_PLAN.md.

All six classical bridges completed8CPU-s total; separate electronic and finite sums are available once actual quantum results finish. No demonstrated LanM discrimination yet. Current runner/seed tests23pass; comparison tests4pass. PQQ unchanged.


### Archived electronic-population diagnostic, Sep28
Four original-H completed states conserve charge and retain positive printed frontier gaps. Population schemes disagree even on the sign of the metal-charge difference; no oxidation-state or repaired-state failure diagnosis follows. No new molecular calls. See electronic_diagnosis_20260928/POPULATION_REPORT.md. Four repaired quantum workers still live; no accepted repaired comparison yet.

## First repaired full-region endpoint accepted — 28 September

Hans8DQ2La1220332/collector1220333 completed and passed collection:−5728.001479272300Eh,25SCF cycles,580,500allocatedCPU-s including collector. Physical gradients assembled with exact repaired classical partner. Glu91 total angular gradient decreases3.46524→1.35121kcal/mol/rad; this is one-metal repair response, not a La/Dy result. See force_assembly_20260928/REPAIRED_LA_REPORT.md and pinned projection JSON. Dy origin1220323, Mex pair1220336 and FNR recovery1220342 remain live; summary1220344 waits. No Glu91 quantum submission, no default change.

## Mex repaired pair accepted — 28 September

Worker1220336/collector1220337 completed La/Dy, energies−5636.761993646614/−5642.249821148190Eh. Physical gradients assembled with exact classical partners. Cost1,344,450allocatedCPU-s; no GPU. Comparison3/6accepted, Hans contrasts still null. Read repaired_exchange_20260928/MEX_REPORT.md and MEX_RESULT.json. Raw cross-element contrast is not affinity. Only quantumworkers1220323(DQ2Dy) and1220342(FNRLa/Dy) remain live;1220344 summary still pending. No Glu91 endpoint launched.

## Repaired Hans8DQ2 pair accepted — September28 08:39PDT

Dy1220323/collector1220324 completed normally with native gradients and field gradients. Cost1,501,026allocatedCPU-s (13402×112+2),zeroGPU. Accepted comparison now4/6: electronic Hans8DQ2-minus-Mex exchange+11.6984892964kcal/mol, finite electronic/classical+11.6898571756. This retains conditional relative La preference for this source; FNR remains null, so the earlier source reversal is not resolved. No affinity or broad discrimination claim.

Repaired physical Glu91 gradients: La1.3512088726,Dy11.3711440764kcal/mol/radian; differential10.0199352038,normalized2.7169243037kcal/mol/angstrom. The declared small donor-motion test remains justified after repair. No displaced quantum calls launched yet. Next: prepare same-metal admitted-origin seeds for the already-declared ±1degree four-cell test, check time before10PDT cutoff, and execute if feasible. FNR recovery1220342 and collector1220343 remain live/pending; final summary1220344 waits. Do not repeat completed origins.

Actual comparison: workspaces/metal_environment_response_20260926/repaired_exchange_DQ2_complete_v1.json. Force result: force_assembly_20260928/REPAIRED_PAIR_PROJECTION.json, generated by PROJECT_REPAIRED_PAIR.py using exact paired coordinates/source IDs. Complete hybrid derivative qualification and solvent remain absent.

## Glu91 response launched — September28 08:42PDT

Worker1220359/collector1220360 executes the four prescribed ±1degree Glu91 endpoints from matched accepted repaired La/Dy orbital seeds. Workspace lady_glu91_seeded_v1; manifest97c426e3c490a34ee4d881c85e5147dc05aa390cab4615b8b19466705c400c88. Four86MPI tasks,344CPU high-memory node,exclusive mem0,normal priority. No origin repeats. Completion/health watchers armed; exactcommands/PIDs in SUBMISSION.json. Source geometry and states unchanged from glu91_motion_v1. Preparation dry-run passed. Read GLU91_SEEDED_EXECUTION_PLAN.md for purpose and limitations.

Next collect actual energies/gradients and match classical terms before comparing the common origin/minus/plus pool. Do not interpret candidates as populations or a minimum; compare energy derivatives with repaired-origin projected gradients. FNR1220342 still runs separately;1220344 aggregates repaired origins only, not this new motion. Preserve partial results if shutdown interrupts; do not blindly relaunch after restart. PQQ unchanged.

## Glu91 collection and classical bridge — September28 08:48PDT

Classical four-state job1220361 completed3CPU-s/16 component queries; no quantum or origin repeats. Its B−A response+0.033739411802kcal/mol is identical for both metals. Full bridges: workspaces/metal_environment_response_20260926/glu91_classical_20260928. Agent report/tests committed79362d9; no active agent jobs.

Root common-pool collector glu91_response_20260928/COLLECT.py/CONFIG.json executed against real accepted origins: four displaced states correctly unavailable, no pool fabricated. Summary1220362 waits afterany1220360, then writes workspaces/metal_environment_response_20260926/glu91_response_final_v1. Completion watcher details in glu91_response_20260928/SUBMISSION.json. No chemistry from this collector. Numerical analysis policy GLU91_ANALYSIS_POLICY.md frozen before displaced results. Current quantumworkers1220342FNR and1220359Glu91 remain live. GitHubDNS intermittently prevents pushes; inspect git origin before assuming all latest commits remote.

## Glu91 finite response complete — September28 09:16PDT

All4seeded endpoints1220359/collector1220360 and common-pool summary1220362 completed. La selectsorigin,Dy−1degree; finite Dy−La poolchange−0.1739401kcal/mol. Directional derivative residuals0.00948La/0.00859Dy below frozen0.1kcal/mol/angstrom diagnostic tolerance. Different donor response is supported, but the small finite shift does not resolve the source reversal or prove improved discrimination. Dy remains nonstationary at the selected boundary; do not extrapolate harmonic correction or automatically widen search. Read glu91_response_20260928/REPORT.md and RESULT.json. Exact incremental cost697,991allocatedCPU-s including16classicalqueries,zeroGPU. No quantum repeats.

Only quantumworker1220342Hans8FNR remains live, collector1220343 and original-pair summary1220344 pending. Completed donor work requires no restart; its wake alerts may arrive later. No new molecular stage planned before interpreting FNR. PQQ unchanged; goal remains open. Current all-six donor result is workspaces/metal_environment_response_20260926/glu91_response_final_v1/RESULT.json.
