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

