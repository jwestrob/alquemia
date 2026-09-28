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
