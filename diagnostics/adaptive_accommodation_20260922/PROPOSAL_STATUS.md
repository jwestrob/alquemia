# Adaptive angular pilot: ready, not executed by preparation agent

Eight real endpoint tasks are prepared for the frozen four consumed contexts.
Manifest SHA-256:
`b5c4218caee58458cc98ee48b0bbcf165bebb88615b08750905de8b66fb5a998`.
Workspace: `workspaces/adaptive_accommodation_20260922/proposals_v1/`.

The common four-mode choices reproduce the archived selector exactly. Source
coordinates, endpoint charges, preparation and mappings remain unchanged.
Ca/La mapping files have distinct paths but identical decoded contents; an
initial path-equality preflight was corrected before output preparation or any
scientific call. Original map pins remain separate and preserved.

The existing SLSQP constraint and analytic Jacobian enforce final maximum source
heavy displacement0.8 Angstrom, alongside angular bounds+/-0.8 rad. Intermediate
infeasible requests are retained with their actual extents and model-call status.
No penalty, coordinate clipping, new start or score selection is added. Read
[the exact plan](PROPOSAL_PLAN.md).

Seven real-fixture tests pass in **15.449 s**, zero skips
([receipt](PROPOSAL_TESTS_v1.txt)). They check all8 source maps, analytic constraint
derivatives against geometry-only finite differences, physical rejection despite
a valid angle box, inactive atoms, corrupted selection rejection, missing-candidate
accounting and the actual incomplete-pool execution gate. Independent review also
passed all4 source probes with maximum derivative error1.10e-9 Angstrom²/radian.
No optimization or molecular engine is mocked as a successful scientific result.

The real completed shared-pool pilot1209840/1209845 passes the runtime prerequisite;
its actual collection hash is
`b391d59352ff3aa2675092ad60071ef620c006f36686e473da0641adc30b75df`.
`PREFLIGHT.json`, `POOL_GATE.json` and `collection_unrun_v1.json` retain preparation
status; the latter has eight unavailable candidates and all new solvent/score
fields null. Parent owns launch and the common-pool expansion.

**New model calls/allocations by this preparation: zero.** No baseline or shared
executor changed. The native warm worker and physical Kinematics are reused.
No accuracy, unconstrained minimum, thermodynamic correction or entropy claim.

## Exact operations

From the repository, replay preparation into a fresh directory if needed:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/adaptive_angular_proposals.py prepare \
  --source workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json \
  --diagnostic workspaces/adaptive_accommodation_20260922/force_projection_v1.json \
  --agreement diagnostics/adaptive_accommodation_20260922/PROPOSAL_PLAN.md \
  --output workspaces/adaptive_accommodation_20260922/proposals_replay_v1
```

The prepared original manifest is the reviewed launch target:

```bash
ADAPTIVE_ROOT=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
sbatch --parsable --chdir="$ADAPTIVE_ROOT" \
  --output="$ADAPTIVE_ROOT/diagnostics/adaptive_accommodation_20260922/proposals_%j.out" \
  --error="$ADAPTIVE_ROOT/diagnostics/adaptive_accommodation_20260922/proposals_%j.err" \
  "$ADAPTIVE_ROOT/diagnostics/adaptive_accommodation_20260922/run_proposals.sbatch" \
  "$ADAPTIVE_ROOT/workspaces/adaptive_accommodation_20260922/proposals_v1/manifest.json" \
  "$ADAPTIVE_ROOT/workspaces/nikasha_shared_pool_20260922/pilot_v2/after_solvent_0_1209845.json"
```

The wrapper uses the immutable implementation snapshot and writes
`after_proposals_JOBID.json` even if execution fails. That collection exports actual
candidate MACE receipts and physical provenance for the separately owned common
pool. It does not launch GFN2 or select a score. Explicit dry-run/collection:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/adaptive_accommodation_20260922/proposals_v1/implementation/adaptive_angular_proposals.py \
  dry-run --manifest workspaces/adaptive_accommodation_20260922/proposals_v1/manifest.json
```
