# MACE proposal / composite selection operations

Declared scope: [PROPOSAL_PLAN.md](PROPOSAL_PLAN.md). Baseline/default unchanged.
All commands run from the repository root. Existing scientific executors and
the separate full-composite pilot are independent of this experiment.

## Prepared input and preflight

The immutable 30-context/60-start manifest is
`workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json`.
Its `Q0_AUDIT.json` verifies all 120 archived primary GFN2 outputs, recipes,
states and actual executed coordinates; all 60 endpoint origins are available.
Known labels appear only in subsequent reporting. The 28 known reference origin
contrasts reproduce the released composite results exactly using archived values.

The evaluated q0 XYZ is pinned for both PLM cases. For 07ab only, it differs
from the torsion design's origin by at most 2.22e-16 Angstrom. Both source pins
and this difference are retained under the established 1e-12 Angstrom mapping
tolerance. No q0 calculation or coordinate rewrite is performed.

Fresh MACE q0 replay uses the existing torsion-control tolerance of
0.01 kcal/mol; actual differences remain in endpoint receipts. This check is
stricter than the declared 0.10 kcal/mol composite-selection scale and does not
change the optimizer, source set or energy-selection rule.

Six real-artifact tests pass, without molecular calls or skips:
[PROPOSAL_TESTS_v1.txt](PROPOSAL_TESTS_v1.txt). Both scheduler wrappers pass
`bash -n`. Preflight checks actual archived selection algebra, method/state,
missing-data handling and complete input coverage; GPU execution is a separate
scientific operation.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/accommodation_nonlinear_20260920/proposals_v1/implementation/accommodation_proposals.py \
  validate --manifest "$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json"
```

## Execute the declared experiment

Run each submission once; consult saved submission receipts before repeating.
First, create proposals on one warm H200 allocation:

```bash
sbatch --parsable \
  --output="$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/proposals_%j.out" \
  --error="$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/proposals_%j.err" \
  diagnostics/accommodation_nonlinear_20260920/run_proposals.sbatch \
  "$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json"
```

After successful proposal execution, the script creates the actual finite
`proposal_GFN2/manifest.json` and its runner dry-run receipt. Failed proposals
are explicit unavailable entries and create no solvent task. An exactly unchanged
proposal can reuse its compatible q0 outputs. At most 120 new GFN2 calls exist.

```bash
sbatch --parsable \
  --output="$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/solvent_%j.out" \
  --error="$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/solvent_%j.err" \
  diagnostics/accommodation_nonlinear_20260920/run_proposal_gfn2.sbatch \
  "$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json"
```

Each wrapper collects completed/failed status even if its executor exits nonzero,
preserving the executor's exit code. Every model evaluation saves actual energy,
forces, geometry, active gradient and cache identity. Endpoint results retain
optimizer success, raw derivatives, bounds and physical checks. Solver receipts
and failed artifacts stay visible; an unavailable proposal never masquerades as
a selected origin.

## Collect and interpret

The CPU wrapper writes `final_<job>.json` and Markdown. For a separate immutable
collection snapshot:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  workspaces/accommodation_nonlinear_20260920/proposals_v1/implementation/accommodation_proposals.py \
  collect --manifest "$PWD/workspaces/accommodation_nonlinear_20260920/proposals_v1/manifest.json" \
  --output workspaces/accommodation_nonlinear_20260920/proposals_v1/separate_collection_v1.json
```

Report all 25 canonical references, three consumed crystals and two unlabeled
PLM contexts separately. Selected geometry minimizes only the actual two-member
candidate set, with the frozen 0.10 kcal/mol decrease requirement. No continuous
composite minimum, entropy, binding-free-energy or PLM accuracy claim follows.
Old bands remain an explicit developmental transfer check. Native validation or
broader fold transfer requires a separately recorded manifest.
