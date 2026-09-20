# Fixed primary225 proposal transfer

Implements [the frozen continuation](FOLD_TRANSFER_PLAN.md), using the unchanged
proposal algorithm and selection expression from [the completed pilot](PROPOSAL_REPORT.md).
The production scorer is unchanged. All outputs are development on consumed
calibration-protein structural repeats, with every source retained.

The new adapter only prepares exact archived inputs and joins results. The existing
proposal engine still performs optimization, energy selection and raw collection.
Its original 30-context default remains checked separately. The explicit fold scope
requires the pinned 225 primary IDs, 208 supported mappings, 416 paired endpoints,
415 available composite origins and at most 832 new proposal solvent singlepoints.
The known missing origin prevents optimization for that endpoint. There are no
fresh q0 calculations, extra poses, failed-SCF retries or new DFT tasks.

## Preparation and validation

Run from the repository root. The output directory is immutable and must not exist.

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/accommodation_fold_proposals.py prepare \
  --design workspaces/accommodation_nonlinear_20260920/fold_maps_v1/design.json \
  --comparison workspaces/accommodation_goal_20260920/folds_v1/comparison_v1.json \
  --reference workspaces/accommodation_nonlinear_20260920/proposal_calibration_v1/REFERENCE.json \
  --agreement diagnostics/accommodation_nonlinear_20260920/FOLD_TRANSFER_PLAN.md \
  --output workspaces/accommodation_nonlinear_20260920/fold_proposals_v1 \
  --cpu-python /groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  --gpu-python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/mace_hybrid_20260916/software_v1/venv/bin/python

/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python scripts/accommodation_proposals.py validate \
  --manifest workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json
```

The manifest pins its source files, model, software, primary numerical recipe,
fixed new reference and original bands. `Q0_AUDIT.json` records every availability
and failure. Origin coordinate replay admits only the existing 1e-12 Å mapping
roundoff tolerance; actual archived energy coordinates remain exact pins.

## Execution

The existing `run_proposals.sbatch` takes the absolute proposal manifest path.
It uses one H200, 32 CPUs and 200000 MiB host RAM, performs all finite MACE starts,
and creates actual-input primary-GFN2 tasks only for successful proposals. It
automatically collects even if the proposal executor exits nonzero.

The resulting `proposal_GFN2/manifest.json` is an index, with **zero executable
tasks**. Its four `shard_0` through `shard_3` manifests have disjoint, contained
input/output paths and independent executor locks. Each has a real `PREFLIGHT.json`.
Submit `run_fold_proposal_gfn2.sbatch` with the absolute main proposal manifest and
one explicit shard index. Each job uses 64 CPUs, eight concurrent eight-rank
endpoints and 128 GiB. Never execute the master as a solver task manifest.
Each wrapper collects partial/failed statuses regardless of solver exit code;
the first failure is preserved and never replaced with an origin success.

## Final collection and fixed comparison

After all four solvent jobs terminate, this command collects without rerunning chemistry:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/implementation/accommodation_proposals.py collect \
  --manifest workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/manifest.json \
  --output workspaces/accommodation_nonlinear_20260920/fold_proposals_v1/final_selection.json
```

`accommodation_fold_proposals.py compare` requires explicit `--manifest`,
`--selection`, `--dft` and `--output` paths. The DFT path must be an actual saved
collection from the preserved fold calculation; no implicit live polling is done.
It produces JSON and Markdown with all225, La100 and Ca125 outcomes; complete
La4/Ca5/equal-mean-of-medians descriptors; every100 three-of-four La subset;
old/new proposal bands; and unchanged native-core/context, composite and DFT
comparators. Pairwise complete-case summaries prevent missing proposals from
making an arm appear more accurate. Original scores remain present even when
the corresponding proposal is unavailable. No calibration code runs in collection.

## Verification

`tests/test_accommodation_fold_proposals.py` uses actual pinned preparations,
archived energies and a staging-only copy of actual solver inputs. It checks
population limits, origin reuse, all available original contrast algebra,
prelaunch failure denominators, unchanged original30 results, disjoint shard
containment and actual receipt lookup. It invents no successful scientific output.
Scientific integration results and receipts are reported separately after execution.
