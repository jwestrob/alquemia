# Intact-chain OMOL operations

Read INTACT_CHAIN_PLAN.md for the fixed inputs and interpretation. Production
is unchanged. Energy-only output explicitly has no forces. All task manifests
use the existing immutable runner, receipts and typed cache validation.

## Core equivalence bridge

The four-task manifest was prepared and passed the copied runner dry-run.
Job1200807 uses `intact_core_v1`; its exact sbatch argv is saved in
`intact_core_v1/submission.json`. These commands inspect without new inference:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_core_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/intact_core_v1/manifest.json
python -m unittest discover -s tests -p test_mace_omol_intact.py -v
```

## Intact numerical qualification

After the actual core collection passes, prepare the14-task ALPHA_1F6S gate:

```bash
python scripts/mace_omol_intact.py --collection workspaces/mace_omol_20260917/benchmark_v1/collection_job_1200799.json --preparation workspaces/mace_global_benchmark_20260916/prepared_v1/preparation_manifest.json --agreement diagnostics/mace_omol_20260917/INTACT_CHAIN_PLAN.md --core-qualification workspaces/mace_omol_20260917/intact_core_v1/collection_job_1200807.json --output workspaces/mace_omol_20260917/intact_qualification_v1
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_qualification_v1/implementation/mace_hybrid.py dry-run --manifest workspaces/mace_omol_20260917/intact_qualification_v1/manifest.json
```

Preparation refuses overwrite. Do not repeat preparation for an existing
manifest; inspect its receipt and live job before recovery. Core bridge success
is an engineering check, not predictive validation. The next16-task comparison
is conditional on the14-task numerical gate, irrespective of score direction.

Qualification manifest SHA256:
`b1609f0b81c73592881102ca123bce39f78d529ba23b64373d26d9e83d491488`.
Preparation ran the same validator as dry-run and passed. Job1200808 executes
this manifest; exact allocation/submission arguments are in its `submission.json`.
The executor repeats validation before inference and collects on exit.

```bash
squeue -j 1200808 -o '%.12i %.24j %.10T %.10M %.24R'
python scripts/mace_hybrid.py collect --manifest workspaces/mace_omol_20260917/intact_qualification_v1/manifest.json --output workspaces/mace_omol_20260917/intact_qualification_manual_collection_v1.json
```

Manual collection computes no new scientific outputs. It preserves missing
endpoints as unavailable; it cannot make an incomplete numerical gate pass.
