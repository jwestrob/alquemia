# Scaffold inventory commands

Read the completed result (no calculations):

```bash
cat diagnostics/scaffold_environment_20260922/REPORT.md
cat diagnostics/scaffold_environment_20260922/RESULT.json
```

Reproduce source/parameter/term inventory only, in a new explicit output directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/scaffold_environment_inventory.py \
  --config diagnostics/scaffold_environment_20260922/INPUTS.json \
  --output workspaces/scaffold_environment_20260922/independent_inventory_v1
```

The command rejects existing output directories, unknown source chemistry, changed
input hashes, mismatched paired maps and missing real cap bonds. It instantiates
and serializes protein parameters; it creates no OpenMM Context and calculates
no energy or force. Large atom/term/XML inventories remain under `workspaces/`.

Actual-fixture checks on the preserved successful inventory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  -m unittest discover -s tests -p test_scaffold_environment_inventory.py -v
```

No energy/force executor or new scheduler submission is part of this branch.
The optional nine-point protein diagnostic is a separate proposal in REPORT.md.
