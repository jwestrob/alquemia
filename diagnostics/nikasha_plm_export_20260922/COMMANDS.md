# Existing PLM export commands

Run from the Nikasha repository root. Every data path is supplied by the pinned
[input inventory](INPUTS.json); no source editing or hidden environment data path
is required.

## Read the completed export

```bash
head -n 3 workspaces/nikasha_plm_export_20260922/export_v1/proteins.tsv
python -m json.tool workspaces/nikasha_plm_export_20260922/export_v1/EXPORT.json
```

The first table has all 176 outcomes, including missing scores. `EXPORT.json`
defines all 14 relational tables, source hashes, original DFT method/bands and
distinct transcript measures.

## Reproduce to a new directory

```bash
export NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$NIKASHA_PY" scripts/nikasha_plm_export.py \
  --inventory diagnostics/nikasha_plm_export_20260922/INPUTS.json \
  --output workspaces/nikasha_plm_export_20260922/export_recheck_v1
```

This reads existing outputs only. It performs no scoring, normalization,
correlation, phylogenetic or biological analysis. An existing destination or
changed source hash fails explicitly; choose a new destination for another
export. The original `export_v1` remains intact.

## Real-fixture checks

The normal Python environment already has pytest; the scientific CPU environment
does not. No new package installation is needed.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python -m pytest -q tests/test_nikasha_plm_export.py
```

Tests use the actual pinned PLM source tables and explicitly corrupted copies for
failure cases. They do not fabricate scientific results.
