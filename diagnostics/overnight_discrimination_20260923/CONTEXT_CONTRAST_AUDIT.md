# Pocket membership mainly changes the starting solvent contrast

This read-only decomposition covers **all55** newly evaluated three-source
context/source pairs and their ten-fold precision counterparts:41 structures
from11 consumed reference proteins. It uses the actual saved MACE and native
GFN2 matrices, the existing operational selections and unchanged decision bands.
No geometry, energy, reference or label was generated or modified.

For each representation, separate the origin contrast from the change produced
by accommodation. Then take three-source minus ten-fold differences. The
absolute solvent-component change exceeds the absolute native-MACE change in
36/55 pairs, both at the origin and in the selected pool. Median absolute origin
changes are0.671MACE,1.034solvent and1.037total model kcal/mol; median absolute
change in differential accommodation is0.422. These are descriptive correlated
observations, not a significance test or attribution of a unique physical cause.

The two168-atom A8R3S4 sources responsible for lost triple decisions illustrate
the dominant mechanism:

| Source | Origin contrast shift | Selected contrast shift | Native part of selected shift | Solvent part | Change in differential accommodation |
|---|---:|---:|---:|---:|---:|
| La sample1 |+8.690|+8.422|+0.564|+7.858|−0.268|
| La sample3 |+7.703|+7.201|+0.759|+6.442|−0.502|

All values are model kcal/mol. The existing accommodation moves these smaller
contexts across the Ca band because their initial margin is reduced; it does
not newly produce a larger differential relaxation. This supports testing the
representation and numerical solvent response before restricting valid motions.

Large membership shifts also occur in correctly classified controls. MMOL2048
sample1 changes by−12.199, mostly solvent; P38539 sample4 changes by−9.856 in one
smaller context. Therefore A8 is not evidence for a universal shift toward La.
Molecular cavity/composition and electronic-solver response are not separated by
this algebra. The root cause remains conditional on the ongoing numerical test.

## Reproduction

The completed artifact is
`workspaces/context_contrast_audit_20260923/run_v1/{RESULT.json,rows.csv}`.
The script checks125 total pairs, all55 changed pairs,225 prior sources, archived
collection hashes and selected-energy closure against both comparison ledgers.
Its actual full-data run passed. There were zero new molecular evaluations.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/audit_context_contrast.py \
  --triples workspaces/union_triple_transfer_20260923/COMPARISON_v1.json \
  --tenfold workspaces/slsqp_precision_transfer_20260923/COMPARISON_v1.json \
  --output workspaces/context_contrast_audit_20260923/replay_v2
```

The output directory must be new; the script never overwrites the completed run.
