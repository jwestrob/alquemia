# Four-cell expanded-region electronic-component harness

`scripts/metal_environment_partition.py` prepares, validates, executes and collects
exactly Ca_A, Ca_B, La_A, La_B using each endpoint's own 63-atom XYZ. B therefore
uses its moved QM HG1, unlike the original fixed-core six-cell adapter. The
external 9,078-charge A/B fields are identical. Inputs are original native
r2SCAN-3c DefGrid3 TightSCF EnGrad with DoEQ false. The named resource layout is
four workers ×86MPI ranks; the existing allocation-aware renderer derives memory.
No second runner or six-cell adapter is introduced.

Collection reports each metal's B−A, La response minus Ca response, and that
contrast's change from the actual pinned small-region collection (approximately
+0.352795817 kcal/mol). It never subtracts absolute energies between regions.
The original source energies and receipts are revalidated, not reconstructed
from the rounded number above. Partial/invalid matrices return unavailable.

HG1 gradients are also projected onto the actual increasing-angle arc tangent,
in kcal/mol/Å and kcal/mol/radian; force is the negative gradient. Original MM
index2412 and expanded QM index57 represent the same physical hydrogen motion.
They remain derivatives of different approximate electronic-component energies,
not full hybrid forces. No force admission tolerance or biological gate is added.

Every collection explicitly retains numerical qualification as unestablished
because the original rigid-motion gate failed. Coarse partition sensitivity may
be informative; small differences remain unresolved. Complete QM/MM, affinity,
classification, dynamics and optimization remain unavailable.

Five real-artifact unit tests pass: paired endpoint geometry/field, intentionally
corrupted B coordinates, corrupted state, real archived response algebra with
missing-cell nulls, and four-cell preparation/dry-run plus missing-task rejection.
No synthetic molecular results or scientific executable were used in these tests.
Root owns submission; this implementation task submitted no jobs.

## Exact preparation command

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/metal_environment_partition.py prepare \
  --inputs workspaces/metal_environment_response_20260926/preparation/expanded_region_v1/INPUTS.json \
  --source-collection workspaces/metal_environment_response_20260926/reference_scout_v1/FINAL_COLLECTION_1219207.json \
  --agreement diagnostics/metal_environment_response_20260926/PARTITION_DIAGNOSTIC_PLAN.md \
  --output workspaces/metal_environment_response_20260926/partition_v1 \
  --workers 4 --mpi-ranks 86
```

Then use the manifest's frozen implementation for `dry-run`, `execute` or
`collect`, passing `--manifest` and an optional new `--output` path. Preparation
and collection outputs refuse overwriting. The independently frozen plan governs
interpretation; numerical failures are not relaxed by the harness.
