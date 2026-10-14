# Actual crystal preparation bound; scoring unrun at freeze

The actual Slurm4093 preparation admits chain A and retains chain B as
unsupported: its mapped catalytic Asp partner is 3.511 angstrom away, beyond
the unchanged 3.5 angstrom gate. No cutoff change or source rescue is included.

Chain A has the exact Q9 local state key
`dbdead29335579d38decb6e1b85ae542a7b16390e8f98e90fb01225309aca561`.
The separately validated remote terminal completion is pinned in its actual
admission. The chemical interpretation remains conditional oxidized PQQ;
crystal redox has not been established.

Frozen native config:
`workspaces/exaa_source_packing_20261014/crystal_native_v1/CONFIG.json`

SHA256: `b35eb5a75698fd574ba3cda8d6cc7474699bbad824da66db63a334787951cf73`.

Root execution from the repository root:

```bash
sbatch \
 --output=workspaces/exaa_source_packing_20261014/crystal_native_v1/run_%j.out \
 --error=workspaces/exaa_source_packing_20261014/crystal_native_v1/run_%j.err \
 diagnostics/exaa_source_packing_20261014/crystal_scoring/RUN.sbatch \
 "$PWD/workspaces/exaa_source_packing_20261014/crystal_native_v1/CONFIG.json"
```

This performs the existing Slurm physical preflight and actual two-fixture
qualification before any new source evaluation. Two unchanged optimizer starts
are planned for chain A; chain B receives no molecular calls. The two-chain
denominator remains explicit. No independent model test or search has executed
as part of metadata binding.

Root may queue the unchanged-solvent continuation after the actual native job:

```bash
sbatch --dependency=afterany:NATIVE_JOB_ID \
 --output=workspaces/exaa_source_packing_20261014/crystal_solvent_%j.out \
 --error=workspaces/exaa_source_packing_20261014/crystal_solvent_%j.err \
 diagnostics/exaa_source_packing_20261014/crystal_scoring/SOLVENT.sbatch \
 "$PWD/workspaces/exaa_source_packing_20261014/crystal_native_v1/MATRIX_INPUT.json" \
 "$PWD/workspaces/exaa_source_packing_20261014/crystal_solvent_v1"
```

Leave the new solvent workspace absent until its Slurm job creates it from the
actual native export. Both stages request seven CPUs, one node, standard queue;
neither specifies memory or time overrides. Root owns completion wakes and
actual scientific interpretation. An absent or interrupted native result does
not trigger a new model run.

Frozen code SHA256:

- `NATIVE.py`: `b6ba562ccadd96026d7c3425fe7f41fad526a207d936aaba7c847292c6a784b1`
- `SOLVENT.py`: `34763e03719b9dc2f4b1d138bc0eb4e72ab5863f7ae7c5d3dbe4e7d8356335b1`

Independent static review by Curie and root found no integration blocker. This
is distinct from an executed scientific or numerical test. Root's prior
MATRIX saved-artifact test4090 remains the existing algebra evidence.
