# Both-chain preparation ready; scientific execution unrun

Manifest: `workspaces/exaa_source_packing_20261014/crystal_v2/manifest.json`

SHA256: `2d134b424085751ac92938e881ff1436f4b39ed06622eedc2ba0b4c06cda7ae7`.

Root submission from repository root:

```bash
sbatch --nodelist=node-64-768g-5 \
  --export=ALL,MANIFEST="$PWD/workspaces/exaa_source_packing_20261014/crystal_v2/manifest.json" \
  --output="$PWD/workspaces/exaa_source_packing_20261014/crystal_v2/slurm-%j.out" \
  --error="$PWD/workspaces/exaa_source_packing_20261014/crystal_v2/slurm-%j.err" \
  diagnostics/exaa_source_packing_20261014/crystal/RUN.sbatch
```

One CPU, one node, standard partition, no memory/time override. Root must arm the
normal completion wake for the actual job. Terminal artifact:
`crystal_v2/execution/FINAL_PREPARATION.json`; two rows, including failures.
No molecular scoring/optimization is included. The root can inspect actual
ADMISSION/source maps and conditional-state caveats before any later scoring.

Metadata-only `freeze` ran on login; no scientific modules, tests or coordinate
analysis ran there. The actual Slurm job performs real source/schema, coordinate,
state, fragment, map, spectator and mature-construct checks on both chains. The
preliminary v1 was never submitted; its source/plan are preserved in `unrun_v1`.
V2 adds the already-established exact PQQ name-map aliases and pre-preparation
import-pin validation. It changes no chemical or geometric admission threshold.

Detailed redox assignment in the primary article remains unavailable because its
PDF returned HTTP403. Current input explicitly tests the existing oxidized-PQQ3−
hypothesis; no experimentally verified redox equivalence is claimed.
