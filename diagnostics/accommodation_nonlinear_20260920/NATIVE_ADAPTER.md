# Native candidate adapter: ready for the final pilot collection

`scripts/accommodation_nonlinear_native.py` connects the existing nonlinear
candidate records, native ORCA runner and analytic-gradient parser. No new
scientific calculation or candidate preparation ran during implementation.

**Five real-artifact tests pass**, zero skips, in 1.513 seconds. They replay actual
native water-gradient drivers and unit conversion; current-coordinate projection
using archived GFN2 gradients; native PLM work/sign algebra; the eight original
DFT state records; and the exact EnGrad-only input change. Malformed-input checks
use explicitly corrupted copies of real records. The water-gradient fixtures used
TightSCF historically; their parser evidence does not change the new native recipe.

The actual origin inventory contains **8/8 completed native q0 endpoints**,
including both PLM La origins. All have TolE=1e-6 and TolG=5e-5. The 1H4I and
4MAE contexts contain 154 and 202 atoms and match their source coordinates
exactly. The 168-atom second PLM mapped origin has the previously observed
2.22e-16 Å arithmetic difference from its copied source; the actual executed
coordinates are retained. Inventory and exact source/receipt pins:

`workspaces/accommodation_nonlinear_20260920/native_adapter_preflight_v1/origins.json`

## Operations and safeguards

- `origins`: read the original expanded native references and explicit missing
  statuses. No molecular calculations.
- `prepare`: require the actual final immutable pilot collection and terminal
  execution receipt. Preserve all eight statuses and select only each optimizer's
  final candidate. Missing candidates remain unrun. Create at most two fixed
  four-task runner manifests, plus a combined collection manifest.
- `dry-run`: check candidate coordinates, physical mappings, inputs, original
  receipts and runner compatibility.
- `execute`: accept only a prepared shard and require at least one Ca/La pair
  with both cheap minima qualified. Use the existing runner; this command does
  not submit an allocation. Parent coordination still precedes submission.
- `collect`: retain all eight states, unavailable energies, native/component
  work, Ca-minus-La work differences, actual SCF-tolerance changes, and native
  versus cheap residual torques. No new accuracy threshold or native minimum
  qualification is invented.

Each input is exactly the original native header plus EnGrad:

```text
! r2SCAN-3c NoAutostart CPCM(Water) DefGrid3 EnGrad
```

Collection requires actual analytic SCF/XC/CPCM, dispersion and gCP output
markers, and an ECP gradient for La. It rejects numerical differentiation,
checks the complete gradient's energy/order/coordinates, and projects through
the cap-aware Jacobian at the actual candidate q. It does not reuse the q0
projection. All inactive physical modes are identified explicitly.

## Runnable checks now

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
-m unittest discover -s tests -p 'test_accommodation_nonlinear_native.py' -v

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
scripts/accommodation_nonlinear_native.py origins \
--source workspaces/accommodation_torsion_20260920/prepared_v3/design.json \
--controls-manifest workspaces/second_shell_20260919/prepared_v2/manifest.json \
--torsion-manifest workspaces/accommodation_torsion_20260920/prepared_v3/dft/manifest.json
```

`prepare --help` exposes all required paths: `--collection`, `--execution`,
`--agreement`, `--controls-manifest`, `--torsion-manifest`, and `--output`.
No source-code edits or shell defaults select inputs.

Full preparation and runner preflight remain **unrun until the real final pilot
collection exists**. No synthetic candidate or successful native output was used
to claim end-to-end readiness. The final collection will provide the concrete
prepare/execute/collect commands under the separately frozen
[native validation plan](NATIVE_VALIDATION_PLAN.md). No additional DFT has been
submitted, and production/default behavior is unchanged.
