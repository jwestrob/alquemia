# Native candidate adapter: eight candidates prepared, execution gate failed

`scripts/accommodation_nonlinear_native.py` connects the existing nonlinear
candidate records, native ORCA runner and analytic-gradient parser. All eight
actual final candidates are now prepared in two runner-compatible shards.
Both four-task dry-runs pass. **No native validation was submitted:** none of
the four cases has both cheap endpoints qualified, so the frozen execution gate
fails. Both shard execution calls were checked to refuse before the runner starts.

**Seven real-artifact tests pass**, zero skips, in 2.955 seconds. They replay actual
native water-gradient drivers and unit conversion; current-coordinate projection
using archived GFN2 gradients; native PLM work/sign algebra; the eight original
DFT state records; the exact EnGrad-only input change; and runner path containment
using an actual completed 1H4I Ca candidate; and final-collection identity/gating.
Malformed-input checks
use explicitly corrupted copies of real records. The water-gradient fixtures used
TightSCF historically; their parser evidence does not change the new native recipe.

The actual origin inventory contains **8/8 completed native q0 endpoints**,
including both PLM La origins. All have TolE=1e-6 and TolG=5e-5. The 1H4I and
4MAE contexts contain 154 and 202 atoms and match their source coordinates
exactly. The 168-atom second PLM mapped origin has the previously observed
2.22e-16 Å arithmetic difference from its copied source; the actual executed
coordinates are retained. Inventory and exact source/receipt pins:

`workspaces/accommodation_nonlinear_20260920/native_adapter_preflight_v1/origins.json`

## Shard layout repair, before scientific execution

The first adapter placed task files under the common output directory while its
runner manifests lived one level deeper in shard directories. The existing
runner correctly requires every task input, XYZ, output and gradient artifact to
remain inside its own manifest directory. Parent review caught this incompatibility
before any native validation submission.

Tasks now live under `shard_N/tasks/CASE__METAL/`, alongside their owning shard
manifest. Each shard keeps its own `execute.lock` and execution events. The
combined inventory references the children; no scientific input bytes changed.

Both shard-location preflights pass using unchanged bytes from the actual
completed 1H4I Ca nonlinear candidate, staged solely to check layout. The test
also confirms that a corrupted manifest pointing outside its shard is rejected.
No final pilot collection was fabricated, no endpoint output was created, and
these staging manifests are not execution candidates. Actual staging receipt:

`workspaces/accommodation_nonlinear_20260920/native_adapter_layout_v1/PREFLIGHT.json`

## Genuine final-collection preparation

Final source: `workspaces/accommodation_nonlinear_20260920/collection_final_v1.json`,
with `pilot_v2/execution_1204162.json`. Final reporting adds donor-contact and
curvature-component descriptions beside the immutable endpoint records. The
adapter permits exactly those two report-only fields and still requires every
original candidate field to match its terminal receipt. The first prepare attempt
rejected those additions before creating outputs; no calculation was involved.

Actual prepared inventory:
`workspaces/accommodation_nonlinear_20260920/native_prepared_v1/manifest.json`.
Its `PREFLIGHT.json` records both runner passes and explicit execution refusals.
`collection_unrun_v1.json` retains all eight candidates with native energies,
gradients and work **unavailable**, not zero. No endpoint output, gradient or
execution event was created. Actual analytic-driver checks on these candidates
remain unrun; compatibility was tested on archived native gradient outputs only.

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

The prepared inventory can be checked without executing chemistry:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
workspaces/accommodation_nonlinear_20260920/native_prepared_v1/implementation/accommodation_nonlinear_native.py \
dry-run --manifest workspaces/accommodation_nonlinear_20260920/native_prepared_v1/manifest.json
```

Do not submit these shards under the frozen [native validation plan](NATIVE_VALIDATION_PLAN.md):
its physical gate failed. Readiness of the runner does not change that scientific
decision. No synthetic candidate or successful native output was used. Production
and the separate proposal model remain unchanged.
