# Executable adapter proposal

This is an implementation handoff, not an executable DFT launcher. Root requested
a bounded inventory while the derivative scout runs. The eight exact source
states are already available; no geometry generator or new energy engine is
needed. Do not pass the response task IDs to the frozen derivative `QUALIFY.py`:
its `IDS`, single-atom displacement assumptions and collector denominator are
specific to the eight derivative fixtures.

## Reuse points

| Function / artifact | Role in the narrow new adapter |
| --- | --- |
| `mex_dft/QUALIFY.py::inputs` with its original pinned config | Replay the native method, origin state, assets and allocation-aware runtime; no new origin call |
| `mex/MEX_CARTESIAN_POOL.py` and `mex/MEX_CARTESIAN.py` saved admission | Replay common5, actual broad search terminals, current full physical state and both-metal common7 evidence |
| Eight saved `trial`, `mapping`, `source_inputs`, `xyz`, `pointcharges` pins | Copy existing current geometry/environment exactly; no `update_geometry` call or source-coordinate reconstruction |
| `mex_dft/native/run_orca_task_manifest.py::run_manifest` | One endpoint and 12 MPI ranks per native allocation |
| Native runtime renderer and reference SCF parser already pinned by the DFT config | Exact state/recipe, allocation-derived MaxCore, native convergence and normal-termination evidence |
| Existing `metal_environment_reference.read_pcgrad`, `affordable_response.read_engrad` and `metal_environment_force_assembly.scatter` | Read analytic QM/PC gradients and scatter using the archived current cap Jacobians and physical IDs |
| Actual current `MM_result` and `ALPB` receipts in each inventory cell | Reuse energy and gradient after exact physical/ledger/charge replay; no recalculation |

The frozen derivative `admit_sp` admits energy/state/SCF evidence but does not
admit the analytic gradient arrays. A response energy-gradient admission must
also use the existing gradient-native reference path or a narrow versioned
wrapper. Preserve the energy/state/SCF gates and add explicit QM/point-charge
gradient, unit and atom-order checks. Record the actual `EnGrad` task type.

## Proposed four-stage interface

```text
RESPONSE_REFERENCE.py --stage prepare --config CONFIG.json --output FRESH_PREPARATION
RESPONSE_REFERENCE.py --stage preflight --config BOUND_CONFIG.json --output FRESH_PREFLIGHT
RESPONSE_REFERENCE.py --stage endpoint --task La_terminal__La --config BOUND_CONFIG.json --output FRESH_ENDPOINT
RESPONSE_REFERENCE.py --stage collect --config ACTUAL_RESULT_PINS.json --output FRESH_COLLECTION
```

All four stages execute within root-owned Slurm allocations. These are proposed
commands; `RESPONSE_REFERENCE.py` has not been implemented or submitted.

### prepare / preflight — zero molecular calls

1. Require actual common7 12981 complete and both broad searches admitted. Match
   the frozen alias order and all eight POINT/physical pins in `INVENTORY.json`.
2. Replay native approximate-model evidence only as source/ledger provenance.
   Recover the complete original-source state, water/spectator identities, full
   covalent graph and bound-state admissions from the frozen Cartesian layer.
3. Reuse the eight current archived core/field/cap mappings. Verify exact physical
   arrays, core XYZ, field charge/coordinate order, current boundary anchors,
   target element/state, MM ledger and full solvent charge/coordinate arrays.
   Geometry equality must include all 2,049 atoms, not just Q209.
4. Copy the original same-metal PBE0 `EnGrad` input bytes to a new task directory,
   and copy each selected core/field byte stream to its expected filenames.
   Bind the same native executable, basis/ECP/auxiliary assets and spin semantics.
5. Produce a full eight-task preparation manifest and eight explicit cell slots.
   Before each endpoint execution, require actual derivative scout collection
   plus root's measured tractability record. An execution manifest contains only
   its one named task. Missing evidence fails closed with no molecular call.
6. Pin the two full physical displacement vectors and the intended force-error
   reporting/precision policy before the first response result. This stage must
   not claim those directions have already passed independent DFT checks.

### endpoint — one energy plus analytic gradient

Reuse MM and solvent from the admitted same-metal paired POINT. The new ledger is
only embedded PBE0 plus those two saved terms. Preserve the current physical maps
and attach them to native task/receipt metadata; do not copy origin environment
provenance onto displaced coordinates. Preserve native failure artifacts and
every original input. No automatic restart or retry route is part of this API.

Admit native QM and point-charge derivatives with the existing parsers, cap
Jacobian scatter, units and sign. Record physical total gradient, component
energies/gradients, full native state, SCF checks and runtime. A native energy
convergence ending does not erase failed density criteria or prove stability.

### collect — saved data only

Require exact task/result links and the complete eight-cell denominator. Produce
within-metal work comparisons, fixed four-geometry paired pools, components,
force projections and uncertainty flags using the frozen declaration. Partial
results remain explicitly partial. No new MM/ALPB or reference computation is
launched by the collector. Actual scheduler cost is joined separately without
counting reused components as new calls.

## Remaining implementation gates

- Native derivative scout is still running; no response numerical gate is met.
- The existing gradient receipt/admission helper must be pinned in the new
  adapter and independently reviewed against the actual `EnGrad` schema.
- The source-only cache inventory establishes no reusable response PBE0 cell;
  a new exact archived result must pass the full same-state native admission
  before reducing the eight-call denominator.
- Root must bind a finite scout/runtime decision. No implicit eight-job launch
  or attempt to infer an ETA from the running pair is included.
