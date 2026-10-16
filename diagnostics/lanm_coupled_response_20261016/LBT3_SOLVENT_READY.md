# Whole-LBT3 solvent/gradient adapter readiness

16 October 2026. Source preparation and static review only; this note records
no new whole-peptide GFN2 execution or scientific result.

## Frozen source

- `LBT3_SOLVENT.py`: SHA256
  `e3a906103b2982aed5f99b4bc01f4fd02fdf0986f8611ba1863987a8df3a4373`.
- `RUN_SOLVENT.sbatch`: SHA256
  `4839581a86e0ef3676d3f97ff7bcfee89dc5453c1c5178ab68d663d58d9d4f97`.

CLI: `python LBT3_SOLVENT.py --config CONFIG.json --output FRESH_DIRECTORY`.
The wrapper takes absolute script, config and fresh output paths, in that order.
Root owns the pinned config, actual Slurm submission and completion wake.

Required config fields: `protocol_id`, `repository`, `implementation`, `plan`,
`preparation`, `mace_result`, `qualification`, `backend_helper`, `gradient_helper`,
`common_helper`, `python`, `pins`, `runtime`, `workers: 4`, and `tolerances`.
All artifact fields use path/SHA256 pins. Optional `mace_forces` pins actual12500.
Runtime is the exact `{executable, parameters, runtime_libraries, environment}`
from an accepted standalone task in qualification12488.

Fixed tolerances:

```json
{
  "charge_closure_e": 0.00001,
  "xtb_json_stdout_hartree": 0.0000001,
  "engrad_energy_hartree": 0.00000001,
  "engrad_coordinates_A": 0.000001,
  "xtb_json_population_e": 0.0000001
}
```

The final named tolerance checks two outputs of the same endpoint, JSON partial
charges and the full charge file. It changes neither backend parity acceptance
nor any scientific score.

## Scope and contracts

The exact four whole256-atom origins from preparation12469 yield eight calls:
7CCO/7CCN × La/Dy × vacuum/ALPB water. Charge−2 and GFN2 effective singlet/UHF0
are explicit; physical La singlet/Dy sextet metadata is retained separately.
The qualified standalone settings and runtime are reused without native ORCA
calls. Actual12488 was read selectively: all eight endpoints, four absolute
energy/population comparisons and two solvent-component comparisons passed.
That small-core result is implementation evidence, not whole-peptide force or
affinity validation.

The adapter calls the existing process helper and analytic `read_engrad` parser.
It checks actual method/state/settings, charge closure, finite full populations,
normal termination, gradient units, order, energy and frozen coordinates. Each
gradient is pinned before parsing and verified afterward. A numerical gradient
or missing/unsupported analytic artifact remains unavailable. The effective
subprocess environment reproduces the declared inherited thread/locale/PATH
controls from qualification and its standalone overrides.

Outputs preserve every endpoint and all four solvent corrections. When12500 is
provided, exact source/state/model/version/energy checks precede force reuse:
`F_composite = F_MACE - (g_water - g_vacuum) * hartree_to_eV / bohr_to_A`.
No extra classical energy is added. Source contrasts, solvent effects and the
common two-source matrix/minima retain both geometries and explicit availability.
No Mex composite, affinity or experimental population is inferred. The later
declared directional derivative check remains necessary before optimization.

## Verification

Python AST and shell syntax checks passed without importing or running the
adapter. Independent child static review found no blocking schema, sign, unit,
source-state, environment or pool-selection issue at the hashes above. Actual
whole-source, force12500 and qualification schemas were inspected read-only.
No scientific computation or submission was performed by these source/review
agents; execution remains root-owned. After independent review, root requested
one accounting-only field rename to `application_allocated_CPU_seconds` so the
application's elapsed-time estimate is distinct from exact scheduler allocation
cost. The final hash above includes only that additional rename; scheduler
receipts remain the source of actual allocation cost.
