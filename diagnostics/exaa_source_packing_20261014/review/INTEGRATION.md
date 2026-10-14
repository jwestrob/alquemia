# ExaA endpoint integration review

14 October 2026. Read-only code/metadata inspection; no scientific execution.
Root owns implementation and Slurm submission. This is an interface review,
not a claim that the new packing or crystal sources have passed preparation.

## Recommended existing route

Reuse warm CPU OMOL, standalone xTB6.7.1/tblite0.7.0 endpoints, and the qualified
native xTB ALPB replay. Root elected to retain full per-cell water-diagonal and
cross-backend qualification: twelve packing XYZ/state cells require 24 inputs
per SCC backend (48 SCC calls), plus 24 native component evaluations and twelve
new MACE evaluations before exact deduplication. Two archived MACE fixture
evaluations per worker are additional actual work. No optimization is needed
for the packing branch. These counts are plans, not completed calculations.

This is **standalone xTB**, not ORCA-native GFN2. Its pinned published-parameter
file is already the corrected lineage used by the frozen candidate. Do not
silently switch its charges to ORCA or interpret a loadable model as parity.
If ORCA is ever invoked, route through
`scripts/guarded_native_gfn2_dispatch.py` and the full-file/shell-aware
`scripts/native_gfn2_parameters.py` contract; there is no native-default fallback.
No ORCA calls are planned here.

## Concrete reusable interfaces

1. Clone the scientific/runtime/fixture fields of
   `workspaces/ca_coverage_recovery_20261012/two_ca_scoring_q9_v1/lane0/worker_manifest.json`
   into a new finite task manifest. Its pinned `CPU_BATCH.py` provides
   `WarmCPU(manifest,m,attempt,lane)`,
   `cell(directory,request,transport,model)` and
   `native_cell(task,result,result_pin,xyz_pin,origin_label)`.
   A request is `{manifest,task_id,xyz,charge,multiplicity}`. Root tasks can be
   arbitrary pinned XYZs; both qualification tasks must be included. The worker
   checks graph charge, multiplicity, atom order, `omol` head and unchanged
   parameters, and resets the calculator for every request. Do not call the
   old `validate()` directly: it hardcodes the historical 106-source population.
   Existing `two_ca_scoring/NATIVE.py:execute` shows the finite per-source and
   two-fixture qualification pattern. Its gates remain 0.01 kcal/mol energy and
   0.001 eV/angstrom maximum force error.
2. `scripts/context_factorial_gfn2.py:execute_point(task,manifest,base,runtimes)`
   runs unchanged vacuum/water xTB and independent tblite. Tasks include ID,
   pinned XYZ, metal, charge, multiplicity, valence-electron/atom counts, `uhf=0`,
   medium, reuse/reference fields. Runtime and settings must come from the
   qualified corrected manifests, not fresh defaults.
3. Gate each vacuum/water pair with
   `expanded_transfer_inventory/REEVALUATE.py:component(base,vacuum_task,water_task)`.
   This verifies raw outputs and requires both backend states/convergence,
   solvent-component agreement within 1e-6 Eh and populations within 1e-5 e.
   **Do not gate on `execute_point`'s legacy absolute `accepted` flag:** the
   preserved absolute-v1 offset is a separate failed test. Do not erase it.
4. Reuse frozen candidate `implementation/NATIVE.py:run_task` after that gate.
   It replays water charges against printed native components, then evaluates
   vacuum charges in the identical cavity. `run_task` alone does not enforce
   xTB/tblite component or population agreement. Its `G_qvac_hartree` is the
   candidate solvent term; do not add the SCC water-minus-vacuum energy as well.

Native `task_preflight` requires: `xyz`, metal/state/electron/atom fields,
`state_signature` with endpoint charges, `source_preparation`, `maps`,
`scientific_input_key.runtime`, `uses[{MACE,MACE_eV}]`, and
`endpoints.{vacuum,water}` containing `point`, `point_manifest`, `task`,
`execution`, `stdout`, `charges`, `executed_xyz`. Exact eight-decimal charge
files, task state, commands and water setup are checked. Keep the existing
precision-derived bounds, 1e-9 Eh native allowance and 1e-6 Eh informative
ceiling. No invented shell populations or higher charge precision.

## Runtime identities

- OMOL model SHA256:
  `9b64b4fd5153ca578c694abc57806d8111050de6ff652e695c9b525bc4d36469`;
  float64, vacuum, nonperiodic, multiplicity convention, no extra dispersion.
  Model Python is `workspaces/mace_hybrid_20260916/software_v1/venv/bin/python`;
  optimizer Python is `/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python`.
  Reuse the worker manifest's exact calculator/source-module pins.
- `workspaces/orca_parameter_audit_20261012/context_gfn2_v1/RUNTIME.json`
  pins `/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb` (SHA256
  `481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f4`),
  its complete `param_gfn2-xtb.txt` (SHA256
  `f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545`),
  isolated tblite0.7.0 environment and all qualified dynamic dependencies.
- `workspaces/pede_solvent_response_20261013/frozen_distribution_candidate_v1/canonical_v1/BUILD.json`
  pins the already-built `replay` binary (SHA256
  `3d6a9d23aa1d1152c368b21f7f6d55a135279596b8d01a3a1b3ebc42617e75b6`)
  and linked libraries. CLI is `replay core.xyz charges ATOM_COUNT`.
  Preserve native ALPB/P16/water/gsolv, 298.15 K, 230-point grid, no ion screening.
  The source manifest pins all settings. The old `gsolv` token warning coexists
  with the explicitly printed correct effective reference; do not silently
  change command syntax in this comparison.

## Caching, failures and interpretation

`cell` reuses only exact request/model identities and does not repeat a started
call missing a terminal receipt. `execute_point` similarly refuses a preexisting
partial point directory; keep those visible. The native `run_task` creates new
output directories and is not an automatic resume facility. Cache identity must
bind source/state/maps, complete geometry, parameters/runtime and protocol, not
only a human-readable source name. Geometry deduplication cannot merge different
charges or proton inventories.

Run the seven-thread warm worker and parent scientific work sequentially in its
allocation. Independent SCC/native cells use one thread each, bounded by actual
allocated CPUs (current QOS at most12/job); do not overlap twelve SCC workers
with an active seven-thread model in a twelve-CPU allocation. All scientific
imports/tests/execution remain on Slurm.

Selection retains the common pool order origin, adaptive_Ca, adaptive_La; exact
ties preserve that order and operational selection requires work strictly below
-0.1 kcal/mol. Keep mathematical minima and independently available origin
results. For larger contexts compare within-context deformation works and
paired differences, not absolute totals or inherited calibrated class bands.
Crystal accommodation is conditional on actual admissible preparation and uses
the existing four-mode proposal/complete cross-scoring route; do not substitute
the fixed-packing endpoint scorer for that search.
