# Remaining LBT hybrid derivative checks

This finite extension is declared before its energies or derivatives are read.
It changes no electronic model, state, partition, solvent, MM ledger or numerical
gate. It does not authorize optimization. Root owns execution.

## Required evidence and denominator

The bound hybrid evaluator is `EVALUATOR_V2.py`, SHA256
`686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661`.
It retains the original model and numerical derivative gates. Its separately
reviewed `RECOVERY_V2.md` changes only derived cap-Jacobian replay after the
zero-call 12975 failure and actual-host diagnosis 12984. Exact task coordinates,
states, charges and raw component-gradient replay remain unchanged. The exact
parent configuration and cache identity bind that helper and recovery evidence.
Its nine-cell 7CCO/baseline/La scout must have passed both derivative directions
at both declared steps, and its eight-origin matrix must be complete. Both
parents must have the same model identity, exact evaluator configuration and
raw component lineage. The only origin-configuration addition is its pinned
scout result. All nine scout cells and eight origin rows are replayed from saved
raw components before any new call.

The new denominator is exactly three cases, in this order:

1. 7CCO / baseline / Dy.
2. 7CCO / expanded / La.
3. 7CCO / expanded / Dy.

Each case reuses its admitted origin, then evaluates two directions, two step
sizes and two signs: **three reused centers and 24 displaced composites**.
This is a 7CCO qualification only; no derivative qualification of 7CCN is
implied. Missing or failed slots remain explicit. A backend/admission failure
stops new calls, preserving attempted failures and untouched slots. A derivative
gate failure retains the complete declared panel; it causes no retry, revised
direction, changed step or wider tolerance.

## Directions frozen by physical identity

Use the same source frame and original 256 real coordinates for both metals.
The saved baseline scout vectors are authoritative and must replay the frozen
evaluator's deterministic construction.

- **Boundary anchor:** baseline uses omitted `A/2//CA`, retained `A/2//C`.
  That omitted atom lies inside the expanded real electronic region, so the
  expanded case uses its first declared boundary: omitted `A/1//CA`, retained
  `A/1//C`. In either partition, normalize the vector from retained to omitted,
  choose the Cartesian axis with the smallest absolute dot product (first axis
  on exact ties), and normalize their cross product. Only that omitted real
  anchor moves. The cap and shifted field follow the unchanged mapping helper.
- **External charge:** the original scout's `A/1//N` remains an external,
  nonzero-charge backbone atom in both partitions. Reuse its exact physical
  index and targetward unit vector for all three cases. Verify its identity and
  external membership in both actual maps; do not substitute another atom if
  this check fails. Its partition-specific shifted charge remains unchanged.

These choices use source IDs and partition boundaries, never force or energy
outcomes. The adapter writes `FD_DIRECTIONS.json` before new evaluation. It
verifies common physical inventories/coordinates and exactly one moving real
atom per direction. No metal-dependent direction, cap coordinate or atom
selection is introduced.

## Unchanged finite-difference and state rules

Displace the selected real atom by `+h*v` and `-h*v`, with unit Cartesian `v`
and **h = 0.0002 and 0.0001 Å**. Rebuild all local real/cap coordinates, shifted
external charges and their coordinates, native MM coordinates and whole fixed
FF ALPB coordinates together through the unchanged evaluator. The central
analytic derivative is `dot(gradient, v)`, not force. Compute central energy
differences in native units before conversion, through the existing
`difference()` function.

For every component and the total at both steps, retain the original gate:

`abs(FD - analytic) <= 0.005 + 0.0001 * max(abs(FD), abs(analytic))`

in kcal/mol/Å. Report signed residuals and step-halving changes. No additional
step-agreement gate or force-convergence tolerance is introduced. Components
are local MACE vacuum, embedded-minus-isolated GFN2, the full MM ledger, and
whole fixed FF ALPB. Native point-charge gradient printing and all original
solver, charge, map and state gates remain unchanged.

MACE uses physical La multiplicity 1 and Dy multiplicity 6. GFN2 retains UHF 0,
effective multiplicity 1 and its published spin-independent f-in-core model.
All local and whole charges remain -2. Passing proves only these derivatives
of this model at the specified cells, not electronic-state accuracy, global
force qualification, a stationary minimum or experimental preference.

## Reuse, costs and execution

The existing evaluator is imported unchanged from the exact configuration
used by the admitted origins. Restoring the centers uses zero model, GFN, MM
or ALPB calls. Its full-composite cache and global-ALPB cache remain unchanged;
this adapter adds no component engine or new cache protocol. Exact admitted
parent ALPB cells populate that existing cache only after their binary, actual
coordinates, symbols, charge vector and native gradient receipt are checked.
Thus the baseline scout's compatible displaced global-solvent cells can also
be reused by the expanded partition; local partition components cannot.

Maximum new work is 24 MACE forwards, 48 standalone GFN2 calls, 24 MM
configurations / 96 MM group queries and 24 ALPB calls. Exact repeated ALPB
coordinates between the two regions may reduce the last count; report actual
calls and cache hits rather than assuming savings. The local electronic and MM
components remain partition-specific. Retain receipts for all attempted cells,
native child CPU time, wall time, model loading and parent provenance separately
from the reused parent costs. Use one warm MACE model, then the existing two
single-thread GFN workers, in the root-owned 12-CPU / one-GPU allocation.

CLI: `--stage preflight|run --config CONFIG.json --output FRESH_DIRECTORY`.
Preflight requires the same admitted parents and performs no new molecular
calls or GPU model initialization. Both stages require Slurm. Config fields
are `protocol`, `implementation`, `specification`, `plan`,
`evaluator_implementation`, `evaluator_config`, `scout_result`, `origins_result`.
The latter two must be actual terminal receipts; placeholders cannot admit a
run. Failed evidence is preserved and no alternate parent is selected.
