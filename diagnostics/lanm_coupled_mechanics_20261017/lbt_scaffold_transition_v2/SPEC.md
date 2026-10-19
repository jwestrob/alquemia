# Scaffold transition audit V2: exact force-representation replay

## Unchanged scientific question and scope

This is a technical correction to the pinned V1 SPEC and ROOT_SCOPE_v2. The same four exact 7CCO256 before/after × La/Dy cells,187 terminal-touching bonded interactions per state,36-context/816-query ceiling, fixed parameters/coordinates, full denominator, unavailable before-Dy whole force and after-Dy hybrid, report and interpretation remain unchanged. All V1 MM energy/force closure thresholds remain exactly1e−8 kcal/mol and1e−5 kcal/mol/Angstrom. Root owns the one-CPU Slurm submission; no new molecular method or extra query is introduced.

## Preserved failure and concrete source diagnosis

V1 source79d311d8, config84060852 and actual job14785 result626252ba remain untouched. That result failed `Whole comparator/ALPB force replay differs` before any MM context/query. The combined gate incorrectly demanded equality between two historically distinct floating-point representations of the same native ALPB gradient. The original source review missed this operation-order difference.

The original producers are explicit:

- `REFERENCE.py::admit_alpb`, line243: stored reference gradient `(g * HA_TO_KCAL / BOHR_TO_A)`.
- `MODEL_COMPLETION.py::package`, lines259–261: stored whole-model solvent force `-g * (HA / BOHR)`, whole MACE force `raw_force_eV_A * EV`, then their sum.
- The actual full collector `CONTINUATION.py::force_report`, lines278–286 (SHA06e40d37), independently repeats that whole-model arithmetic and copies its admitted composite into the reference force table.

Metadata-only inspection confirms all three available stored whole composite arrays exactly equal their corresponding full-reference force-table arrays. The actual before-La first ALPB component is printed as +0.5506655383113366 in the whole force array, while the stored reference gradient prints −0.5506655383113365. No residual was calculated on login. The source formulas explain why the invalid bitwise cross-representation equality cannot serve as an admission test; this is not evidence of changed physics or inaccurate forces.

## Narrow correction

V2 pins the actual original collector, whole-model implementation and original constant source. Values remain EV23.06054783061903, HA627.509474 and BOHR0.529177210903. For each available whole force array:

1. Verify the same physical IDs, force sign/units and native MACE receipt.
2. Reconstruct `mf = raw_MACE_force * EV` and `sf = -raw_ALPB_gradient * (HA/BOHR)` using the original collector's operation order.
3. Require exact array equality for stored MACE force, stored whole ALPB force, stored composite versus `mf + sf`, and composite versus the original full-reference force table.
4. Separately require exact equality of stored reference solvent gradient versus `raw_ALPB_gradient * HA / BOHR` using its original producer order.
5. Preserve the exact shared native ALPB receipt/vector and the existing reference/hybrid component checks. Neither admitted reference force table nor whole-model force array is replaced or recomputed as a different representation.

No tolerance is added or relaxed. Original failed gates/results remain failed. All source-only checks are AST/hash/metadata checks; actual replay and MM analysis remain Slurm-only. A fresh V2 output directory is mandatory, with the same per-cell partial-failure policy and no retry. The existing V1 launcher is reused unchanged. CONFIG and READY bind every prior/source input and the exact command.
