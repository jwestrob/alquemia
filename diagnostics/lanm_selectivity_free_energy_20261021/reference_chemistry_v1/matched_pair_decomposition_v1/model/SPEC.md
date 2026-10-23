# Literal extracted-pair model energies

Evaluate exactly 34 four-atom fragments from the shared reference PANEL: each
literal W8/W9 water with La/Dy, at unchanged coordinates, with the exact failed
17676 selected water keys and parameter file. Reuse four selected cluster
energies. No geometry rotation, parameter change, fitting, gradient, finite
difference, new reference, donor or held-out cell is included.

Copy the selected source key/prm bytes unchanged. Extract atom indices supplied
by the shared panel from its parent typed XYZ, renumber only indices and bonds,
and retain coordinate tokens. Require their numerical values to equal the
shared ordinary XYZ/geometry exactly. Each graph is isolated metal plus O–H,H.
The same-metal W8/W9 selected keys/prm must match byte-for-byte.

Use the pinned native Tinker analyze executable with `system.xyz E`. Adapt only
the analyze-energy/component part of the existing RUN_PANEL parser; preserve
version, finite-value and nonconvergence checks. Total/component sum closure
must be within 1e-6 kcal/mol, a printed-output consistency check. No force
program is run. Thirty-four programs maximum, four single-thread workers in
one four-CPU Slurm task, existing 300-second per-program helper timeout,
default allocation memory/time. A first native/admission failure stops new
starts; preserve all 34 slots, completed siblings and failure receipts.

Bare model ions have exactly zero energy under this particular vacuum key.
The pinned executed route is `analyze.f:539` → `analysis.f:38–69,216–227`
→ `ehal3/empole3/epolar3`. `analysis` initializes every component to zero;
bonded terms have no graph entries. The nonperiodic double loops start at
`ehal3.f:171,199`, `empole3.f:190,232`, and `epolar3.f:165,222`, so no
distinct pair exists for one ion. It has no field or induced dipole. There is no Ewald,
periodic, solvation, restraint or external-field term in these exact keys/prm.
This analytic convention is explicit evidence, not a missing-result default.

For either method C=Dy9−La9−Dy8+La8. For the model, P is the same signed sum
of its 34 dimer energies minus the explicitly zero Dy−La own-ion term; B=C−P.
Report kcal/mol and component sums. The root collector joins reference values;
this adapter does not change fitting targets or infer a unique physical cause.
If any slot fails, model P/B is unavailable while C and all admitted slots remain.

Preparation replays four saved cluster analyze outputs, validates all literal
maps and selected inputs, and binds sources/binaries/shared panel. It performs
zero native calls. Runtime requires a separate exact root release and passing
preparation; source preparation alone does not authorize execution. All outputs
are fresh directories. Historical fit, ledgers, source and CP results stay intact.
