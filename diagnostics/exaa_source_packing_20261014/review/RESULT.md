# Integration review and actual saved-artifact test

**No remaining static integration blocker identified.** Root's new fixed-pool
bridge preserves the qualified MACE/standalone-GFN2/native-ALPB methods and
frozen numerical gates. It runs MACE and seven independent one-thread solvent
workers in separate phases of the same seven-CPU allocation. The wrapper now
uses `--ntasks=7`, matching its runtime allocation check.

Actual Slurm job **4090** completed successfully on one CPU. All nine real
target matrices reproduce archived per-metal work, same-geometry paired work,
operational and mathematical selections within the declared 1e-8 kcal/mol
arithmetic tolerance. All 54 intentionally incomplete copies of those fixtures
produce unavailable pools, with no classification or successful fallback.
The analysis measured **0.330972 seconds**; scheduler elapsed rounds to zero
seconds and is not a claim of zero resource cost. **Zero molecular calls.**

The optional test of an actual marginal minimum retained as origin was
explicitly **unrun**: no such example exists in the declared saved panel. No
energies were fabricated to construct one. The strict `< -0.1 kcal/mol` rule
was inspected directly against the unchanged scorer; the actual regression
tests cover the selections that occur in these real matrices.

The review also checked the two preparation interfaces. Crystal v2 uses an
element-, graph- and donor-checked CCD/PQQ atom-name alias map, preserving actual
coordinates and existing state admission. Its proposed oxidation state remains
conditional. Packing preparation exposes complete candidate-specific maps and
paired states accepted by the endpoint adapter. These static interface checks
do not themselves establish new-source chemical or energetic validity.

The bridge explicitly gates xTB/tblite solvent-component and population parity
before unchanged native water-diagonal/vacuum-charge replay. The preserved
absolute-energy parity failure is not converted into a pass. New packing
contexts receive no inherited classification bands. Any missing required cell
makes the complete packing pool unavailable; raw successful origin artifacts
remain separately inspectable.

Exact receipt and source pins are in `RESULT.json`. Root owns actual new
endpoint execution and its independent completion records. This review does
not claim an improvement in discrimination or a measured new-site runtime.
