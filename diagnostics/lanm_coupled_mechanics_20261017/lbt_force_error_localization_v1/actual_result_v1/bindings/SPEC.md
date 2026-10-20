# Localize shared and metal-dependent force discrepancies

## Question and inputs

Job15720 completed the eight existing whole-peptide force comparisons. Individual
MACE-minus-PBE0 force-error RMS values are about3.5–3.6kcal/mol/Å, whereas matched
Dy-minus-La error RMS values are0.476–0.727. Which represented atom categories and
residues carry the shared discrepancy and the remaining metal-dependent error?

Use the exact complete15720 RESULT
`f42cffcf5717a988ed569db925b840f02dfb8bd3b4c04750b77e33bc231c4790`,
its eight pinned force reports, and original inventory
`e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957`.
Retain all four source/geometry pairs and all256 atoms in each. Original7CCO
before/after and7CCN before/step80 coordinates remain unchanged; step725 and
pending15545 are outside this diagnostic.

## Fixed arithmetic and reporting

For each pair, use the saved **vacuum** vector errors
`e_La = F_MACE,La − F_PBE0,La` and `e_Dy = F_MACE,Dy − F_PBE0,Dy`.
Define their common part as `(e_La + e_Dy)/2` and their metal difference as
`e_Dy − e_La`. The shared solvent cancels from model-versus-reference errors;
this diagnostic does not validate solvent physics.

Retain all full vectors and per-atom norms. Report the full-system and group
atom-vector RMS, maximum, sum of squared norms and each group's fraction of
the full-system squared norm. Use these exhaustive, mutually exclusive groups,
in this precedence order:

1. The selected metal, ID`metal`.
2. All hydrogen atoms, from existing element metadata.
3. Backbone heavy atoms named`N`, `CA`, `C`, `O` or`OXT` in the final ID field.
4. Remaining protein heavy atoms.

All hydrogens remain together; no covalent donor assignment is inferred. Also
report every residue using the existing ID prefix before its final atom field,
with the metal separate. Show the ten largest common-part and metal-difference
atom norms for each pair, breaking ties by physical ID. The full atom table is
authoritative; the ten-atom display is a fixed summary, not a selected scoring
region. No distance cutoff, new geometric region or favorable pair is selected.

Preserve the eight original endpoint RMS values and four original matched-force
metrics. Reproduce the saved difference vectors and their norms. Check the
arithmetic identity
`||e_La||² + ||e_Dy||² = 2||common||² + 0.5||difference||²`.
Use absolute1e−10 and relative1e−12 tolerances solely for saved floating-point
arithmetic comparisons. They do not replace any original admission rule or
define a new physical-accuracy threshold.

## Execution and interpretation

One saved-only Slurm job on one CPU, standard_expand with scheduler defaults.
Verify file hashes, all eight available rows/four pairs, exact256 physical IDs
and matched metadata before calculation. Require finite full vectors; preserve
explicit unavailable slots and failure evidence if a check fails. No scientific
packages, model loads or electronic engines are needed. Zero new MACE, DFT,
ALPB, GFN, MM, geometry, optimization or parameter-fitting calls.

This localizes a force-error distribution. Squared-force fractions are not
energy contributions, and a large atom error does not identify a unique faulty
interaction. PhysicalDy6 MACE versus frozen-f/effective-singlet PBE0, native
SCF/gradient/state limitations, conditional termini and the vacuum-dianion
assumption remain. No threshold fitting, production change, closed-branch
restart, natural transfer or preference claim follows.
