# Saved local force comparison

Declared before reading the new force values. Root owns the saved-only Slurm
launch. No molecular calls, optimization, geometry changes or new force gate.

Compare the four exact matched cores, 7CCO/7CCN × La/Dy, using all 48 atoms and
the fixed original Cartesian axes. Reuse admission of the original MACE12471,
GFN12825 and four PBE0 endpoint receipts through the versioned local-reference
energy comparator. Require the completed local MACE force result to retain
the same native model, physical state, source/XYZ/cap-map pins and existing
1e-7 eV energy parity gate. This gate checks replay only, not force accuracy.

Read each actual core map to partition atom indices into target metal 0,
other mapped real atoms, and synthetic cap atoms. Also retain the all-atom
summary. Do not infer caps from element names, omit inconvenient atoms, align
the two different sources or select axes from DFT outcomes. Compare MACE and
PBE0 only within the identical source/metal coordinate frame.

MACE supplies force in eV/Å. The admitted PBE0 parser supplies gradient in
kcal/mol/Å, converted to force by negation and division by
23.06054783061903 kcal/mol per eV. For explicit native-unit provenance, verify
the saved engrad array with the unchanged parser: Hartree/bohr gradients
multiply by 627.509474 / 0.529177210903 to give kcal/mol/Å. Then negate for
force. Keep the native gradient, both converted force vectors and the signed
MACE-minus-PBE0 difference for every core atom. Report the target's three
Cartesian components directly.

For each predetermined group report the number of atoms, each method's mean
force vector, the mean signed difference vector, and RMS atom-vector difference
`sqrt(sum_i ||F_MACE(i)−F_PBE0(i)||² / N)`, in eV/Å. These are descriptive
quantities, with no pass threshold, ranking, fitting or biological decision.
These are raw cluster-atom forces, including separately identified artificial
caps; they are not a cap-chain-rule assembly of whole-peptide physical loads.
Missing force data leaves that method/cell comparison unavailable. Compatible
GFN forces may be included only if the exact saved core receipts already
contain them; no force regeneration or whole-peptide substitution is allowed.

Keep all four source/metal slots and physical/map identities. Preserve native
PBE0 SCF residuals and warnings, including the distinction between admitted
energy convergence and unqualified force consistency. MACE's physical La1/Dy6
input and PBE0/GFN frozen-f effective singlets remain different electronic
treatments. Agreement or disagreement is a conditional local diagnostic;
it does not establish physical Dy-state accuracy, force qualification, affinity
or a transferable protein preference profile.
