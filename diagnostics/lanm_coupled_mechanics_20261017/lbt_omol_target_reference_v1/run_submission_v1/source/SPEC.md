# One central reference using the published OMol25 recipe

21 October 2026 UTC. Root declares this finite experiment under the active
coupled-mechanics PLAN and Jacob's autonomous method-development authorization.
It follows source verification in TARGET_RECIPE_EVIDENCE.md. The separate PBE0
directional test retains its failed strict gate; this is a different method.

## Question and fixed input

How much of the model/reference force discrepancy remains when the comparator
uses the published electronic recipe underlying MACE-OMOL?

Use exactly the whole 256-atom 7CCN La step725 terminal already consumed in
reference15545: identical physical atom IDs, coordinates, protonation, charge−2,
physical multiplicity1, ECP46 and 1030 explicit electrons (515 alpha/515 beta).
No atom selection, alignment, optimization or new preparation of the molecule.
Reuse the saved MACE and PBE0 vacuum forces and the previously frozen full-vector
direction `u=(F_MACE-F_PBE0)/norm(F_MACE-F_PBE0)`.

This tests a published-recipe **port to ORCA6.1.1**, not exact original ORCA6.0.0
training-label reproduction. The backend, internal libraries, defaults and
electronic solution may differ. This consumed, conditional vacuum dianion is
not a blind benchmark or a validated physical ground state.

## Consequential settings

- One central `EnGrad`, `wB97M-V def2-TZVPD RIJCOSX def2/J DEFGRID3`.
- Complete external `def2-tzvpd.bas` from FAIRchem commit
  `d4e10d0f4e87593922fda1d1ffb7dd28e504cadd`, SHA256
  `7389c56f0cff9e76098b0757a34d1abdd79c48e9249597cea074a0106a3462a6`.
  Preserve its exact coefficients; use the recipe's built-in auxiliary/ECP
  selection and inspect actual native assignments.
- `UKS`, beta HOMO/LUMO rotation `{514,515,20,1,1}`; `NoUseSym`, `DIIS`,
  `NOSOSCF`, `NormalConv`; `%scf Convergence Tight maxiter 300` with
  `THRESH 1e-12 TCUT 1e-13` as in the published generator.
- Explicit `PAtom` and `NoAutostart` in a fresh directory preserve the documented
  unseeded default. No old PBE0 GBW or density is reused.
- Post-SCF VV10 as provided by `wB97M-V`; no D3/D4 or `SCNL`. Preserve DEFGRID3
  and the version's nonlocal-grid defaults; capture actual headers.
- Follow the public singlepoint wrapper's no-NBO/no-NPA convention, keeping
  population/spin/density evidence useful for this scout. No point charges,
  continuum solvent, MM, GFN or other correction is added to its vacuum energy.
- Pinned installed ORCA6.1.1 executable and existing generic native runner;
  a separately reviewed method-specific admission replaces old RHF/PBE0 gates.
- One CPU zero-call preparation, followed only after actual review by one
  24-MPI-rank native allocation, `banfield/standard/standard_expand`, one
  numerical thread per rank. Default scheduler memory/time, allocation-derived
  MaxCore, no GPU. The public wrapper's 12-process default is an execution
  setting; 24 ranks here provide default ~96GB memory headroom for the diffuse
  UKS scout. Actual costs and peak resources must be recorded.

## Execution and admission

Preparation validates source/hash lineage, exact physical identity, saved
electron counts, complete external basis and deck, and makes zero engine calls.
Root reviews that result before releasing exactly one native invocation.
Fresh output directories and a root execution decision prevent accidental
duplicate work. Failed or unavailable results remain immutable; there is no
automatic restart or settings recovery.

Native admission requires unchanged inputs/runtime pins, normal termination,
native SCF convergence, one finite energy and complete finite 256-atom gradient,
correct charge/counts/UKS/ECP46 and the actual declared method/basis/dispersion.
Inspect actual rotation/guess, orbital and auxiliary basis assignments, XC/COSX
and nonlocal settings. Preserve native residual criteria, density-integral,
frontier-gap and S² evidence whether favorable or unfavorable. Missing or
unrecognized evidence remains unavailable/pending saved-output review; it does
not authorize rerunning chemistry merely to repair a parser.

Native admission does not qualify force accuracy, electronic stability, spin
purity or physical accuracy. Compare full vacuum forces by RMS/max difference
and their projection on the old frozen direction. Report old PBE0 and new-recipe
discrepancies together, with explicit force conventions and units. Do not treat
absolute energies across these different methods as deformation work. No new
MACE, ALPB, GFN, MM, Dy or optimization calls occur in this experiment.

## Decision after the measured scout

A smaller model/recipe force discrepancy supports the old Hamiltonian mismatch
as part of the discrepancy. A persistent difference motivates model/target
investigation but still combines backend-version, electronic-solution and
derivative uncertainties. Neither result uniquely attributes functional,
basis or spin effects, validates shared ALPB, or establishes metal preference.

Only after actual method/state admission and measured resource review may root
separately release the two existing `x±0.01u` fixed singlepoints using this
recipe. If released, their derivative rule remains
`abs(F_FD-F_analytic) <= 0.005 + 0.0001*max(abs(F_FD),abs(F_analytic))`
kcal/mol/Å, with no automatic step/grid/tolerance sweep. Those two calculations
are not released by this specification or by central scheduler success alone.
There is no new optimizer, source search, production change or preference claim.
