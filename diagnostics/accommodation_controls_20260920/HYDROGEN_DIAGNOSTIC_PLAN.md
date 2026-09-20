# MMOL1770 hydrogen minimizer diagnostic v1

Frozen2026-09-20 before replays. Parent approved this contained three-preparation
diagnostic under Jacob's overnight authorization. No scoring or default change.

## Evidence and question

Of17 failures,16 have ions about25 Å from PQQ; only MMOL1770 La sample4 has
collapsed H. Four preselected prepared controls have no sub-0.45 Å pairs.
The failed sample has1,759 such H/H pairs and4,685/4,687 H–parent distances
within0.002 Å of the installed OpenMM initial1.000 Å placement radius.
All hydrogens have exactly one heavy-atom topology bond.

Question: does the installed CPU50 hydrogen minimizer silently fail to move
its initialization, and can a numerical preconditioner resolve this while
preserving the exact same potential, atom inventory and fixed heavy atoms?

## Finite tasks and immutable sources

Existing `folds_v1/source_preparation/cases/`:

1. `mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-4` (failed):
   exact pinned protonator replay with a passive minimizer reporter; capture
   original System XML, initial/final coordinates, forces, energies, iteration
   count and topology. Compare replay to the existing output.
2. `mmol_1770-pqq-la_model__conditioned_La__seed-1_sample-3` (canonical success):
   the identical passive replay and comparison; an implementation control.
3. On task1's captured initial System/coordinates only, run the following
   deterministic numerical preparation, retaining all names/bonds/H inventory.

No other sources, RNG seeds, pH, force field or chemical states are tested.
Use pinned PDBFixer1.12/OpenMM8.5.1, pH7, RNG20260914, CPU1, forcefield=None.
Each independent source gets1CPU/8GiB on standard,memory. No GPU/DFT/MACE.

## Numerical repair declared before results

The System is exactly the one installed OpenMM generated: generic H-parent
harmonic bonds, oxygen-angle terms, and its repulsive nonbonded potential.
All mass-zero existing atoms remain exactly immobile. Only added H may move.
This is an inexpensive preparation potential, not a protein/metal energy model.

Starting from the captured native initialization, take deterministic gradient
descent steps, normalizing each movable atom's force to a maximum displacement
of **0.002 nm=0.02 Å** per step. Accept a step only if the original total
potential satisfies Armijo decrease `Enew <= Eold - 1e-4 sum(F dot step)`;
halve displacement up to20 times if needed. Stop this preconditioning once
no movable-H pair lies below the existing0.45 Å rejection threshold, or after
200 accepted steps. A nonfinite energy/force or failed descent is a failure,
not a reason to modify the potential or seeds.

Then invoke the existing native minimizer at its original tolerance1kJ/mol/nm
and50-iteration maximum, with an iteration reporter. No atom can move more
than **2.5 Å** from its original H position in the accepted final output;
this bound is a geometry guard, not a fitted restraint or extra energy term.

## Acceptance and interpretation

Require exactly unchanged heavy coordinates; identical atom names, bonds and
protonation; all H with one recorded heavy parent; no sub-0.45 Å pairs; all
H-parent distances within0.8–1.5 Å; finite energy/forces. Retain actual final
movable-force RMS and iteration evidence: passing geometry is not evidence
of force convergence. Passing repair is a numerical preparation result only.

Record all failures and receipts. Source code is inspected, not patched.
No changes to the all250 campaign, replacements, classifications or historical
artifacts. A later repaired protocol would require its own explicit version.
