# Scaffold feasibility — fixed three-source inventory

Parent assignment, 22 September: determine whether an available standard-protein
parent can supply only the mechanical environment omitted from the compact MACE
context. Fixed sources: 1H4I, 4MAE, PQQSEQ_83440678cbbd658047c9 from the original
30-context archive. No new molecular energy/force calls, parameter fitting, H
normalization, optimization, dynamics, Hessians or scoring.

Inventory exact source-to-context/cap mappings and actual ff19SB bonded/nonbonded
term supports. Reuse the already recorded OXT repairs for the two crystal parents
only; keep every source coordinate/proton. Old whole-protein H normalization is
explicitly excluded. Parameter assignment and serialization create no OpenMM
Context and evaluate no energies/forces. Unknown chemistry fails explicitly.

Determine whether boundary subtraction is actually available. If not, report the
missing cap/reference and environmental interactions without manufacturing a score.
A fixed parent-force diagnostic may be proposed to root, but is not executed here.
This branch is separate from missing-metal-translation tests. No claim that a
0.8 Angstrom proposal boundary diagnoses a unique physical cause.
