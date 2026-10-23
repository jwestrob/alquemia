# Water parameter correction and periodic mapping review

This saved-source correction supersedes the earlier description of the small-cluster water model as a clean renumbering into the protein namespace. It changes no prepared coordinate, parameter file, frozen report or execution release.

## Verified water difference

The small-cluster `amoeba09.prm`, line 967, has
`ureybrad 35 34 35 -7.60 1.5337`; the retained protein `amoebabio09.prm`, line 1425, has
`ureybrad 74 73 74 -7.60 1.5537`.
These are corresponding H–O–H Urey–Bradley records: their equilibrium H–H distances differ by **0.0200 Å**. The corresponding water van der Waals, permanent multipole, polarization, bond and angle records match after the explicit type/class renumbering in the periodic mapping inventory.

Both prepared exchange legs use the same unchanged **protein BIO-2009 water**. Their within-cycle water convention remains matched. Transfer from the small-cluster model is nevertheless not an identical full Hamiltonian, particularly for flexible-water sampling. Rigid translation leaves a water molecule's internal distances unchanged, so this internal term does not itself supply a rigid-translation force; that limited observation does not establish equilibrium or dynamical equivalence. No bonded term was changed to hide the difference.

## Actual class and type coverage

Independent inspection of the pinned protein parameters and both La box files agrees with the prepared periodic mapping inventory:

| Source SI family | Protein destination | Bound recipients |
|---|---|---:|
| Amide O class 55 | Class 5: types 5, 11, 153 | 19 |
| Amide C class 53 | Class 3: types 3, 9, 152 | 19 |
| Carboxylate O class 403 | Class 31: types 143, 161 | 10 |

The hypothetical VDWPR mapping is by **class** and therefore covers every represented member, including backbone, side-chain and terminal amide carbonyls. The corresponding POLPAIR mapping is by **type**, covers all listed oxygen types, and does not apply the acetate rule to protein water H type 403. No carboxylate-carbon pair override has been invented. Water O/H are protein types 402/403, classes 73/74. The boxes contain 3,894 bound and 3,980 aqueous waters. All proposed donor pair records remain unselected; this is an inventory, not physical transfer qualification.

The reviewed inventory preserves the distinct protein polarizabilities and permanent multipoles. Their differences from the small donor models still require explicit transfer qualification even when the chemical-family mapping is correct. Metal types/classes 900 and 901 are unused in the pinned protein parameter namespace.

Exact evidence pins, inspected source records and independent box counts are in `WATER_PARAMETER_REVIEW_v1.json` in the matching paired-leg workspace. This review made zero molecular calls and released no simulation.
