# Matched water force sensitivity: all-electron versus frozen-f ECP

**The reference recipe materially changes these fixed-geometry water forces.**
The recovered Dy result and three existing matching results replay from their raw
outputs and gradients. All four have exactly the same coordinates and charge.
The comparisons below are descriptive; no new accuracy threshold is introduced.

Define outward displacement as rigid translation of all three water atoms along
the metal→oxygen unit vector. The reported slope is
`dE/dR = −Σwater F·u`, in kcal/mol/Å. Positive slope means an inward radial force
on the water.

| Reference recipe | La slope | Dy slope | Dy−La slope |
|---|---:|---:|---:|
| All-electron scalar ZORA, native PBE0-D4 | +1.548472 | +36.306189 | +34.757717 |
| Frozen-f ECP, native PBE0-D4 | −4.167242 | +32.664979 | +36.832222 |
| All-electron minus ECP | +5.715715 | +3.641210 | −2.074505 |

La's radial force changes sign between recipes. The individual atomic-force
difference has vector RMS4.055804 forLa and2.880571 forDy; the paired difference
has RMS1.539178. Partial metal cancellation leaves a finite change in the paired
response. Full Cartesian differences and rigid-body responses are retained in
the comparison receipt.

This does not isolate an f-electron effect. The recipes simultaneously change
core treatment, scalar-relativistic Hamiltonian and metal/ligand/auxiliary bases;
Dy is an explicit sextet in the all-electron calculation and an effective singlet
in the frozen-f calculation. Neither electronic-root choice nor spin–orbit
accuracy is established. The all-electron Dy force also lacks an independent
displacement test after its substantial final-COSX energy update. Do not replace
the consistently used provisional ECP development targets with these numbers
without a separately declared qualification.

Absolute energies are retained only for raw lineage. Raw different-element or
all-electron-minus-ECP total-energy subtraction is not a binding contrast. This
four-complex comparison has no matched all-electron fragment cycle, so it does
not supply a balanced binding energy or exchange free energy.

## Existing W8/W9 ECP energy-stage check

The requested one-time read of all four original ECP outputs found these final
COSX corrections, in kcal/mol:

| Configuration | La | Dy |
|---|---:|---:|
| W8 | +0.0284833 | +0.0209450 |
| W9 | −0.0228928 | −0.0146624 |

None approaches1kcal/mol. Their contribution to the balanced
`(Dy−La)W9−(Dy−La)W8` contrast is+0.0157687kcal/mol, far below the observed
−4.9512304kcal/mol model error. These saved corrections provide no evidence for
an analogous large energy-stage effect in that coordination screen. They do
not independently qualify the W8/W9 derivatives, and no diagnostic correction
was subtracted from any reference or fitted target. No new DFT run or blanket
halt follows from this metadata check.

## Evidence

Matching workspace `paired_comparison_v1/RESULT.json`:
SHA`c09e6906d7d6e11e2834179031b0c134601f2079417ad83abe6c58108a346ab8`.
The receipt contains all four source/input/output pins and separate SCF/energy-stage
metadata. `ECP_COORDINATION_STAGE_REVIEW.json`:
SHA`8c478786b53ddcf53deda1f74b9b46ce1433a29481c699e34ea50b5846b0e7e1`.
Both analyses used only saved results; zero new molecular calls.
