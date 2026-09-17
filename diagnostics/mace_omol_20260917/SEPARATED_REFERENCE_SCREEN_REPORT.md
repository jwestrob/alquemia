# An isolated reference alone cannot repair the ordering

All10 existing intact bound endpoints were verified against their receipts.
For a charge-consistent separated reference, the common protein term cancels
between Ca and La. A shared ionic/aquo offset also cancels between proteins.
The resulting relative ordering must therefore match the raw bound contrasts:

| Comparison | Raw bound difference, kcal/mol | Required direction |
|---|---:|---|
| XoxF−MxaF | −623.299274 | Fails |
| alpha1F6S−GGR | +483.352723 | Passes |
| alpha6IP9−GGR | +467.590724 | Passes |

These values have no quantitative affinity interpretation. They show that
computing a new common reference cannot fix the PQQ reversal. The matched
coordination descriptor's useful ordering depends on the subtraction beyond
a mere zero shift; its separate charge-conditioning weakness remains.

No new molecular energy was evaluated. No isolated-ion energy, absolute
separated score, decision band or missing correction was supplied. This is
consumed development algebra, not independent validation. Full unrounded
results, exact script and timing are under
`workspaces/mace_omol_20260917/separated_reference_screen_v1/`.
