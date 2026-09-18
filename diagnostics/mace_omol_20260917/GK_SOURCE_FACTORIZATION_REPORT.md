# GGR source-self shift is dominated by projected-charge changes

The exact symmetric decomposition of the verified2FW0-minus2FVY GK source-self
contrast gives:

| Factor | Contribution to Ca-minus-La difference, kcal/mol |
|---|---:|
| Projected charge distribution | +18.452792731 |
| Pairwise source distances | -0.372799829 |
| Effective Born radii | +0.049028048 |
| **Sum** | **+18.129020950** |

All source identities and formal charges match. Each endpoint uses all eight
combinations of the actual charge vectors, distance matrices and effective
radii, with marginal changes averaged over all six factor orderings. Both
unmixed energies reproduce the independently verified native source-only
energies; endpoint and paired closure pass1e-7kcal. Mixed combinations are
explicit algebraic counterfactuals, not computed quantum states or measurements.

This identifies a sensitivity to the projected charge distribution in this
particular solvent kernel. It does not distinguish actual electronic response
from CHELPG fitting instability, or prove that the point-source GK functional
responds correctly to either. It also shows that changes in effective radii
between these two structures contribute little to this particular discrepancy.
The common30A metal-radius saturation is therefore not a sufficient causal
explanation, although the validity of the common approximation remains open.

No new quantum, MACE, charge fit, solver or biological score. Analysis took
1.833442wall seconds; full CPU/RSS, raw combinations, per-atom charge/radius
changes and source pins are in the actual result. Baseline and failed hybrid
remain unchanged. [Compact result](GK_SOURCE_FACTORIZATION_RESULT.json).

Next: investigate documented native CHELPG sampling controls on the same
saved densities. The standalone installed utility accepts onlyGBW/density,
so numerical refinement requires an explicitly verified route; do not silently
swap charge representations or modify binary wavefunctions.
