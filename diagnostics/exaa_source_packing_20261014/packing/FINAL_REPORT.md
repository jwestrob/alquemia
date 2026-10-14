# Added packing changes individual works, with little metal-specific effect

The complete fixed-geometry test does **not support omitted local packing as
the main explanation** for ExaA two-Ca sample3's unusually large La accommodation
response. Both sources and all12 new endpoint cells are available. Adding the
declared surrounding fragments makes movement less favorable for both metals,
but most of that change cancels in their difference. Neither source changes its
selected geometry under either mathematical or operational selection.

## Actual comparison

Every number below is an energy change from the **same context and metal's own
origin**, in kcal/mol. Negative work means a lower energy than that origin.
The differential is selected Ca work minus selected La work; it is the
accommodation contribution, not the full discrimination score or an affinity.

| Source | Selected work | Original192 atoms | Expanded249 atoms | Expanded − original |
|---|---|---:|---:|---:|
| sample0 control | Ca | −7.624675 | −6.402300 | +1.222375 |
| sample0 control | La | −12.934776 | −11.465679 | +1.469097 |
| sample0 control | Ca − La | +5.310100 | +5.063378 | −0.246722 |
| sample3 target | Ca | −9.698819 | −7.247497 | +2.451321 |
| sample3 target | La | −26.390560 | −23.564450 | +2.826110 |
| sample3 target | Ca − La | +16.691741 | +16.316953 | −0.374789 |

For both sources, Ca still selects `adaptive_Ca` and La selects `adaptive_La`.
The same choices follow mathematical minima and the unchanged0.1kcal/mol
operational origin-retention rule. The other candidate geometries were also
evaluated; their full works and differences remain in the comparison result.

The target's selected MACE work becomes less favorable by+3.293185 for Ca and
+3.371314 for La. The frozen-charge solvent work partly offsets those changes
by−0.841864 and−0.545204, respectively. Consequently, the target's differential
change comprises only−0.078128 from MACE and−0.296660 from solvent. In the control,
the corresponding differential changes are−0.285135 and+0.038413. The large
target La response therefore survives this context change; it is not simply
removed by charging the saved motion for the newly represented surroundings.

## What was qualified

Preparation4089 retained the exact saved three-geometry pools and common
fragment membership, with complete source/cap/bond checks described in
[REPORT.md](REPORT.md). Endpoint job4091 completed all12 cells. Both saved MACE
runtime fixtures passed with zero energy residual and force residuals below
1e-13eV/Å. Every geometry/metal's standalone xTB6.7.1 versus tblite0.7.0 solvent
component and population gates passed. All four native water-diagonal component
checks passed at every endpoint; the vacuum-charge rounding bounds remained
informative. The actual candidate is OMOL vacuum energy plus native ALPB at
GFN2 vacuum charges, not the old self-consistent solvent subtraction.

Comparison4092 verified66 saved-output algebra checks. Its maximum residual
was6.40e-11kcal/mol against the frozen1e-8 tolerance. There are zero missing
sources or endpoint cells in this declared2-source/12-cell matrix. This is
numerical/implementation qualification, not independent biological validation.

The endpoint job took680.802s on seven allocated CPUs, recording4765.611
allocated CPU-s. Actual receipts contain12 scientific MACE endpoints plus two
qualification evaluations,24 geometry/metal/medium inputs evaluated by both
GFN2 backends (48 SCC calls), and24 fixed-charge native ALPB evaluations. No new
optimization, DFT, fold or GPU computation was used. The separate saved-output
comparison took0.250s elapsed. Allocation cost does not claim measured utilization.

## Decision and limits

Close the narrow hypothesis that this particular omitted local packing dominates
the problematic accommodation response. Further indiscriminate context expansion
is not supported by this result. The already-authorized comparison using real
1FLG source structures is the more useful next discriminator between a modeled
source-geometry problem and behavior persisting in experimental source geometry.

This finite test concerns existing motions in two consumed structures. It does
not exclude other scaffold mechanisms or establish how a fresh search in the
larger context would behave. The additions include peptide and aromatic chemistry
alongside Leu485, so the response cannot be attributed uniquely to steric packing
or one residue. No absolute total energies were compared across context sizes,
no reference bands were transferred to the249-atom context, and no repaired
classification or production promotion is claimed.

Machine-readable exact values, gate evidence and immutable artifact pins:
[FINAL_SUMMARY.json](FINAL_SUMMARY.json). Authoritative comparison:
`workspaces/exaa_source_packing_20261014/packing_comparison_v1/execution/RESULT.json`
(SHA256 `664a35cf2fc117f77acadcf45f5e2a8fa5409a83952f27874e987dafa010ac3c`).
