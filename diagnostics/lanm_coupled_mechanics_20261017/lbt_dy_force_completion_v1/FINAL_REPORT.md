# Most of the force discrepancy is shared between La and Dy

20 October 2026. Job15720 completes the three missing Dy force requests and
all eight endpoint comparisons. The MACE–PBE0 discrepancy is substantially
smaller after taking the matched Dy-minus-La force difference. Cancellation
is incomplete and is not a validation of metal preference.

## Complete fixed panel

RMS means the root mean square of atom-vector magnitudes over all256 atoms;
all force values are kcal/mol/Å. Error means MACE minus PBE0. Both models use
the same saved per-metal ALPB term, which cancels from their error.

| Source and geometry | La error RMS | Dy error RMS | Dy−La error RMS | Dy−La maximum atom error |
| --- | ---: | ---: | ---: | ---: |
| 7CCO before | 3.511987 | 3.529411 | 0.501500 | 2.704859 |
| 7CCO after | 3.516438 | 3.531992 | 0.476265 | 2.259919 |
| 7CCN before | 3.513854 | 3.549831 | 0.726763 | 6.242802 |
| 7CCN after | 3.574975 | 3.609827 | 0.563471 | 4.148825 |

The matched error is computed from full vectors:
`(F_Dy − F_La)_MACE − (F_Dy − F_La)_PBE0`.
It is not the difference between two RMS values.

For context, the reference composite Dy−La force RMS is2.411517/2.435107
for7CCO before/after and4.563925/2.187166 for7CCN. Thus the remaining
metal-dependent discrepancy is smaller than the reported reference contrast,
but cannot be treated as negligible or converted into an affinity error.
Full per-atom vectors, maxima and original frozen projections are retained.

The close Dy deformation-energy agreement previously reported does not imply
close absolute forces: all four Dy force-error RMS values remain about3.5–3.6.
Conversely, using the full absolute error as the uncertainty in the metal
contrast would miss the substantial cancellation actually observed here.

## Execution and provenance

All eight saved-input checks passed before model initialization. Exactly three
native MACE force calls were made, on7CCO/before_Dy and7CCN/before_Dy/after_Dy.
Every new energy reproduces its original saved energy exactly (reported
delta0.0eV), passing the original1e−7eV replay gate. Five MACE force arrays and
all eight DFT/ALPB gradient sets were reused. No original POINT was overwritten.

The job used one warm H200 model, loading in6.749189s, and took10.067909s
application wall time. Slurm records12s on12CPUs and one GPU:
**144 allocated CPU-seconds and12 allocated GPU-seconds**. Allocation cost is
separate from measured CPU utilization. No new DFT, ALPB, GFN, MM, geometry or
optimization calls occurred.

Result SHA256:
`f42cffcf5717a988ed569db925b840f02dfb8bd3b4c04750b77e33bc231c4790`.
Original execution is under
`workspaces/lanm_coupled_mechanics_20261017/lbt_dy_force_completion_v1/run/`.
Source/config are `d77cc40e`/`b1655fb8`; inventory is `e2c196d1`.
The separate raw supplements and their task identities are linked directly
from each force report. Original energy/work comparisons remain unchanged.
The [complete execution archive](actual_result_v1/ARCHIVE.json), SHA`e1eecdbe`,
preserves raw files and receipts with verified copy hashes.
[Independent actual review](ACTUAL_REVIEW.md), SHA`3d3d36dd`, finds no
admission/provenance blocker and retains all qualification limits.

## Interpretation and next decision

This completes the declared force-coverage test. The subsequent
[saved-only localization](../lbt_force_error_localization_v1/FINAL_REPORT.md)
completed as15797: shared discrepancy is broadly distributed across peptide
heavy atoms, while the remaining metal-dependent discrepancy is concentrated
more strongly in7CCN side-chain oxygens. That is an error-distribution
diagnostic, not an attribution to a unique interaction or a fitted correction.
The separately running step725 La reference15545 still controls conclusions
about reference forces at the model-converged endpoint.

The pairs are two consumed deformations of the same peptide:7CCO follows the
hybrid-La motion,7CCN the original80-step whole-MACE-La motion. They are not
independently equilibrated metal states or independent proteins. PhysicalDy6
MACE and effective-singlet frozen-4f9 DFT remain different approximate
representations. Native SCF, derivative/state, vacuum-dianion, terminal chemistry,
ECP/no-SOC and post-added solvent limitations remain. No independent solvent
validation, numerical force-accuracy qualification, natural transfer, production
promotion or preference claim follows.
