# Whole peptide reference for the 7CCN deformation

All four native PBE0 calculations and the saved comparison are complete. The reference supports the large energy lowering obtained when the whole peptide can move. MACE overstates the La lowering by 5.81 kcal/mol and the difference between the two metals' lowering by 5.04 kcal/mol. The terminal forces also disagree enough that further MACE optimization needs a local direction check.

## Fixed comparison

This is the same 256-atom, Lu-conditioned 7CCN peptide at its existing start and 80-step La-search endpoint. Both metals use the identical coordinate pair and the same saved fixed-charge ALPB correction. No new geometry, optimization, solvent calculation or MACE call was introduced. The endpoint was admitted by physical checks but did not meet the original force convergence criteria.

The frozen [specification](SPEC.md) and [inventory](INVENTORY.json) preceded the new reference results. This is a second deformation of the same peptide studied in 7CCO, not independent protein validation.

## Deformation work

Values are after minus before, in kcal/mol. Negative values mean lowering along this particular deformation.

| Quantity | La | Dy | Dy minus La |
| --- | ---: | ---: | ---: |
| PBE0 vacuum | -107.812398 | -74.584530 | +33.227868 |
| MACE vacuum | -113.621610 | -75.352475 | +38.269135 |
| Shared fixed-charge ALPB | -24.025982 | -28.154146 | -4.128164 |
| PBE0 plus ALPB | **-131.838380** | **-102.738676** | **+29.099704** |
| MACE plus ALPB | -137.647592 | -103.506621 | +34.140971 |
| MACE minus PBE0 work error | -5.809212 | -0.767945 | +5.041267 |

Only the Dy work error lies within the original 1 kcal/mol reporting resolution. That resolution is a discrepancy descriptor, not an admission gate, uncertainty estimate or fitted tolerance. The shared solvent term cancels from the model-versus-reference error: these discrepancies are in the vacuum energy response under the two models.

The +29.10 kcal/mol differential means this **La-selected deformation** lowers the La energy more than the Dy energy. It is not a metal preference, binding energy or comparison of two independently equilibrated ensembles.

## Forces and the remaining mechanics question

All values below are kcal/mol/Angstrom. RMS is the RMS of full atom-vector magnitudes over all 256 atoms.

| Composite force measure for La | Before | After |
| --- | ---: | ---: |
| PBE0 plus ALPB maximum | 77.299293 | 12.826316 |
| PBE0 plus ALPB RMS | 22.009411 | 4.204686 |
| MACE plus ALPB maximum | 81.733284 | 9.739698 |
| MACE plus ALPB RMS | 21.424237 | 1.987798 |
| MACE minus PBE0 maximum error | 8.128667 | 7.660828 |
| MACE minus PBE0 RMS error | 3.513854 | 3.574975 |
| PBE0 plus ALPB force along frozen deformation chord | +40.285425 | -1.716526 |
| MACE plus ALPB force along frozen deformation chord | +38.510768 | +1.580140 |

The relaxation greatly reduces the reference residual force, but does not produce a reference minimum. At the final point, the two models exert opposite force components along the frozen before-to-after chord. That chord is not necessarily the optimizer's next direction; its sign alone cannot justify either extending or stopping the actual optimizer step.

Both Dy reference force arrays are available. The corresponding archived MACE requests were energy-only, so both Dy model-force comparisons remain unavailable. No replacement forces were synthesized.

## Decision

Retain the evidence that broad motion addresses a real energetic strain in this represented system. Keep the original force convergence criteria and preserve all incomplete searches. The previous hybrid-model and electronic-base candidate failures remain closed; this panel does not repair those models.

The subsequent [saved direction audit](../lbt_7ccn_direction_audit_v1/REPORT.md) reconstructs the **direction the unchanged optimizer would produce from the existing 7CCN La terminal checkpoint if continuation were permitted**. Both local energy derivatives are negative: MACE -0.647635 and PBE0 -0.744240 kcal/mol per dimensionless step. The completed optimizer stopped at its 80-step budget before proposing this direction. The audit uses saved data only and chooses no alternative direction. It supports preparing one separately declared 12-step continuation scout; it does not establish full force accuracy or guarantee descent over a finite step.

## Native admission and limitations

The unchanged Slurm collector re-admits all four raw native outputs, task/configuration identities and full 256-atom gradients. Native convergence took 15, 18, 17 and 17 cycles for before-La, after-La, before-Dy and after-Dy respectively. Printed residual criteria do not all pass; the original native convergence semantics and all flags remain preserved. After-La reaches native gradient-check convergence with only its orbital-gradient residual within the printed threshold. Small-gap warnings remain.

The vacuum dianion, conditional free termini, frozen-f Dy representation and post-added solvent correction remain the model assumptions. SCF convergence does not establish electron binding or ground-state stability. Analytic gradient admission is not independent derivative qualification. Physical-state and force-consistency qualification remain false. No affinity, natural LanM transfer or production promotion follows.

## Execution and preservation

| Job | Work | Elapsed seconds | CPUs | Allocated CPU seconds |
| --- | --- | ---: | ---: | ---: |
| 14793 | Saved preflight | 3 | 1 | 3 |
| 14795 | Before La native reference | 23,440 | 12 | 281,280 |
| 14880 | After La native reference | 26,722 | 12 | 320,664 |
| 14881 | Before Dy native reference | 23,858 | 12 | 286,296 |
| 14882 | After Dy native reference | 24,232 | 12 | 290,784 |
| 15173 | Saved four-cell comparison | 8 | 1 | 8 |
| **Total** | Four molecular calls | | | **1,179,035** |

These are allocation costs, not measured CPU utilization. Collector 15173 performs zero new molecular calls. Result `fff95c2f235947acb24e91991915b320731992f48c2b3d6c9ca3912a6e927f38`, complete force tables, configuration, submission, code and accounting are preserved in [full_result_v1](full_result_v1/EVIDENCE_PINS.json). Native endpoint evidence is retained in the scout and endpoint review directories. Original source results and failed branches remain intact.
