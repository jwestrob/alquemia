# Fixed scaffold transition audit

Four fixed reference cells; both metals use the same La-before and La-proposed after coordinates.
All energies below are after minus before (kcal/mol); forces are kcal/mol/Angstrom.
MM interaction energies are whole-interaction contributions, not residue energies or quantum bond energies.

## MM work

| Group/class | La | Dy |
|---|---:|---:|
| retained_bonded | -10.682960 | -10.682960 |
| MM_LJ | -6.180039 | -6.180039 |
| MM_Coulomb | -3.189036 | -3.189036 |
| QM_MM_LJ | -3.229524 | -3.228347 |
| total | -23.281560 | -23.280382 |
| HarmonicBondForce | -1.406824 | -1.406824 |
| HarmonicAngleForce | -2.981914 | -2.981914 |
| PeriodicTorsionForce | -3.874648 | -3.874648 |
| CMAPTorsionForce | -2.419574 | -2.419574 |

## Terminal-touching interaction work

The union counts every interaction once. Per-terminus sets can overlap; overlap is listed explicitly.

| Set / class | Count | La | Dy |
|---|---:|---:|---:|
| union / total | 187 | -9.532816 | -9.532816 |
| union / HarmonicBondForce | 34 | -1.777162 | -1.777162 |
| union / HarmonicAngleForce | 61 | -2.035376 | -2.035376 |
| union / PeriodicTorsionForce | 90 | -3.173561 | -3.173561 |
| union / CMAPTorsionForce | 2 | -2.546718 | -2.546718 |
| Phe1 / total | 133 | -7.378008 | -7.378008 |
| Phe1 / HarmonicBondForce | 23 | -1.413262 | -1.413262 |
| Phe1 / HarmonicAngleForce | 41 | -2.290653 | -2.290653 |
| Phe1 / PeriodicTorsionForce | 68 | -2.335136 | -2.335136 |
| Phe1 / CMAPTorsionForce | 1 | -1.338957 | -1.338957 |
| Ala17 / total | 54 | -2.154808 | -2.154808 |
| Ala17 / HarmonicBondForce | 11 | -0.363900 | -0.363900 |
| Ala17 / HarmonicAngleForce | 20 | 0.255277 | 0.255277 |
| Ala17 / PeriodicTorsionForce | 22 | -0.838425 | -0.838425 |
| Ala17 / CMAPTorsionForce | 1 | -1.207760 | -1.207760 |
| both_termini_overlap / total | 0 | 0.000000 | 0.000000 |
| both_termini_overlap / HarmonicBondForce | 0 | 0.000000 | 0.000000 |
| both_termini_overlap / HarmonicAngleForce | 0 | 0.000000 | 0.000000 |
| both_termini_overlap / PeriodicTorsionForce | 0 | 0.000000 | 0.000000 |
| both_termini_overlap / CMAPTorsionForce | 0 | 0.000000 | 0.000000 |

## Fifteen largest absolute interaction works

Descriptive ranking by the larger absolute La/Dy work; all 187 interactions remain in RESULT.json.

| Class / native index | Atoms | La | Dy |
|---|---|---:|---:|
| CMAPTorsionForce / 0 | A/1//C, A/2//N, A/2//CA, A/2//C, A/2//N, A/2//CA, A/2//C, A/3//N | -1.338957 | -1.338957 |
| HarmonicAngleForce / 0 | A/1//N, A/1//CA, A/1//C | -1.241135 | -1.241135 |
| CMAPTorsionForce / 14 | A/15//C, A/16//N, A/16//CA, A/16//C, A/16//N, A/16//CA, A/16//C, A/17//N | -1.207760 | -1.207760 |
| HarmonicBondForce / 4 | A/1//CA, A/1//N | -0.632817 | -0.632817 |
| PeriodicTorsionForce / 587 | A/16//C, A/17//N, A/17//CA, A/17//CB | -0.544007 | -0.544007 |
| HarmonicAngleForce / 6 | A/1//H3, A/1//N, A/1//CA | -0.531561 | -0.531561 |
| HarmonicBondForce / 21 | A/1//H3, A/1//N | -0.524121 | -0.524121 |
| PeriodicTorsionForce / 586 | A/16//C, A/17//N, A/17//CA, A/17//CB | 0.504824 | 0.504824 |
| HarmonicBondForce / 247 | A/17//C, A/17//CA | -0.399335 | -0.399335 |
| PeriodicTorsionForce / 681 | A/16//C, A/17//CA, A/17//N, A/17//H | -0.300050 | -0.300050 |
| PeriodicTorsionForce / 568 | A/16//CA, A/16//C, A/17//N, A/17//H | -0.293354 | -0.293354 |
| PeriodicTorsionForce / 10 | A/1//H3, A/1//N, A/1//CA, A/1//CB | -0.279695 | -0.279695 |
| PeriodicTorsionForce / 7 | A/1//H2, A/1//N, A/1//CA, A/1//CB | -0.278055 | -0.278055 |
| PeriodicTorsionForce / 28 | A/1//CB, A/1//CA, A/1//N, A/1//H | -0.277615 | -0.277615 |
| PeriodicTorsionForce / 11 | A/1//H3, A/1//N, A/1//CA, A/1//HA | -0.275588 | -0.275588 |

## Force offsets and changes

Errors are model minus PBE0 plus the identical fixed-charge ALPB. Change is after minus before.

| Metal / model / quantity | Full256 vector RMS | Phe1-direction projection | Deformation-chord projection |
|---|---:|---:|---:|
| La / reference_PBE0_plus_fixed_ALPB / force_before | 3.666509 | 3.341027 | 3.638430 |
| La / reference_PBE0_plus_fixed_ALPB / force_after | 8.921877 | -40.851355 | -2.572905 |
| La / reference_PBE0_plus_fixed_ALPB / force_change | 7.815341 | -44.192382 | -6.211335 |
| La / whole_MACE_composite / force_before | 0.921311 | 0.798547 | 2.466204 |
| La / whole_MACE_composite / error_before | 3.511987 | -2.542480 | -1.172226 |
| La / whole_MACE_composite / force_after | 7.773379 | -44.464198 | -1.076003 |
| La / whole_MACE_composite / error_after | 3.516438 | -3.612843 | 1.496902 |
| La / whole_MACE_composite / force_change | 7.831780 | -45.262745 | -3.542207 |
| La / whole_MACE_composite / error_change | 0.530091 | -1.070363 | 2.669128 |
| La / hybrid_composite / force_before | 9.121235 | 51.063006 | 19.928221 |
| La / hybrid_composite / error_before | 10.141916 | 47.721979 | 16.289791 |
| La / hybrid_composite / force_after | 0.901487 | -1.895426 | -0.329568 |
| La / hybrid_composite / error_after | 8.850041 | 38.955928 | 2.243337 |
| La / hybrid_composite / force_change | 9.255258 | -52.958432 | -20.257789 |
| La / hybrid_composite / error_change | 4.463637 | -8.766050 | -14.046454 |
| Dy / reference_PBE0_plus_fixed_ALPB / force_before | 4.401643 | 3.345638 | 3.873822 |
| Dy / reference_PBE0_plus_fixed_ALPB / force_after | 9.228823 | -40.850094 | -2.518261 |
| Dy / reference_PBE0_plus_fixed_ALPB / force_change | 7.821246 | -44.195732 | -6.392084 |
| Dy / whole_MACE_composite / force_after | 8.227376 | -45.113217 | -1.120943 |
| Dy / whole_MACE_composite / error_after | 3.531992 | -4.263123 | 1.397318 |
| Dy / whole_MACE_composite / transition unavailable | — | — | — |
| Dy / hybrid_composite / force_before | 9.611968 | 51.064870 | 20.194774 |
| Dy / hybrid_composite / error_before | 10.162397 | 47.719232 | 16.320952 |
| Dy / hybrid_composite / transition unavailable | — | — | — |

Before-Dy whole-MACE forces and after-Dy hybrid forces/work remain unavailable. After-Dy MM is a component-only diagnostic.
Two endpoints do not measure a Hessian or separate equilibrium offsets from stiffness uniquely. No parameter fit, accuracy, affinity or preference claim follows.
Inherited limitations: conditional terminal chemistry, vacuum dianion, frozen-f representation, SCF evidence, nonstationary endpoints, and no independent numerical qualification of these whole-reference forces.

Actual MM diagnostics: 36 contexts; 816 energy/full-force queries.
