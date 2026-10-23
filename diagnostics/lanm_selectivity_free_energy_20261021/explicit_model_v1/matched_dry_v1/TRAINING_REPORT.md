# Saved training admission and controls

**All 28 training cells are admitted and replayed; the Dy fit remains on hold.** Preparation used one CPU thread for6.03s and made zero molecular calls. It checked the exact saved AMOEBA native energy/gradient outputs, original ECP-PBE0 parser results, input decks, coordinates, states, atom types and covalent graphs. The self-target objective is zero and a central equal/opposite internal pair force has zero rigid response. Reserved validation outcomes were not opened.

Dataset: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/matched_dry_v1/training_v1/RESULT.json`, SHA256 `99ea27f86871b505dc622c1391a6321381f7430d3ad4a86a72e34ab9a629f6bb`. The dataset contains all source/raw pins, arrays, fourteen matched pairs, complete connected-body definitions, printed native stopping evidence, objective values and unavailable-qualification flags. Source and configuration are pinned inside it.

## Fixed La already fails the independent control

| Control | Category-balanced rigid RMS | Maximum body component error | Gate |
| --- | ---: | ---: | --- |
| La water dimer+aquo8 | 6.163 | 16.143 | Fail |
| La donor fixtures | 8.776 | 39.674 | Fail |
| Initial Dy water | 1.462 | 6.156 | Fail maximum |
| Initial Dy donors | 10.999 | 29.233 | Fail |

Force units are kcal/mol/Å; body torques are divided by the declared1Å scale. Frozen thresholds were RMS≤2 and maximum≤5. Five of six La hydrated radial-work errors exceed the separate1kcal/mol threshold. A Dy-only fit cannot repair this fixed La control or qualify a paired potential. The root decision therefore holds fitting while the La/reference discrepancy is diagnosed.

## Coarse reference derivative screen

The exact moving coordinates translate the entire ligand at ±0.15Å, with metal and waters fixed. Central analytic slopes are minus the sum of ligand forces projected onto that translation axis. All six individual and three paired comparisons are retained.

| Family | La analytic / secant | Dy analytic / secant | Dy−La analytic / secant |
| --- | ---: | ---: | ---: |
| Mono-acetate | 45.994 /44.199 | 63.064 /61.878 | 17.070 /17.680 |
| Bi-acetate | 37.956 /34.921 | 72.341 /70.202 | 34.385 /35.281 |
| Acetamide | 14.853 /13.139 | 31.729 /30.655 | 16.877 /17.516 |

Eight of nine pass the frozen coarse screen. Acetamide-La differs by1.713kcal/mol/Å against its1.485 screen threshold. Its endpoint analytic slopes are −5.024 and24.284, so curvature over this large displacement matters. The flag remains unchanged: this is neither small-step derivative qualification nor evidence that the analytic force is defective. No row was removed, no threshold changed, and no new displacement was evaluated.

## Interpretation

Targets are consistently **uncorrected** ECP-PBE0 energies and gradients. CP/CC results remain separate method-sensitivity evidence. The published La calibration reference differs; native admission and exact replay do not establish the physical accuracy of either model. Full atomic forces remain available, but only net body force/torque enters the proposed loss, preventing stiff shared intramolecular distortions from dominating an ion-parameter fit. Existing24 reserved validation cells remain geometry-only and unopened. No production sampling or parameter promotion follows.
