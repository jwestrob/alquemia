# Matched local/global MACE decomposition

Status: complete. Baseline unchanged.

All values are kcal/mol, R = E(Ca) - E(La). Global minus core uses matched corrected H coordinates.
Archived DFT uses its original core geometry; no matched corrected-H DFT or valid hybrid is available.

| Case | Archived DFT R | Archived MACE+GB R | Corrected-H MACE+GB R | Full−core R contribution | H effect |
|---|---:|---:|---:|---:|---:|
| ALPHA_1F6S | -405431.282721 | -405470.387923 | -405465.167715 | -32.805847 | 5.220208 |
| ALPHA_6IP9 | -405434.266180 | -405482.340983 | -405478.291014 | -33.435726 | 4.049969 |
| GGR_1GLG | -405416.723475 | -405459.073003 | -405458.227926 | -18.486083 | 0.845077 |
| PQQ_1H4I | -405412.263104 | -405450.037314 | -405449.041553 | 16.182648 | 0.995761 |
| PQQ_4MAE | -405381.691215 | -405439.354512 | -405437.986752 | -12.500322 | 1.367760 |

## Predeclared differences between sites

```json
[
  {
    "environment_contribution_difference_kcal_mol": -28.68296967697097,
    "higher_expected": "PQQ_4MAE",
    "local_states": {
      "archived": 10.682802344905213,
      "global_H": 11.054801915190183
    },
    "lower_expected": "PQQ_1H4I"
  },
  {
    "environment_contribution_difference_kcal_mol": -14.319763683946803,
    "higher_expected": "ALPHA_1F6S",
    "local_states": {
      "archived": -11.31491991475923,
      "global_H": -6.939788108342327
    },
    "lower_expected": "GGR_1GLG"
  },
  {
    "environment_contribution_difference_kcal_mol": -14.949643190891948,
    "higher_expected": "ALPHA_6IP9",
    "local_states": {
      "archived": -23.267980301287025,
      "global_H": -20.06308794557117
    },
    "lower_expected": "GGR_1GLG"
  }
]
```

Positive local differences support the predeclared qualitative ordering. A negative environment
difference erodes that ordering; it includes changed cavity, direct interactions and learned
charge response together. Components do not establish a unique cause.

The archived DFT versus MACE+GB discrepancy includes both electronic and solvent approximations.
Do not mix archived DFT with corrected-H low-level subtraction and call it a computed hybrid.
All cases are consumed development data. Alpha structures are one observation; PQQ class
association is separate from direct affinity evidence. No threshold, aquo gauge or absolute class.


## Interpretation and actual cost

The matched local MACE+GB PQQ contrast is +11.054802 kcal/mol. The global
contribution changes it by −28.682970, reversing the direction. Alpha−GGR
already has the wrong local ordering (−6.939788/−20.063088); global additions
worsen it by −14.319764/−14.949643 kcal/mol. Thus one failure arises when adding
the tested global contribution; another already exists in the local descriptor.
Neither supports blindly adding this environment correction to the DFT baseline.
The observed hydrogen-preparation effects are0.845–5.220kcal/mol in local MACE+GB.
No inference about the uncomputed DFT hydrogen effect is substituted for data.

Jobs1200717/1200719 completed20MACE+20GB calls, no failures, no new DFT.
Allocated GPU time152+27=179s; allocated core time2864s; actual CPU time
200.719+23.425=224.144s. Four real source/parser/recovery tests passed2.539s.
The original baseline scores and coordinates remain unchanged. New protocols
mace_polar_1m_local_core_H_comparison_v1 and
mace_polar_1m_local_frozen_obc2_H_comparison_v1 remain research descriptors.
No valid hybrid score, calibrated class or structural response correction exists.

Next: inspect the already saved short-range learned energy component as a
separate whole-protein descriptor. This tests whether useful local chemical
information survives without the unstable long-range charge-response terms.
It is not a claim that an energy component is a complete physical binding energy.
