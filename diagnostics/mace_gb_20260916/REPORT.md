# Frozen MACE monopole / OBC-II solvent check

Status: complete. Numerical checks pass: True.

Baseline unchanged. Consumed 1H4I development inputs; zero new MACE/DFT evaluations.
G is GB reaction energy only. No duplicate direct Coulomb or CPCM contribution.
Raw descriptor R = MACE vacuum R + G(Ca) - G(La); all quantities below in kcal/mol.

```json
{
  "direct": {
    "large": {
      "GB_Ca_minus_La_kcal_mol": -901.5888791488014,
      "MACE_vacuum_R_kcal_mol": -404414.35789542325,
      "S_kcal_mol": null,
      "calibrated_class": null,
      "descriptor_R_kcal_mol": -405315.9467745721
    },
    "medium": {
      "GB_Ca_minus_La_kcal_mol": -84.0550594094293,
      "MACE_vacuum_R_kcal_mol": -405338.2747411182,
      "S_kcal_mol": null,
      "calibrated_class": null,
      "descriptor_R_kcal_mol": -405422.3298005276
    }
  },
  "checkpoint_disagreement": {
    "GB_Ca_minus_La_kcal_mol": -817.5338197393721,
    "MACE_vacuum_R_kcal_mol": 923.916845694941,
    "descriptor_R_kcal_mol": 106.38302595552523
  }
}
```

## Numerical checks

```json
[
  {
    "energy_difference_kcal_mol": 3.9818590408913224e-05,
    "force_difference_max_eV_A": 4.3367611946010243e-07,
    "name": "native_custom_La",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 4.603263630542642e-05,
    "force_difference_max_eV_A": 3.1398410382799113e-07,
    "name": "native_custom_Ca",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 7.275957614183426e-12,
    "force_difference_max_eV_A": 2.9021878270810847e-11,
    "name": "medium_La_rotation",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_kcal_mol": 0.0,
    "force_max_eV_A": 0.0,
    "name": "medium_La_identity",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 2.7284841053187847e-12,
    "force_difference_max_eV_A": 2.6289071614060688e-11,
    "name": "medium_Ca_rotation",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_kcal_mol": 0.0,
    "force_max_eV_A": 0.0,
    "name": "medium_Ca_identity",
    "pass": true,
    "status": "computed"
  },
  {
    "difference_kcal_mol": 4.547473508864641e-12,
    "name": "medium_contrast_rotation",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 3.092281986027956e-11,
    "force_difference_max_eV_A": 2.6354285073693084e-11,
    "name": "large_La_rotation",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_kcal_mol": 0.0,
    "force_max_eV_A": 0.0,
    "name": "large_La_identity",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 2.9103830456733704e-11,
    "force_difference_max_eV_A": 2.8800423851338408e-11,
    "name": "large_Ca_rotation",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_kcal_mol": 0.0,
    "force_max_eV_A": 0.0,
    "name": "large_Ca_identity",
    "pass": true,
    "status": "computed"
  },
  {
    "difference_kcal_mol": 1.8189894035458565e-12,
    "name": "large_contrast_rotation",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 0.0,
    "force_difference_max_eV_A": 0.0,
    "name": "medium_La_repeat",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 0.0,
    "force_difference_max_eV_A": 2.413347299778934e-13,
    "name": "medium_La_translate",
    "pass": true,
    "status": "computed"
  },
  {
    "energy_difference_kcal_mol": 3.50145537595381e-05,
    "force_difference_max_eV_A": 1.4412916893302707e-07,
    "name": "medium_La_reference",
    "pass": true,
    "status": "computed"
  }
]
```

## Limits

Frozen monopoles omit dipolar and self-consistent solvent response. Metal cavity factors are unvalidated.
GB forces hold charges fixed; they are not gradients of the combined descriptor.
No compatible aquo reference, calibrated S/class, or predictive-accuracy claim.
Recorded attempts: 19; successful tasks: 19/19.
Whole-job costs are recorded separately from evaluation timings.

## Measured cost and judgment

Successful job 1200700: 19/19 calculations, 46 allocated GPU seconds,
736 allocated core-seconds, 51.460 actual CPU seconds. Including failed
CUDA13 context job1200695: 54 GPU seconds, 864 allocated core-seconds,
64.419 actual CPU seconds, 23 completed solver calls and one failed context.
The isolated CUDA12 environment fixed execution without a physical-model change.
Install/preparation wall and CPU time were not fully instrumented; they are
one-time engineering work and are not included in those allocation totals.

Warm full GPU evaluations were about0.09s; full Reference evaluation4.903s.
GPU context/process initialization dominates. Maximum worker RSS303916KiB;
GPU memory was not measured by this OpenMM worker. No MACE/DFT reevaluation.

Numerical credibility: all declared checks pass. Scientific informativeness:
medium/large raw R disagreement falls923.916846 to106.383026kcal/mol; solvent
accounts for a substantial difference but residual model dependence remains
large. This is not evidence of a correct biological classification. Affordability:
the solver contribution is inexpensive next to the measured MACE endpoint pair.
Keep the baseline; pursue the solvent descriptor on multiple real proteins.

Rerun/collect commands are in COMMANDS.md. Scientific limitations above remain
in force, including unvalidated metal cavity and frozen monopole response.
