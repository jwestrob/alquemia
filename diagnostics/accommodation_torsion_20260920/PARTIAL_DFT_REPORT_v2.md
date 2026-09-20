# Partial native adjudication: first large PLM Ca response

**Native DFT confirms the large Ca energy response to the compressed PLM Asp
motion.** This is new physical evidence beyond the small 4MAE control. It does
not yet establish the metal-selectivity shift, because the corresponding La
energies remain pending. No new calculation, retry, method or selection change
was introduced to obtain this snapshot.

Source: `PQQSEQ_83440678cbbd658047c9`, original sample0 expanded context. Positive
extra-Asp chi2 rotation relieves the pre-existing short metal–oxygen contact;
negative rotation compresses it further. Coordinates outside the prescribed
carboxylate motion, charges, water inventory and internal carboxylate geometry
remain fixed.

| Angle | Native DFT Ca work | Native OMOL Ca work | Solvent-transfer work | Composite Ca work | Composite − DFT |
|---|---:|---:|---:|---:|---:|
| −0.2 rad | +30.4742 | +31.0748 | −0.0843 | +30.9906 | +0.5164 |
| +0.2 rad | −22.3929 | −23.3049 | +0.2861 | −23.0188 | −0.6259 |

Units: kcal/mol relative to each model's same-source Ca origin. Native DFT uses
the predeclared r2SCAN-3c/CPCM(Water)/DefGrid3 expanded context. Both directions
and the large magnitude agree. Most of the useful response is already in OMOL;
the solvent term modestly improves agreement in both directions. Thus the Ca
compression-relief signal is not merely an artifact of the cheap model.

This confirms only two Ca displacement works in one previously inspected PLM
structure. No biological label is inferred. The hypothesized stronger La
stabilization, its differential magnitude and any classifier benefit remain
unresolved. No minimum, angle extension, new bands or thermodynamic correction
is established. The much smaller 4MAE differential-work errors (0.137–0.139
kcal/mol) reported in partial_v1 are unchanged.

## Actual completion and failure state

- **7/16 native endpoints complete:** the four 4MAE displacements, plus this PLM
  Ca origin and its two predeclared displacements. All seven receipts confirm
  normal termination and SCF convergence.
- This yields **6/12 available displacement-work comparisons**: four control
  and two PLM Ca. All nine pending native endpoint energies remain null.
- Four endpoints have active outputs and five have not started. No native DFT
  failure marker is observed. La's active SCF output is not treated as a result.
- The pre-existing two low-level SCF failures at the unrelated 1H4I La Glu point
  remain preserved and unavailable. No retry or result substitution occurred.
- There is no execution issue requiring intervention based on this snapshot.
  Final allocation cost remains pending while the existing job runs.

## Immutable artifacts

- `prepared_v3/dft/partial_collection_1203771_v2.json` under
  `workspaces/accommodation_torsion_20260920/` retains all 16 endpoint statuses.
- `workspaces/accommodation_torsion_20260920/partial_DFT_result_v2.json` is the
  unchanged frozen analyzer's output, with unrounded works and missing La
  selectivity fields.
- [PARTIAL_EXECUTION_v2.json](PARTIAL_EXECUTION_v2.json) records the live execution
  checkpoint; partial_v1 and all original scientific artifacts remain unchanged.

The next evidence is the already-running La adjudication. No additional task
was submitted and production remains unchanged.
