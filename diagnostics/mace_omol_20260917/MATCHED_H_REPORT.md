# Matched H preparation: boundary consistency improves; ordering still fails

The normalized-H vacuum hybrid passes all 13 numerical checks and the GGR
partition test, but only **2/4** predeclared affinity-order comparisons.
The all-four accuracy gate fails. Baseline/default and historical results remain
unchanged; no production promotion or broad accuracy claim.

| Alpha structure | GGR core | Alpha minus GGR, kcal-scale | Gate |
|---|---|---:|---|
| 1F6S | Extended | +4.329631810 | Pass |
| 1F6S | Connected | +4.211135326 | Pass |
| 6IP9 | Extended | -0.805098984 | Fail |
| 6IP9 | Connected | -0.923595467 | Fail |

The fixed criterion remains each difference >0.02. Both alpha structures and
both GGR representations remain in the denominator. These are two consumed
biological groups with qualified cross-study directions, not four independent
or prospective observations.

## Transferable finding

The same already established source-H normalization was applied to both the
quantum and learned cores, matching the frozen whole-protein inputs. The GGR
connected-minus-extended hybrid shift falls from **1.829806242 to
0.118496483 kcal-scale**. Its normalized components are DFT vacuum
-23.951208420 and learned core -24.069704903; the common whole term cancels.
This improves consistency between these two representations. It does not
establish general partition invariance or aqueous affinity accuracy.

The two alpha structures still differ by 5.134730793 kcal-scale. A source
structure effect now exceeds the tested GGR partition sensitivity. Neither
structure is discarded. Normalization improves one alpha/GGR comparison and
worsens the other; favorable classifications were not used to choose H geometry.

Protocol: `masked_omol_matched_normalized_H_vacuum_hybrid_v1`.

```
H_M = E_DFT,vacuum(normalized_core,M)
      + T_mask(normalized_full,M) - T_mask(normalized_core,M)
R_H = H_Ca - H_La
```

This expression has no aqueous reaction field, compatible aquo reference or
calibrated absolute zero. It remains a mixed quantum/learned descriptor.
Relaxation, entropy and predictive-improvement fields remain unavailable.
The archived original-H failure is a separate immutable result.

## Preparation and execution

Transferred 21/47/16/18 existing normalized source H coordinates into
GGR extended/connected and alpha 1F6S/6IP9. Every heavy coordinate and synthetic
sigma cap stayed exactly unchanged. Donors, charge, multiplicity, protonation,
assembly and water inventory stayed fixed. Source mapping, charge parity,
paired coordinates and actual whole-result cache receipts pass. All duplicate
whole matches agree; selection was lexical, independent of energy ordering.

Eight native ORCA 6.1.1 vacuum r2SCAN-3c analytic-gradient endpoints completed
job1200950, with eight matching masked-OMOL core energy calls in job1200951.
Six exact normalized whole results were reused. No failed scientific call,
new whole inference, numerical gradient, optimization, solver or training.
Gradients are retained; no old-H coordinate Jacobian is applied to them.

Measured allocation costs: **751 seconds on 64 CPUs** for DFT, 48,064 allocated
core-seconds and 43,768 reported CPU-seconds; **67 seconds on one A5000/16 CPUs**
for MACE, 1,072 allocated core-seconds and 83.861 reported CPU-seconds.
MACE model evaluation totals 3.856576160 seconds; peak GPU allocation
802,472,448 bytes. Slurm peak batch RSS: DFT 15,534,988 KiB, MACE 434,244 KiB.
Local preparation/dry-run/report/test resource receipts are recorded separately.
Hardware differs from the previous DFT job; this is not a matched speedup test.

Five distinct new real-fixture tests and six existing native regressions pass.
The actual-output test was explicitly skipped before execution, then passed in
3.973 seconds against the completed DFT and MACE receipts. Four preparation
checks took 5.157 seconds; native regressions took 12.286 seconds. No scientific
executable was replaced with synthetic successful output.

## Artifacts and recommendation

[Commands](MATCHED_H_COMMANDS.md), [declared plan](MATCHED_H_NORMALIZATION_PLAN.md),
[compact result](MATCHED_H_RESULT.json). Under `workspaces/mace_omol_20260917/`:
`matched_H_prepared_v1`, `matched_H_quantum_v1`, `matched_H_mace_v1`,
`matched_H_report_v1/result.json`, `matched_H_cost_v1.json` and terminal receipts.
The model-engineering subtotal V16 includes eight new MACE calls; DFT costs
remain separate and are not counted as model forwards.

**Retain the baseline.** The new hybrid is affordable on these cores and passes
this numerical/partition check, but its predictive usefulness is insufficient.
The broader MACE goal remains active. The vacuum model's missing solvent and
the persistent source-structure sensitivity are distinct unresolved issues;
this result does not identify either as the unique cause of the failure.
