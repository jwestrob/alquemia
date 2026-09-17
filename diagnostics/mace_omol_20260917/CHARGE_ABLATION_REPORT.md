# Charge-feature ablation: five-case development pass

Job 1200828 completed all 42 declared forwards without failures. The modified
model passes core native/batched equivalence, whole-chain reproducibility,
rotation and detached-distance checks, all three fixed relative directions,
and the GGR disconnected-sodium test. No thresholds were fitted on these cases.
All cases were previously consumed; two alpha structures are one qualified
affinity comparison. PQQ functional class is a separate evidence stratum.

| Case | Modified descriptor, model kcal | Native R_coord, kcal/mol |
|---|---:|---:|
| ALPHA_1F6S | 36.729109306 | 99.868407432 |
| ALPHA_6IP9 | 28.099819820 | 90.803132693 |
| GGR_1GLG | 23.987421458 | 62.714844813 |
| PQQ_1H4I | 17.579955566 | 101.614125678 |
| PQQ_4MAE | 61.909134503 | 110.647960039 |

Margins: XoxF−MxaF +44.329179; alpha1F6S−GGR +12.741688;
alpha6IP9−GGR +4.112398 model kcal, all above the frozen +0.02 criterion.
The GGR sodium perturbation changes its descriptor by -1.07384e-8 model kcal
(tolerance 0.1). Charge-feature masking imposes spectator invariance; this is
an implementation consistency check, not independent evidence of electrostatics.

## What changed and what it establishes

The same 100M weights and whole-chain inputs are used, with the raw total-charge
embedding replaced by zeros before the native joint projection. Spin, physical
charge, electron count, protonation, explicit waters, coordinates and all learned
parameters remain unchanged. Actual molecular receipts verify the charge input,
zero feature for every atom, unchanged parameters, and separate model units.
These outputs are learned descriptors, not quantum energies for those physical
charges. They do not supply an affinity zero, a solvent correction or a force.

The useful five-case ordering survives removal of the nonlocal charge feature.
This is narrower than broad validation. The original native whole-chain model's
canonical and spectator failures remain recorded; no claim is made that masking
repairs its physical Hamiltonian. The conditional 25-case calibration is the next
test and must keep all cases and the original acceptance rule.

## Numerical checks and cost

14 core checks and 31 alpha checks pass, including energy bookkeeping within
those subsets; all 42 forwards pass native readout accounting. Checks are not
independent biological observations. Eight software/real-fixture tests pass
57.308 s, including actual masked receipt acceptance and rejection of native or
incomplete qualification as authority for extending the changed model.

The job used 798 GPU/allocation seconds (13m18s), 12768 allocated core-seconds,
894.358 reported actual CPU-seconds, 409.63256069645286 s summed model evaluation,
11536478720 bytes peak allocated GPU memory and 2876296 KiB peak process RSS.
This includes qualification and spectator work. Local component/preparation/
test/report receipts are separate under the same workspace, including failed
component checks and the unexecuted superseded v1 preparation. No extra DFT,
solver, training, gradient, trajectory or environment installation occurred.

Full actual receipts and unrounded values: `charge_ablation_report_v1/result.json`
under `workspaces/mace_omol_20260917/`; compact pins in CHARGE_ABLATION_RESULT.json.
The baseline/default and original reference scales remain unchanged. The goal
remains active. Proceed to the conditional canonical extension; do not promote.
