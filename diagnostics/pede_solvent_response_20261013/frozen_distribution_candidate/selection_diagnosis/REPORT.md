# Changed solvent response, not a new geometric failure

Saved-matrix analysis4016 completed on one Slurm CPU in one allocation-second.
It retains all278 records × three policies, including unavailable inputs. All
available reconstructions close within the predeclared1e−8kcal/mol tolerance.
No electronic, native-solvent, folding or optimization call and no recalibration.

**The new abstentions are driven by the changed energy model plus movement of
the frozen Ca boundary. Geometry reselection does not cause these regressions.**
This is an arithmetic attribution, not proof of a unique electronic defect.

The candidate's canonical Ca upper boundary moves by−4.398647kcal/mol, making
the correct-Ca criterion stricter on this raw scale. The La lower boundary moves
by+0.579505kcal/mol. Both follow the previously frozen rule; neither is adjusted
using these transfer outcomes.

For the five newly inconclusive La-conditioned PedE/ExaA sources, changing the
energy expression at the old selected geometries shifts Ca−La upward by
1.9585 to5.6490kcal/mol. Four keep exactly the old selections. ExaA sample2
reselects a lower-Ca-energy candidate, reducing its contrast by1.007315kcal/mol;
this helps but does not restore the correct-Ca band. Two-Ca ExaA sample3 also
keeps its old geometries; its model shift is+1.886036kcal/mol. Therefore rerunning
the same optimizer or repairing a newly changed geometry is not a supported
explanation for these abstentions.

## PedE corrections are not a favorable threshold shift

Both corrections overcome the stricter Ca boundary:

| Two-Ca PedE source | Model effect at old selections | Reselection effect | New correct-Ca margin |
|---|---:|---:|---:|
| sample1 | −15.455739 | 0 | +3.633224 |
| sample2 | −7.930920 | −0.449340 | +0.232757 |

All values are kcal/mol on the declared contrast scale. Positive correct-class
margin reaches that class's band; it is not a confidence interval or affinity.
Sample2 would remain inconclusive if the old selections were held fixed under
the new model. Its correction combines the changed solvent contribution and
reselection from the same already-computed common pool. Its small margin merits
caution; no new basin search or selective source choice was performed.

RepairedMMOL2048 keeps its selections and gains+1.481323kcal/mol in contrast,
against the+0.579505 movement of the La boundary. A0ACCa4 also keeps its selections
and moves+5.783719, improving to inconclusive without reaching the La band.

## Implication

The candidate gives genuine PedE corrections while losing some Ca-control
decisions through combined energy and calibration changes. Retuning the bands
on those controls, selecting the favorable solvent model per protein, or adding
a fitted response fraction would conceal this tradeoff. None is implemented.

The next physical question is how much solute electronic response is reliable
across these different pockets. The current results do not justify a new fold
campaign or broader geometry search as the immediate remedy. Preserve both
versions as separately named development evidence while selecting the next
response model; no blanket production promotion.

Exact result: `workspaces/pede_solvent_response_20261013/selection_diagnosis_v1/execution/RESULT.json`,
SHA256`f746c4473e86a84019318aed92a2a0f2508a41f6ec6dfbb89ffb2cbd33354cb6`.
The adjacent`ROWS.tsv` contains every model/selection/band term, counterfactual
decision and closure residual. `SUBMISSION.json`, `EXECUTION.json` and root
acknowledgement preserve actual commands, timing and completion wake.
