# Fresh W9 radial scoring contract

Pure saved-data `evaluate(data, model, body, anchors)`; no molecular entry point.
Four fixed slots: aquo9_minus_La, aquo9_minus_Dy, aquo9_plus_La, aquo9_plus_Dy.
`data['rows'][gid]` contains exact geometry, input graph, bodies and optional
reference energy/forces. `model[gid]` is an energy/force record or unavailable.
All energies are kcal/mol; forces are kcal/mol/Å. `anchors['model']` and
`anchors['reference']` are keyed La/Dy and contain the exact central W9 E/F.

Individual work is displaced minus central energy within the same metal and
same method. Paired work is Dy work minus La work. Error is model minus reference.
All four individual and two paired work errors must have absolute value ≤1.
The unchanged body helper supplies category-balanced net forces and centroid
torques divided by 1Å. Its RMS is `2*norm(weighted)` and maximum is the largest
absolute raw body component. Each of the six new rigid checks must pass RMS≤2
and maximum≤5. All twelve comparisons are required; no aggregate replacement,
paired cancellation rescue or outcome-dependent threshold exists.

An unavailable cell/anchor leaves its dependent comparisons explicitly unavailable.
Malformed available values are rejected. No missing value becomes zero. Full
availability and all gates passing are necessary for the validation pass. This
is limited transfer against a provisional reference, with no automatic fitting,
physical-model promotion or production release. The central geometry anchors are
consumed; the four displaced cells remain excluded from the fit.
