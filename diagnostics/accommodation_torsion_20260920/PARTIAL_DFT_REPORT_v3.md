# Native La adjudication: large differential accommodation is real

**The first compressed PLM geometry has a large metal-dependent response in
native DFT, and MACE captures its direction and approximate magnitude.** At the
preselected +0.2 rad extra-Asp rotation, DFT lowers La energy by 52.261 kcal/mol
and Ca by 22.393: the Ca−La contrast shifts **+29.868 kcal/mol**. Native OMOL
predicts +28.535; OMOL plus matched GFN2 solvent transfer predicts +26.957.
These are electronic displacement works, not relaxed affinities or biological
metal-specificity labels.

This supports pursuing actual nonlinear donor accommodation. It does **not**
establish a numerical pass for general response accuracy: no such threshold was
frozen. The solvent correction worsens agreement for this differential, and
neither endpoint has a demonstrated minimum in the existing sampled domain.

## Exact question and unchanged physical state

Source: `PQQSEQ_83440678cbbd658047c9`, original sample0, same 190-atom expanded
context and native r2SCAN-3c/CPCM(Water)/DefGrid3 recipe for every DFT endpoint.
The coordinate is terminal extra-Asp307 chi2; the complete carboxylate rotates
with the existing physical mapping. Positive rotation relieves the original
1.766 Å metal–oxygen contact. Coordinates outside the prescribed motion,
internal carboxylate geometry, atoms, waters, charge and protonation stay fixed.
Ca and La use the same nuclear coordinates at each angle. This also clarifies
the shorthand "protein coordinates remain fixed" in the older vault checkpoint:
the explicitly rotated donor is the exception.

For model X and metal M, `W_X,M(q) = E_X,M(q) − E_X,M(0)`.
The reported differential is `ΔR(q) = W_Ca(q) − W_La(q)`; positive means the
displacement stabilizes La more on that model's scale. The unchanged reference
offset cancels. No absolute aquo score or decision threshold is transferred here.

## Endpoint work and component audit

All values below are kcal/mol; unrounded values and source pins are in
[PARTIAL_COMPONENTS_v3.json](PARTIAL_COMPONENTS_v3.json).

| Angle (rad) | Metal | Native DFT | Native OMOL | Solvent-transfer work | Composite | Composite − DFT |
|---|---|---:|---:|---:|---:|---:|
| −0.2 | Ca | +30.474221 | +31.074848 | −0.084270 | +30.990578 | +0.516357 |
| −0.2 | La | +62.632960 | +61.813204 | −2.332100 | +59.481104 | −3.151856 |
| +0.2 | Ca | −22.392882 | −23.304891 | +0.286067 | −23.018825 | −0.625943 |
| +0.2 | La | −52.260708 | −51.840285 | +1.864312 | −49.975972 | +2.284736 |

Both models reproduce all four endpoint-work signs. OMOL endpoint errors are
0.420–0.912 kcal/mol in absolute value. The solvent term improves both Ca works
but degrades both La works. This is an actual limitation of the composite
descriptor here; it cannot be credited as the source of the useful differential
response or silently removed after inspecting the result.

| Angle (rad) | Native DFT ΔR | OMOL ΔR | Solvent ΔR | Composite ΔR | OMOL error | Composite error |
|---|---:|---:|---:|---:|---:|---:|
| −0.2 | −32.158738 | −30.738355 | +2.247830 | −28.490526 | +1.420383 | +3.668213 |
| +0.2 | +29.867826 | +28.535393 | −1.578246 | +26.957148 | −1.332433 | −2.910678 |

Relative to the actual native differential signal, OMOL errors are **4.42% and
4.46%**; composite errors are **11.41% and 9.75%**, respectively. Both preserve
the direction and most of the signal. Errors of roughly 3 kcal/mol would still
matter near a decision band; large-signal agreement is not sub-kcal accuracy.
The much smaller 4MAE control differential errors (0.137–0.139 kcal/mol) remain
unchanged. The two cases do not establish uniform accuracy over geometries.

## Decision and evidence against premature relaxation scoring

The conditional prerequisite in
[NEXT_PHYSICAL_DECISION.md](../accommodation_goal_20260920/NEXT_PHYSICAL_DECISION.md)
is met **for pursuing and testing bounded nonlinear accommodation**: both metal
surfaces respond strongly in native DFT, and the cheap models distinguish their
response with the correct sign and approximate magnitude. It is not met for
claiming a validated relaxation correction or improved classifier accuracy.

The existing evidence imposes specific limits:

- Both composite endpoint profiles still have their lowest sampled energy at
  **+0.4 rad**, the positive boundary. No stationary basin or optimum exists in
  the evidence. Native DFT has checked only ±0.2 rad, not ±0.4 or a minimum.
- Composite central curvature changes from **199.29 to 127.37** kcal/mol/rad²
  for Ca and **237.63 to 118.63** for La between the 0.2 and 0.4 rad intervals.
  A quadratic extrapolation would rely on a stiffness that is not constant over
  even the current sampled range.
- Only the terminal extra-Asp mode is adjudicated on this PLM geometry.
  Coupling to anchor Glu, backbone or other donors is untested. The earlier
  unrelated Glu-point SCF failures remain visible and unavailable.
- The useful native response mainly comes from OMOL; the matched solvent
  correction has a measurable differential error here. This argues for retaining
  component reporting and native checks at proposed stationary geometries, not
  assuming that a more elaborate model improves every physical quantity.
- Neither PLM target has a known specificity label. The one geometry-selected
  known-reference grid test weakened a native OMOL false-Ca call to indeterminate
  but added no correct decisive call; the composite already called it La.
  A common preparation/response rule must earn utility on the known-label
  structural variants before broad scanner claims.

A separately frozen nonlinear experiment can now test whether bounded stationary
behavior exists and whether one common rule helps those known-label inputs.
No new angle, optimizer, band or scientific calculation is introduced here.
Current baseline and all prior experimental records remain unchanged.

## Actual execution and checks

This immutable partial snapshot contains **12/16 completed native endpoints**:
four 4MAE displacements; all six first-PLM endpoints; and the second PLM Ca
origin/−0.2 displacement. There are **9/12 available displacement comparisons**.
The four remaining endpoints have active outputs and unavailable energies.
All 12 completed receipts have normal termination and converged SCF; raw final
output energies equal the frozen collection values exactly. No native execution
failure is observed and no intervention is indicated.

The second PLM's available Ca compression work is **+16.348627 kcal/mol** in DFT
versus **+16.599907** composite (+16.604726 OMOL); its La differential remains
missing. The original two GFN2 failures on the 1H4I La Glu−0.2 point remain
unavailable. Nothing was rerun or replaced.

The frozen collector and analyzer produced:

- `workspaces/accommodation_torsion_20260920/prepared_v3/dft/partial_collection_1203771_v3.json`
- `workspaces/accommodation_torsion_20260920/partial_DFT_result_v3.json`
- [PARTIAL_EXECUTION_v3.json](PARTIAL_EXECUTION_v3.json), including all receipts,
  current statuses and input/output pins.

Read-only checks confirmed 12 completed endpoints, 9 available displacement
comparisons, exact raw-energy extraction and explicit missing second-PLM La
values. **Zero new scientific calls.** Final allocation cost remains pending
while job1203771 continues independently. This snapshot replaces no earlier file.
