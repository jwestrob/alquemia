# Conditional susceptibility when a second water is added

## Question and interpretation

Measure whether the frozen model reproduces the change in directional electrical susceptibility on adding water 3 to the existing W9-derived M+water-1 fragment, for La and Dy. This is a fixed-geometry response diagnostic, not a fit or a new reference target.

The admitted 18413 direct/mutual ±5% probes leave the W9 paired energy error between −4.3495 and −4.3216 kcal/mol. Every sampled trimer energy and collective dipole remains underestimated; dipole magnitudes move only millidebye against deficits of tenths of a Debye. These are local parameter sensitivities, not proof that every damping function fails. The existing μ3 already measures the first field derivative of E3. The new information here is curvature and its change with environment.

At fixed geometry the linear induced model has susceptibility BᵀA⁻¹B. Changing the direct zero-field driving term leaves that susceptibility unchanged; changing the mutual operator A can change it. This statement assumes the same external-field coupling and supported original baseline mode. No field is enabled in the guarded direct/mutual prototype.

For each metal report χ(trimer), χ(dimer), Δχ=χ(trimer)−χ(dimer), and the paired contrast Δχ(Dy)−Δχ(La). This conditional increment includes the added water's intrinsic and pair response and changes in the metal response. It is **not α3**, a full-shell B, an isolated water polarizability, or a unique separation of polarization, exchange, penetration and charge transfer.

## Frozen geometry and exact finite inventory

Use the already consumed W9 far pair, waters **1 and 3**, and its exact M+water-1 dimer. The two M–O distances are 2.50 Å and O–O is 4.683748498798 Å. Keep every original parent coordinate, atom order and ghost center. The lab +z direction is the O–M–O bisector; do not rotate or relax a fixture. Both metals retain their original, individually pinned parent basis and auxiliary assets.

| System | Real centers / electrons | Charge / effective multiplicity | Orbital / Aux-J functions |
| --- | --- | --- | --- |
| La M+W1 | 4 / 18 | +3 / 1 | 434 / 990 |
| La M+W1+W3 | 7 / 28 | +3 / 1 | 434 / 990 |
| Dy M+W1 | 4 / 18 | +3 / 1 | 434 / 1008 |
| Dy M+W1+W3 | 7 / 28 | +3 / 1 | 434 / 1008 |

All four inputs retain 28 parent centers. Physical La singlet/Dy sextet labels remain distinct from the frozen-f effective singlets. No equality or cancellation of La/Dy ghost-basis responses is assumed. A common change of origin cancels in each same-charge trimer−dimer dipole; response derivatives are separately origin invariant. Keep explicit lab origins anyway.

Each system receives z fields **−0.001, −0.0005, +0.0005, +0.001 atomic units**: exactly **16 new reference single points** and four compatible zero-field SP reuses. Set h=0.001 and h/2=0.0005 au. Their physical scales are approximately 0.0514 and 0.0257 V/Å. The larger field is 0.744% of the bare +3 ion field at 2.5 Å; that is a scale comparison, not a proof of linear response. The larger increment improves curvature signal relative to a 0.0001-au illustrative manual example. Actual qualification below remains mandatory.

The original-baseline model inventory is **19 programs**: sixteen matched field moment/energy calls; one zero-field Dy-trimer identity check; and ±h bare-La-ion sign/scale/factor checks. All three model controls precede the sixteen field calls. Keep the exact selected 17676 parameters and original native source/binary. No prototype edits, model gradients, Cartesian FD, coefficient search or fitting.

## Reference recipe and origins

Reuse the qualified common-parent-basis/trimer SP recipes exactly: electronic PBE0, original def2-TZVP orbital and def2/J auxiliary assets, RIJCOSX, DefGrid3, VeryTightSCF, NoAutostart, Guess PModel, MaxIter 300, original restricted state and real-metal ECP46/55. Retain all full-parent GTO/Aux-J centers, original real/ghost masks and DelECP conventions. These four compatible SP zeros already have **D4 off and no EnGrad**; neither is added. Keep original failed printed SCF criteria, all active stopping signals and final-integration stages visible.

Only field/property directives are added: `%scf EField 0,0,F` and `EFieldOrigin 0,0,0`; `%elprop Dipole true` and `Origin 0,0,0`. Both input origins are the metal's lab coordinate origin. ORCA field units are atomic units; explicit property-origin coordinates use geometry units, while printed `Position of electric origin` is Bohr. Use full-precision property dipoles and exact effective-ECP nuclear charges. Preserve old COM-origin properties and their independently reconstructed common-origin values as zero-field evidence.

The native ORCA convention and finite-difference cautions are documented in the [ORCA 6.1.1 finite-field manual](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/finEfield.html) and [electrical-property manual](https://orca-manual.mpi-muelheim.mpg.de/contents/spectroscopyproperties/electric.html). Neither the manual example nor VeryTightSCF is an accuracy guarantee for these cells.

## Estimators and prospective numerical gates

For each system, with E in Hartree and μz in ea0 at the lab origin:

```
χ_E(s) = −[E(+s)+E(−s)−2E(0)]/s²
χ_μ(s) = [μz(+s)−μz(−s)]/(2s)
μ_E(s) = −[E(+s)−E(−s)]/(2s),  s ∈ {h,h/2}
```

All four systems must first pass finite/native/property admission. Then require:

- For each system, |χ_E(s)−χ_μ(s)|≤**0.5 a0³** at both steps, and |χ_E(h)−χ_E(h/2)|≤0.5 and |χ_μ(h)−χ_μ(h/2)|≤0.5 a0³.
- Apply the same four comparisons directly to each metal's Δχ and to Δχ(Dy)−Δχ(La), with **1 a0³** limits. Do not infer their consistency from cancellation or per-system passes alone.
- For each system, require |μ_E(s)−μz(0)|≤**1e−4 ea0** at both steps. Report the corresponding conditional and paired odd slopes, retaining their signed differences and propagated print bounds.
- Report full dipole vectors, individual estimates, both increments and the paired contrast at both steps. Do not replace a failed raw-step check with extrapolation, selective averaging, a larger tolerance or new field points.

These are numerical consistency resolutions, **not physical accuracy bounds, confidence intervals or an equivalence margin**. At h/2 the known model energy print curvature floor is about 0.000128 a0³; induced-moment slope floors are about 0.002754 (trimer), 0.001574 (dimer), 0.004328 (increment) and 0.008655 (paired) a0³. Preserve actual bounds rather than substituting these estimates for runtime admission. Quantum SCF TolE does not bound energy error. Preserve separate SCF/final-COSX/property energies and their field dependence.

A failed odd-slope or curvature check can reflect finite-step nonlinearity, inconsistent energy/dipole stages, or SCF/numerical error. Do not assign its cause automatically. Retain the even dipole response and odd-energy step dependence diagnostically. For example μ_E(h/2)−μ(0) contains β(h/2)²/6. Failure stops the qualified linear interpretation and the staged continuation; it does not erase possible nonlinear information or authorize a retry.

For the original model use its native total field energy without adding another −μF correction. Confirm field conversion and the energy-derivative moment. In general dual polarization uses permanent moment plus half the sums of uind and uinp; use the native moment only after establishing the equal direct/polarization scaling for these complete-water fragments and passing energy/moment consistency. The guarded direct/mutual prototype remains unsupported for external fields.

## Outcomes and scientific limits

1. A reproducible model/reference curvature or increment disagreement establishes an additional response defect beyond the already measured zero-field energy/dipole offsets. It supports testing an environment-dependent response model, but does not locate that defect exclusively on water or identify charge transfer.
2. Consistent curvature despite the known energy/dipole deficits weakens a missing **uniform-field incremental susceptibility in this direction** as their explanation. It leaves short-range field-dependent offsets, other spatial field modes, anisotropy and higher orders unresolved. It does not establish that a field-independent three-body energy correction is sufficient: such a correction cannot repair the measured μ3 deficit.
3. Unqualified field-step or energy/dipole consistency yields no susceptibility claim. Preserve the actual curves and causes as unresolved; no automatic solver/grid/step change.

Uniform fields probe global modes, unlike the inhomogeneous metal field. Far O–O separation reduces the relevance of direct water–water overlap but leaves both 2.5-Å metal–water contacts unchanged; it does not remove metal-mediated short-range effects. Same added water/coordinates improve the paired comparison, but parent basis/Aux-J differences and changes in metal response do not cancel. No full α3 is claimed: its exact two-metal seven-term field panel would need 56 SPs. W8 and other near/far results remain saved context, not additional field endpoints.

## Stages, resources and preservation

Prepare source/inputs only first. Root owns release and submissions. The initial scout is the **complete four-field La-dimer quartet**, inside the sixteen; the remaining twelve fields require its actual passing native and numerical review. Native or admission failure stops further release. Preserve all sixteen slots and every original failed attempt; no duplicate calls or automatic recovery.

Each independent reference uses four MPI ranks, default scheduler memory/time, with MaxCore derived from the actual memory grant. The four matching zero-field SPs took 325.973/261.511 s (La trimer/dimer) and 277.025/307.241 s (Dy), each four ranks. Repeating that cost for four fields estimates **5.208 CPU-hours** for the full reference panel; this is an estimate, not a cap or promised timing. Root reviews actual scout cost before the remaining twelve.

`SOURCE_INVENTORY.json` binds all four zero-field RESULT/input/output/property records, exact geometries and prior independent evidence. No current input, output, coefficient, production setting or held-out outcome is changed. No molecular call was made while declaring this design.
