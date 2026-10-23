# Next mechanics discrimination after 17676

Recommend a fixed metal-polarizability response panel first. This is a parameter-response diagnostic at the unchanged terminal point, not another optimizer continuation or a qualified model.

## Evidence and scope

The terminal result `5001cb5798571245d5dcff7c6afbcc8c482e1525f0e8e09e2c4c1b124eb03832` stopped after 22 steps because a complete trial batch gave no actual improvement. Loss is 0.895803724; 9/10 energy and 13/15 force controls pass. Remaining failures are the W9-minus-W8 Dy−La energy error −4.338284 kcal/mol, La W9 rigid RMS 2.042854, and paired expanded-dimer RMS 2.203824 kcal/mol/Å. No parameter is at its bound.

The last measured Jacobian is exactly at the selected point, rank six, condition 149.17. Its unrestricted linear projection reduces loss by only 1.55e−10; the actual final trials are slightly worse. This supports testing another response direction locally. It does not establish global impossibility of the six-parameter family.

Saved component arithmetic for the selected W9-minus-W8 Dy−La contrast is vdW −5.29009248 plus polarization +0.41760198 = −4.87249050 kcal/mol. Bonded terms and permanent multipoles cancel at printed precision. This is a decomposition within this model, not a decomposition of the QM error into uniquely physical causes. W9 also rearranges the eight other waters; it is not a pure ninth-water insertion. Water polarization still depends on the fixed water model.

## Ranked options

### 1. Metal alpha response — implement this finite diagnostic first

Keep all six selected pair parameters and every geometry fixed. Change only the metal `POLARIZE` alpha token, separately for La (1.6597 Å³) and Dy (0.7280 Å³), by −1%, +1%, −5%, +5%. For each of eight probes evaluate the five consumed fixtures of that metal, reusing the five unchanged opposite-metal cells. Exactly 40 new physical cells, at most 80 analyze/analytic-testgrad programs; no Cartesian FD, fitting, new reference or selected replacement.

Use the unchanged 413 weighted residuals and all 25 controls. Reconstruct the existing six-column Jacobian from its actual saved ±0.005 normalized probes. For each alpha step size add two columns in fractional coordinates `(alpha−alpha0)/alpha0`, retain the same absolute cutoff `s6_max × 1e−6`, and report residualized column rank/norms and nested residual projections. Compare 1% and 5% slopes and center-midpoint curvature; isolated probe improvements are not required for a useful coupled direction. Do not treat nonlinear alpha as the exactly linear epsilon parameter. No optimum or new parameter vector is solved.

Native source gives `pdamp_i=alpha_i^(1/6)` and damping exponent `a_ij r³/sqrt(alpha_i alpha_j)`; induced dipoles also depend explicitly on alpha. This experiment measures **coupled susceptibility and damping-length response**. For an infinitesimal metal-alpha change, increasing its O-pair Thole coefficient by one half the fractional alpha change would cancel the O-screening exponent change; that O coefficient already belongs to the six-column span. The remaining alpha direction can still contain metal susceptibility and fixed H-screening response. A useful new span therefore does not by itself identify which mechanism repairs the error, nor validate the intrinsic alpha against experiment.

### 2. Independent metal–water-H screening — conditional follow-up only

Tinker supports `POLPAIR 37 500/501` by atom type. These currently resolve to native minimum Thole values La 0.25 and Dy 0.139, while O uses independently fitted `POLPAIR 36 500/501`. A separately declared small H-pair damping probe, with alpha and all six O-pair parameters fixed, can distinguish missing H-screening flexibility from the coupled alpha response. Use the same five fixtures per metal and residual projection, preferably the same two step sizes if linearity remains uncertain. Do not change water polarizabilities, permanent multipoles or donor terms.

This option ranks second because it is a localized empirical pair correction; alpha has a direct susceptibility interpretation. If both responses are nearly collinear after projection against the original six columns, these consumed fixtures cannot identify them separately. A compensated-alpha experiment (O and H coefficients scaled by `sqrt(alpha/alpha0)`) could then isolate the explicit susceptibility contribution within the same native functional form, but is not included in the first 80-program panel.

## Source basis and limits

Inspected frozen native `kpolar.f` (assignment lines 434–448 and pair resolution 490–521), `damping.f` (original Thole exponent lines 100–109 and direct-field fallback 192–198), `induce.f` (alpha-times-field lines 232–239 and mutual residual lines 387–408), and the pinned `alchemical_path_v1/PATH.py::resolve_pair` (class-based vdW and type-based POLPAIR). The resolver alone omits alpha-dependent damping length; native sources must remain pinned in the executable diagnostic.

Reference remains the same provisional uncorrected ECP-PBE0-D4 surface. No donor result, reserved orientation/radial outcome, explicit-electron result, or bulk sampling claim is introduced. The Cartesian and fresh-radial validation adapters remain unbound because 17676 failed scientific controls. Further physical terms or a parameter fit require a separate finite declaration after this diagnostic.
