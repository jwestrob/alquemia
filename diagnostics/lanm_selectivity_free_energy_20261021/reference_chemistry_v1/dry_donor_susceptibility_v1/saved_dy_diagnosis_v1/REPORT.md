# Saved Dy field diagnosis: final integration dominates the odd discrepancy

The original ten-cell panel remains **numerically unqualified (14/16)**. Five
saved Dy outputs were examined; no native call, changed threshold or replacement
target was produced. Exact Decimal arithmetic, all raw pins and terminal SCF
excerpts are in `actual_v1/RESULT.json` (`99997dce`) and `RAW_EXCERPTS.json` in
the corresponding workspace. The La control uses the already collected summary
only; no additional La raw calculation or state search was made.

## Direct evidence

All five Dy inputs are identical after removing their EField vector line.
They use charge+2, effective restricted singlet/physical-parent sextet, identical
geometry, basis/ECP, origin0 and grid settings. Reported point counts are exactly
107343 for XC and16804/28624/52765 for COSX grids1/2/3 in every output. D4 is the
same printed −0.005764079 Eh. Property electronic energy matches the native
post-final-integration SCF energy; Final=SCF+D4 closes within the existing token
bound. No field-sign, origin, D4-subtraction or parser energy-stage join defect
was found.

Each SCF stops after37–38 cycles on Energy Check. Energy and orbital-gradient
flags pass; MAX/RMS-density, DIIS and orbital-rotation flags fail. MAX-density
changes are5.36e−6–1.25e−5 versus1e−8; RMS-density3.40e−7–7.24e−7 versus1e−9;
DIIS3.07e−3–3.29e−3 versus1e−8. The final integration then changes the scalar
energy by roughly9.59e−6–9.84e−6 Eh without another SCF iteration. The initial
negative-gap/pre-diagonalization warning occurs smoothly across fields; it is
not evidence of a field-induced state switch or a certified physical state.

For the same five actual jobs, compare their last pre-final-integration SCF
energies to their full-precision final electronic energies:

| Quantity, ea0 | h=.001 | h=.0005 |
|---|---:|---:|
| Final odd-energy slope minus zero dipole | +0.000171507216 | +0.000250652211 |
| Contribution from the final integration change | +0.000127533131 | +0.000230096079 |
| Fraction of observed mismatch |74.36%|91.80%|
| Pre-final odd slope minus zero dipole, diagnostic only | +0.000043974085 | +0.000020556131 |

The nine-decimal printed correction alone bounds its odd contribution to
5e−7/1e−6 ea0. The table uses the higher-precision difference between the last
SCF energy and final property energy; both reproduce that printed correction
within its rounding. The energy curvature contribution from final integration
is −0.012062565 and−0.129107983 a0³. Thus the final stage measurably changes both
derivatives. Subtracting it would put the diagnostic odd residual below1e−4,
but **does not admit a corrected energy surface or rescue the original panel**.

## Truncation and remaining uncertainty

For smooth E(F)=E0−μ0F−χF²/2−βF³/6+…, the leading odd-energy error is βh²/6,
and the even dipole numerator is βh². The latter gives β≈193.16/196.78 au,
predicting odd shifts3.219e−5/8.199e−6 ea0. Its h/half ratio3.926 is near the
leading-order expectation4, whereas the observed odd mismatch ratio is0.684:
the mismatch grows when the field shrinks. Leading smooth cubic response alone
does not explain these data. Higher orders and residual numerical errors are
not uniquely separated by two steps.

This localizes a large **post-SCF final-integration contribution**, not a unique
grid or state error. Identical grid sizes do not make the SCF and final grid
the same, nor eliminate field-dependent quadrature error. Failed density/DIIS
flags remain a possible additional limitation. Energy tolerance is not an
energy-error bound. No residual or stage is discarded.

## Smallest justified prospective numerical-recipe test

Root selected a new **uniform ten-call block**, Dy5 first and matching La5 only
after the actual Dy six-check/native/cost gate. Within every five-cell block,
include the zero and all four original fields. Do not splice one new sign or
zero into the failed original panel, or combine new Dy with old La as one
uniform recipe.

The one consequential change is late SCF COSX grid2→the existing final grid3:

```
%method
 AngularGridX 2,4,4
 IntAccX -1,-1,-1
end
```

Arendt verified these exact controls against the executed Mex consistent-grid
workflow. Its physical grid2/grid3 fingerprints match; that old derivative
test nevertheless failed. It supplies installed syntax/activation evidence,
**not evidence that this recovery will pass**. Do not import Mex's unrelated
`DoEQ false`: the original dry donor input has no method block, and its defaults
must remain unchanged. Keep original grid1, final grid3, XC, geometry, fields,
state, basis/ECP, PBE0-D4, SCF settings, origins, both estimators and all16
numerical checks. Use4MPI/one node/default grants and grant-derived MaxCore.

Admission must verify actual grid2=grid3 fingerprints (radial/angular/pruning,
weights, cutoffs, fitting and point count), unchanged grid1/grid3/XC, all
property/energy/D4 joins and raw stopping flags. Actual sixteen checks determine
qualification; there is no promised repair, subtraction, tighter step, relaxed
criterion or automatic additional retry. This diagnosis authorizes no native
execution by its author; root owns declaration/release/submission.
