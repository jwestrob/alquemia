# Two-call La-dimer response-matrix scout

Qualify an opt-in reporting path on the exact consumed W9 La water-1 dimer.
Root binds a reviewed isolated build and submits one CPU/task, one thread.
Two native `analyze system.xyz EM` calls maximum: export OFF first, then ON.
No reference, field, force, fitted coefficient or altered Hamiltonian is run.
The executable and build/source receipts remain unbound during preparation.

The OFF input is byte-identical to the admitted zero-field dimer. The ON key
adds only `RESPONSE-MATRIX-EXPORT`. Reuse the admitted old zero result and four
La-dimer field outputs from the original 19-call susceptibility panel.
OFF must pass baseline parity before ON starts. Any native/parser/numerical
failure stops the scout; both slots remain represented and no retry occurs.

## Prospective numerical checks

- Native moment/energy/induction parsers remain frozen. The prior baseline
  comparison permits the sum of both print half-units: energy/components
  1e-8 kcal/mol, site dipoles 1e-6 Debye, total printed M 1e-3 Debye. Exact
  equality flags are also retained. Despite `DIGITS 10`, analyze prints
  energies/components to eight decimals, sites to six, M to three; induction
  summaries print ten. ON/OFF outputs must be exactly equal for these values.
- OFF has no export files/announcement. ON has exactly two announced numbered
  files, one per native EM induction solve. Parse the complete V1 schema,
  finite values, exact site/atom/type/name/coordinate/alpha maps and native
  constants. Both files must report exact restoration; external field is off.
- `A=diag(1/alpha)-T` is in Å^-3, Ed/Ep in e/Å², u/v in eÅ. Check both channel
  actions, direct-RHS replay, symmetry, dense-versus-native action, and all
  exported diagnostics. Use a floating arithmetic allowance
  `128*dimension*machine_epsilon*max(1,relevant magnitudes)`; this is an
  arithmetic comparison allowance, not a physical accuracy limit. Export
  fields have 17 decimal digits. No averaging of an asymmetric A is allowed.
- Require the assembled normalized matrix `sqrt(alpha) A sqrt(alpha)` to have
  positive eigenvalues and lambda_min >1e-6 lambda_max, as in the declared
  response design. Report its spectrum/condition and those of A. This only
  qualifies numerical stability on this fixture.
- The native CG summary precedes its final peek update. For its printed RMS q,
  the pre-peek residual norm is bounded by `sqrt(N)*(q+0.5e-10)/DEBYE` in the
  numerical convention used by induce.f. The post-peek bound multiplies this
  by `||I-peek*A*diag(alpha)||2`; retain a roundoff allowance. Resolve PEEK
  from the pinned key/parameter files (native default1). Require actual Au-Ed
  and Av-Ep norms within that bound; they need not be exactly zero. With
  lambda_min(A)>0, propagate residual/lambda_min to bound the difference from
  the direct dense solve, including its measured residual and roundoff.
- Reproduce native polarization energy from `-K/2*u.T*Ep`, K=electric/dielec,
  using the actual native solution from export1 and the actual component's
  half-last-digit. Compare the dense-solution energy using its propagated
  solution-error bound. Map export2's native u to the six-decimal site file.
- Predict all four *existing* field-site vectors using
  `solve(A,Ed+Bz*h*Hartree/(K*Bohr))`. Compare to the old raw site vectors with
  their half-last-digit, post-peek residual/conditioning, and arithmetic
  bounds. Predict the centered moment-slope vector at both steps as
  `sum_sites(solve(A,Bz))*Hartree/(K*Bohr²)`; compare to independently summed
  old site differences, propagating the same two field-solution bounds and
  site rounding. No additional field calculation or a posteriori tolerance
  change is allowed.

The two exports must agree on matrix/RHS/maps within arithmetic allowance;
record exact equality observations. Qualification includes raw source/input
and rebuilt executable/build-receipt pins, and actual native/parent costs.
No synthetic suite or broad native-admission rerun is requested.

## Interpretation

A pass admits a fixed-geometry baseline operator for finite energy/dipole/chi
screening. It does not qualify candidate physics, coordinate gradients, force
fixtures, the original25 development gates, any other fixture, or GPU/periodic
operation. A matrix snapshot supplies no dA/dR, dEd/dR or dEp/dR. Printing and
CG bounds here are numerical bounds under the frozen linear operator; they
are not quantum/reference error bars or physical tolerances.
