# Native field-energy identity before direct quantum-field input

Declared after the failed charge/dipole sampling follow-up. Active-goal
authorization applies. This is a native solver accounting control, not a new
metal score, another source fit or a relaxation calculation.

## Motivation and source evidence

The exact saved-density fields already passed numerical checks and cost only
minutes for all8endpoints. Compact source fits did not meet the full spatial
screen. Before exposing per-site external fields in a native solver, determine
how those fields enter its actual energy; do not alter induced dipoles and
reuse an energy that omits their interaction.

Pinned Tinker CPU source87050685eff8840d312e2a332cc82c33f63c7c3d:

- `induce.f::induce0c` calls dfield0d to form vacuum d/p and solvent d/p fields,
  adds the existing uniform external field to all4arrays, then stores
  `udirp=alpha*fieldp`, `udirps=alpha*fieldps` before solving induced dipoles.
- `epolar.f::epolar0e` computes native vacuum polarization as
  `-0.5*(electric/dielec)*sum(uind*udirp/alpha)` on allowed sites.
- `esolv.f::egk` calls egk0a (permanent plus induced GK terms) and ediff
  (vacuum-to-solvent induced-dipole correction). ediff uses the permanent
  multipoles and native p-scaled fields; it contains no arbitrary per-site
  external-field input. A keyword alone therefore cannot establish general
  external-field energy accounting.

These observations suggest an identity to verify, not an already qualified
energy expression:

```
I_vac = -0.5*C*sum(mu_vac,d dot F_vac,p)
I_GK  = -0.5*C*sum(mu_solv,d dot F_solv,p)
C = native electric / dielec.
```

F values have native units e/Å² and mu values eÅ; C supplies kcal/mol.
Distinct native d/p fields must remain distinct. Do not replace them with one
field by assumption or call a bilinear AMOEBA expression a general variational
functional without checking the underlying operator.

## Real cases and finite calculation inventory

Use exactly the12completed Tinker framework controls from1201015:3real
protein/water frameworks,4states each. Reuse their full unrounded native energies,
induced dipoles, physical coordinates, parameters, masks and receipts. The metal
remains explicitly outside these framework-only controls. No old biological
result or QM-source approximation is used to choose the identity tolerance.

1. Twelve native direct-field queries, one per original state. Run the same
   initializer, mask restoration, Born calculation and native multipole rotation,
   then call maintained dfield0d. Export all4field arrays and actual native
   constants/parameters. No new induced solve, energy(), force or atom motion in
   this operation. Preserve the native source/library unchanged.
2. Six native no-response energy evaluations: untransformed and recorded rigid
   geometry for each framework. Keep permanent multipoles and all cavity/radius
   parameters fixed. Disable polarization and explicitly initialize induced
   arrays to zero before energy(); this defines the constrained no-response
   state, not a fabricated zero energy. Compute the actual permanent/GK/nonpolar
   energy. Masks cannot alter this state when every induced variable is disabled;
   verify compatibility before sharing it among original masks/tolerances.

Thus6newnativeenergy() calls and12direct-field calculations, noSCF/DFT, MACE,
source fit, charge utility, gradient or trajectory. Reuse the isolated pinned
library and existing64CPU standard/memory allocation/exclusions. Record all
preparation/build/query/energy costs and failures; no project time/CPU budget.
The previous12mutual energies cost17.1summed kernel seconds; this is a short
implementation control, not a new production cost regime.

## Frozen verification and outputs

Compare I_vac with each original native ep. Compare I_GK with
`ep_mutual + es_mutual - es_no_response`; verify em is unchanged for matched
geometry. Require absolute accounting residual≤1e−7kcal/mol. This is an algebraic
identity between evaluated linear-in-dipole terms, independent of whether the
stored dipoles reached an exact infinite-iteration solution; the tolerance is
above ordinary double-precision summation error and below scientific scales.

Require original inventory/settings and finite field/radius values. Verify
source-frozen dipoles remain the actual saved zeros; use no reciprocal alpha
at disabled/zero-alpha sites. Verify no-response arrays are zero and no induced
SCF was performed. Compare native field vectors under the already recorded
proper rotation and translation within1e−8e/Å². Report d/p reciprocity diagnostics
at both old convergence levels, with no new pass threshold for approximate
SCF reciprocity. Keep permanent, polarization, GK/nonpolar and contraction terms
separate so a mismatch is diagnosable.

If the identity fails, investigate the actual native term responsible before
changing any energy. If it passes, it supplies an accounting prerequisite for
a separately defined direct-field adapter. It does not resolve QM/MM covalent
boundary charges, source reaction-field representation, damping, full-cavity
subtraction, partition sensitivity, absolute references or predictive usefulness.
No per-site field injection, native-source patch or score is authorized by this
specific calculation manifest; any such next experiment needs its own recorded
definition under the standing autonomy, not a new user permission request.
