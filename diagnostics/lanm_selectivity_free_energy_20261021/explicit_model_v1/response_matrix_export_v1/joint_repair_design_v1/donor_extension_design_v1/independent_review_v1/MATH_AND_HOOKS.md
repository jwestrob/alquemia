# Independent donor-response math and native hooks

This is a prospective source/design review. It does not qualify donor predictions, choose coefficients, extend the running water fit, or authorize native execution.

## Dual-field answer

Let M=A+H, d=E_d+h, p=E_p+h. The existing stationary functional is

```
L=K/2 [uᵀMv − uᵀp − vᵀd],   Mu=d, Mv=p,
U*=−K/2 uᵀp.
```

It is a saddle in independent u/v; numerical stability concerns the symmetric physical operator M. Native dry-acetate direct and energy exclusions differ, so d and p must remain separate. Sharing h between them does not make u=v.

At fixed geometry, take an external scalar field F along t to add λBF to both right-hand sides, with B repeating t at every active polar site and λ=Hartree/(K·Bohr) in the frozen native units. Then

```
μ_conjugate = μ_permanent + Kλ/(2 Hartree) Bᵀ(u+v),
χ = Kλ²/Hartree BᵀM⁻¹B.
```

Thus a field-independent h changes the equilibrium dipole without changing χ. H changes χ and generally also the dipole. Separate drive and stiffness parameters provide distinct directions; they do not guarantee independent control of arbitrary dipole and susceptibility errors. The relevant h directions must reach the observed moment, and the H directions must overlap the field-response mode. Native u-only moment output is insufficient for general donor exclusions.

For a donor isotropic block H_i=β_i I, with x=M⁻¹B,

```
∂χ/∂β_i = −Kλ²/Hartree ||x_i||².
```

A positive block increase can only lower this projected susceptibility locally. Negative increments can raise it, but only while the actual total M remains positive definite and well conditioned. A bound on individual β_i alone is not a total-operator stability guarantee. Small donor response vectors imply weak stiffness leverage. Do not infer an accessible correction from the metal/donor alpha values alone.

## Conservative geometry derivative

At converged stationarity the additional Cartesian energy derivative is

```
δ(∂U/∂R_j)=K/2 [Σ_i u_iᵀ(∂H_i/∂R_j)v_i
                         − Σ_i (u_i+v_i)·(∂h_i/∂R_j)].
```

Forces are its negative. Native permanent-field and mutual-tensor derivative/torque routines still operate using the changed u and v. Include every coordinate dependency in h/H, all metal reactions, and any environment/coordination derivatives. Scalar isotropic H blocks need no extra frame rotation derivative. An anisotropic block defined in a moving donor frame requires derivatives of that frame; calling it diagonal in its local frame does not remove the torque. A radial laboratory h needs its vector Jacobian, even when its magnitude depends only on distance. Rotation invariance, translation invariance and contact regularity are design requirements, not optional force corrections.

## Minimum changes to the existing isolated CPU source

1. Add an explicit donor selector and chemically complete graph guard; the current enabled guard accepts only metal and intact water types. Preserve absent/zero exact bypass, original parameters, alpha/pdamp/POLPAIR, exclusions and unsupported-route guards. Permit the host donor atoms without assigning them unintended h/H.
2. Generalize geometry setup to construct selected donor h, β (or explicit 3×3 H), and their full Cartesian Jacobians. The current scalar `lrbeta(i)` operator hook already acts at every site; its geometry currently populates only the metal.
3. Keep `dfield0a` adding the same h to both fields and `ufield0a` applying −H to each of u/v. Preserve the actual-operator spectrum check, native PCG stopping and exporter restoration/readbacks. A positive unchanged preconditioner does not replace the total-M check.
4. Reduced energy routes already contract the corrected energy field; do not add an independent H energy on top. Analyze requires its existing −K/2 u·h correction. Compare analyze and gradient-side reduced energies on donor fixtures.
5. Generalize `localresponsegradient`: its present H term is hardcoded to `uind(:,lrm)·uinp(:,lrm)` and one metal β Jacobian. Donor H requires the sum over selected donor blocks; adding selector records alone is incorrect. Retain the original permanent multipole frame torques and virial accounting.
6. Preserve conjugate half-sum moment reporting for enabled response. Uniform-field checks must compare the total energy derivative with this moment, without imposing u=v or confusing native COM output with the shared laboratory origin.

Necessary actual qualification is limited to donor off/zero parity, finite nonzero drive/stiffness/mixed readbacks, actual M symmetry/stability and dual solves, energy-route/field conjugacy, and Cartesian plus ligand-rotation derivatives that move local frames and all environment dependencies. These would qualify numerical implementation only. No new test panel is released by this note.

## Identifiability limits

One geometry and one field projection supply a projected moment and susceptibility, not a full response tensor or a decomposition into permanent/induced/charge-transfer physics. Two field sizes and energy/dipole estimators qualify numerical consistency; they are not independent physical targets. A La result cannot identify Dy amplitudes or a transferable donor-type law. Even a successful two-parameter correction at that geometry does not establish uniqueness: coupled metal and ligand blocks can share the same projected response. Existing rotation energies/forces and dipole derivatives constrain additional geometry dependence but do not isolate a unique term. Actual local response directions and a rank/conditioning assessment are needed before claiming any larger family identifiable.

## Final reviewed proposal and admitted evidence

Reviewed Singer's final DECISION.md (`fbd1ff8a`) and SOURCES.json (`cf766589`): one shared radial drive and one isotropic O403 stiffness amplitude per metal, both O403 atoms included with the same rule. No mathematical/source-design blocker remains. In the specified block, ∂β_i/∂c=w_i exp(cw_i)/α_i≥0. The two-observable Jacobian is triangular in drive/stiffness coordinates, and has rank two only if ∂μ/∂a≠0 and ∂χ/∂c≠0. This conditional statement is sufficient; actual response leverage and conditioning are still unmeasured for the donor proposal.

The shape is conservative and well behaved at contact: φ→κ and h∝R as r→0, while f and w are finite and smooth to the derivative order needed here. Both h and H vanish at separation; isotropic H and radial h are rotation covariant, and inclusion of both O403 atoms is permutation symmetric. Native metal–donor gamma/pdamp must be resolved from the original donor parameters. Neither water gamma nor a nearest-O selection is interchangeable with this definition.

Only La susceptibility is presently admitted. The exact La collection/review joins retain all five SPs and six numerical checks: χ_E=.001/.0005 gives47.88238959/47.88101978 a0³; χ_mu gives47.86582969/47.86256267. The admitted model half-step χ_E is54.12571667 a0³, while its projected moment is much smaller than the reference. These two observables support separating driving and stiffness questions; they do not identify oxygen placement or a unique mechanism.

Root reports the full ten-SP panel has ten native admissions but14/16 numerical checks: both Dy odd-energy/zero-moment checks fail. The final design explicitly retains this failure. Dy and paired susceptibility are unqualified for parameter choice or target promotion. This review does not independently re-admit that new full result; Arendt owns its actual audit. The admitted La result remains distinct. No coefficients, fitting loss, alternate field step or threshold are selected here.
