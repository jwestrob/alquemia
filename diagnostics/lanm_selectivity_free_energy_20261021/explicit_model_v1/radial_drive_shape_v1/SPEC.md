# Independent radial drive shape — source and design only

Root declaration `85bad5ca` fixes the scientific question: can one independent radial value/slope basis per metal relieve paired P/B energy versus dimer-force conflict in the existing family? This is not a new collective or stiffness term. No build, molecular evaluation, solver, fit, reference calculation, or periodic integration is released by this bundle.

## Functional and units

Use exactly the qualified vacuum source, damping and units in INPUTS.json. Let R_i = r_Oi − r_M, r=|R_i|, r0=2.50 Å, psi=(r²−r0²)/r0², and retain the original f(r), phi(r)=(1−exp(−kappa r³))/r³ with phi(0)=kappa. Add

    delta h_i = 3 e_M f(r) phi(r) psi(r) R_i.

The coefficient e_M is dimensionless. The new keyword is `POLAR-RADIAL-DRIVE <500|501> 36 e_M`; its two coefficients are appended after the original twelve parameters as `La_radial_drive`, `Dy_radial_drive`. Their prospective development box is [−1,1] each. This is an explicit finite exploration range, not a reference accuracy tolerance or global stability/force qualification. The original twelve bounds and u/v→c/d mapping are unchanged. Separate records make e distinct from vdW epsilon.

Both native direct RHS receive the SAME delta h. A=A0+H, H, T, intrinsic alpha, POLPAIR, overlap f, eta_i, eta, original a/b/c/d/R/epsilon, energy units and solver/spectral thresholds are unchanged. The new term is pair additive, contains no environment sums, and contributes only on water O. It vanishes for a bare ion, no metal, no water, and at infinite separation. The stationary dual energy remains −K/2 u·(Ep0+h_old+delta h); the existing common hooks, full native local-frame torques and updated u/v must be used. There is no extra independent scalar energy or H correction.

## Conservative derivative

With B=3 f phi psi and B'=3[(f' phi+f phi')psi+2r f phi/r0²],

    d(delta h)/dR = e [B I + (B'/r) R R^T],       r>0.

At contact, delta h=0 and its Jacobian is exactly −3 e f(0) kappa I. The implementation uses the original regularized phi/phi' branch and avoids division at r=0. Add this Jacobian to the existing O/O h-derivative block and its negative to the O/metal block. No derivative with respect to another water is introduced. The new explicit stationary gradient is −K/2 (u_O+v_O)^T d(delta h)/dx. Existing H derivatives, native coordinate/frame derivatives and both induced channels still use the new converged u/v. Equal/opposite O/metal derivatives preserve translation; their coordinate-gradient contraction gives the existing symmetric virial. No force or energy is patched after a solve.

At ideal r=r0, f=1, delta h=0 but the Jacobian is 6 e phi(r0) n n^T. Thus the added basis changes radial force at the ideal center without directly changing its energy/moment at fixed old parameters. Literal stored geometry radii are not clamped to r0: small nonzero values from printed coordinates remain. The static susceptibility is unchanged at fixed geometry and old coefficients because A is unchanged, including away from r0. This limitation is a central interpretation condition for the local 14-column comparison.

## Source boundary and exact zero

Only `polar.f` (two stored coefficients) and `induce.f` (new parser and conditional h/J addition) change in a fresh copied tree. The old local-response geometry calculations, energy, gradient, exporter, moments, mutual operator, dense spectrum routine and unsupported-route guards are byte preserved. An absent keyword or exact e=0 skips all new h/J arithmetic. Original zero/absent early return is retained when all original and new coefficients are zero; an e-only case activates the same guarded finite-vacuum path. Nonfinite, duplicate, out-of-box or unsupported pair records fail closed. The extra token check reads a character token so a trailing word cannot be ignored. Native readback is emitted only for nonzero e, leaving zero/absent numerical readback unchanged.

This source-level branch equivalence plus representative actual zero parity is the proposed basis for a separately reviewed immutable cache/Jacobian bridge. It is not a present cross-binary cache admission. Old artifacts remain pinned to their actual original binary. No source-wide equality, native qualification, cache bridge or physical improvement is claimed from artificial tests alone.
