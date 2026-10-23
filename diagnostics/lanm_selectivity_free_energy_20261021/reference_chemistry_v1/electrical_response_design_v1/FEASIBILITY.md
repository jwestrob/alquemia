# Fixed-geometry electrical response: conditional feasibility

## Recommendation and status

A native ORCA finite-field test is feasible on the four consumed W8/W9 La/Dy clusters. The smallest useful **two-step directional** panel has **16 new single points**: four parents × two field signs × two magnitudes, with four existing zero-field energies. Recommend the existing laboratory z direction for a first directional test; retain every Cartesian coordinate. This measures alpha_zz, not an isotropic polarizability. A three-axis extension would require 48 field points and is outside the minimal proposal.

**Source assessment only:** no executable adapter, molecular input, field magnitude, numerical acceptance tolerance, or release has been created. Actual release remains conditional on the pending common-parent-basis review and a prospective field/accuracy declaration. This document does not authorize automatic continuation after a failed numerical check.

## Native convention

ORCA 6.1.1 uses `%scf EField Fx,Fy,Fz` in atomic units; 1 au is approximately 51.4 V/Angstrom. A positive z component is a physical field pointing along +z. Set `EFieldOrigin 0.0,0.0,0.0` explicitly, at the metal position, in geometry input units. The default field origin is already zero, but the default dipole-reporting origin is the center of mass. [Official finite-field documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/finEfield.html).

Request `%elprop Dipole true` with `Origin 0.0,0.0,0.0` for the finite-field dipoles. The manual provides numerical dipole differentiation and an illustrative 0.0001 au increment, while warning that second energy derivatives require careful increments and tight convergence. That example does **not** qualify a step for these clusters. Analytic CP-SCF polarizability is supported, but is a separate calculation, not included here. [Official electrical-property documentation](https://orca-manual.mpi-muelheim.mpg.de/contents/spectroscopyproperties/electric.html).

## Charged-cluster algebra

For one fixed cluster and field F along the chosen direction, use Hartree energies and atomic-unit fields:

```
E(F) = E(0) - mu*F - alpha*F^2/2 + O(F^3)
alpha_E(h) = -[E(+h) + E(-h) - 2*E(0)] / h^2
alpha_mu(h) = [mu(+h) - mu(-h)] / (2*h)
mu_E(h) = -[E(+h) - E(-h)] / (2*h)
Delta alpha = alpha(Dy,W9) - alpha(La,W9)
            - alpha(Dy,W8) + alpha(La,W8)
```

Evaluate alpha_E and alpha_mu at both h and h/2. All four parent comparisons and the paired Delta comparison need their own step checks; individual error cancellation must not be assumed. Dipoles are available from the same field single points, so this crosscheck costs no extra molecular calls. Report full dipole vectors, but a single applied direction supplies only one tensor column.

For origin O shifted to O+a, the scalar potential changes by F dot a. At fixed charge Q, E changes by Q F dot a and mu by -Q a. The symmetric second derivative and dipole derivative therefore remain invariant. Here Q=+3 throughout; the coordinates and origin stay fixed. No geometry relaxation, translation, rotation, solvent response, or nuclear polarizability enters this test.

## Exact reusable authorities and matching

[EVIDENCE.json](EVIDENCE.json) pins the four original RESULT/input/XYZ/output records and the full parent orbital/auxiliary assets through common-basis PANEL `2f6d7f86`. Native executable SHA is `38b5f057452fef275c0a1b98d270ad03d0411dab70453820d1c678bd732f6c83`.

| Parent | Original job | Electrons | Orbital / Aux-J functions | E(0), Hartree |
| --- | ---: | ---: | ---: | ---: |
| W8 La | 17306 | 88 | 391 / 919 | -641.920061602464 |
| W8 Dy | 17307 | 88 | 391 / 937 | -647.397264879816 |
| W9 La | 17308 | 98 | 434 / 990 | -718.291956177106 |
| W9 Dy | 17309 | 98 | 434 / 1008 | -723.770010766971 |

Preserve PBE0, D4, def2-TZVP O/H, exact metal NewGTO/NewAuxJGTO/NewECP, def2/J, RIJCOSX, DefGrid3, VeryTightSCF, NoAutostart, Guess PModel and MaxIter 300. La/Dy cores are 46/55; charge +3 and effective restricted multiplicity 1 remain unchanged, with physical La singlet/Dy sextet labels retained. Each parent keeps its own complete basis; W8 and W9 do not have identical basis sizes. Full real clusters need no new ghost centers. Do not copy the manual example's ReadMOs behavior.

Prospective single points remove only EnGrad, add the declared field/origin and dipole reporting, and use four MPI ranks with actual default-grant-derived MaxCore. No solver, grid, occupation, spin, basis or convergence change is implicit. Native input echo, field/origin, AO/Aux counts, charge/electrons/ECP, restricted state, termination, active stopping signal and all residuals must be retained for each point.

## Numerical qualification is still open

The original zero-field outputs are EnGrad calculations, with post-final-COSX energies. Reusing E(0) is a defensible same-Hamiltonian proposal, **not an observed zero-field SP replay**. Preserve the finite-point SCF/final-integration/property energy stages and their closure; compare the energy curvature with the dipole slope. A stage discrepancy or disagreement leaves the response unqualified. No automatic zero-field repeat or settings repair is authorized here.

Actual raw zero-field records expose limits hidden by the old RESULT's empty `printed_scf_criteria` arrays: W8 Dy stopped on its gradient criterion; the other three stopped on Energy Check. All four failed the printed density/DIIS criteria. W9 Dy's orbital gradient was 5.1099e-6 against 2e-6. Native success and the VeryTightSCF label do not supply an energy-error bound.

Before release, declare h, h/2 and a scientifically useful absolute response resolution. Require finite complete outputs, two-step agreement for every parent and the paired contrast, and energy/dipole response agreement at that resolution. Preserve odd energy slopes, nonlinear response, residuals and stage changes instead of averaging away inconsistencies. Richardson extrapolation may be reported diagnostically; it cannot replace failed raw-step checks. With per-energy error bounded by epsilon, curvature error is at most 4*epsilon/h^2 per parent and 16*epsilon/h^2 for the four-parent contrast. SCF TolE is not that bound. Thus neither the manual's illustrative h nor printed energy precision alone closes qualification. Failed qualification ends this finite test; tightening or another step panel would be a new declared protocol.

Keep electronic and D4 energies separate. At fixed geometry the original geometry/charge-based D4 contribution is expected to be constant under this field perturbation; verify the native values rather than assuming it. Preserve final-COSX corrections and their field dependence. The old zero-field paired correction is not a bound on its field curvature. Old dipoles use metal-dependent center-of-mass origins; do not compare them directly with new zero-origin odd slopes or infer response from their raw difference.

## Finite cost and interpretation

The four historical EnGrad calls used 1.886 runtime-times-allocated-rank CPU-hours. Their startup/SCF/property modules excluding gradients imply 1.220 CPU-hours. Scaling those two historical proxies gives **about 4.88–7.54 CPU-hours for 16 SPs**, or 14.64–22.63 for 48. These are estimates, not bounds: field-dependent SCF convergence may cost more. Use four MPI ranks/default resources per independent call, preserve actual costs and failures, and count any initial scout inside the finite panel.

This tests fixed-geometry response of the specified finite-basis frozen-core ECP/PBE0 Hamiltonian. It does not establish a physical Dy state, complete-basis accuracy, or the missing frozen-core response. Full-cluster alpha and its W8/W9 metal contrast do not by themselves decompose pair versus cooperative polarizability, prove a unique charge-transfer mechanism, or validate/refit a metal polarizability parameter. A direct model comparison would require a separately declared matching external-field observable; no such model calls are included here.

Evidence extraction used one local thread and zero molecular calls. Existing results and inputs remain unchanged.
