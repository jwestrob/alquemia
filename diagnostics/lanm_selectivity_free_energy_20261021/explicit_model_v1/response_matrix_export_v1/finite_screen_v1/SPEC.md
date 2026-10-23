# Finite saved-matrix drive/stiffness diagnostic

Question: do independent short-range drive and environment-dependent metal
stiffness changes offer useful coupled energy/moment/susceptibility leverage on
the consumed data? This is neither a fit nor candidate selection. No molecular
program, reference calculation, force evaluation or radial-shape search occurs.

## Frozen family and finite points

Use the original qualified54 A/Ed/Ep matrices, selected17676 parameters and
literal coordinates. For each metal-water O pair, use the previously declared
normalized S2U shape

`f(r)=[(1+t+t²/3)exp(-t)/((1+5+25/3)exp(-5))]²`, `t=2r/Å`.

At each O add `h=a f(r) 3 r_vector/r³` to **both** native driving fields. Other
sites have no added drive. With `s=sum_O f`, use the fixed phenomenological
`g=2s/(1+s)` and add `H_MM=(exp(l*g)-1) I3/alpha_M` to the metal inverse-response
block. All other A entries, native damping lengths, permanent multipoles and
nonpolarization terms remain unchanged. This changes effective metal response;
it is not a claim that native intrinsic alpha is physically correct. Bare ions
have g=h=0; metal-free46 terms are exactly cached. Both shapes vanish for separated
waters, and their geometry dependence would require derivatives in a future engine.

Compute one global `L=log(2)/max_fixture g` from **all54** frozen geometries,
including both metals. Independently for each metal retain exactly nine points
`a∈{-0.1,0,+0.1}`, `l∈{-L,0,+L}`. This bounds every local inverse-response
multiplier between1/2 and2. Do not choose values using energy/moment outcomes.
There are27 fixtures per metal,486 matrix states including54 baseline states;
reuse the baseline solves at zero. Compose all81 La/Dy point combinations using
their existing state predictions. Do not restrict to matching grid indices.

## Exact energy-conjugate response

For `At=A+H`, `d=Ed+h`, `p=Ep+h`, solve `At*u=d`, `At*v=p`. The stationary
dual functional is `K/2[uᵀAt*v-uᵀp-vᵀd]` with `K=electric/dielec`; these fixtures
have dielec1. Its energy is `U=-K/2*uᵀp`. Candidate total energy is the admitted
baseline total plus `U(candidate)-U(baseline dense)`. This anchoring preserves
the nonpolarization energy; an independent **unanchored** baseline polarization
component replay is mandatory first.

Use reconstructed permanent laboratory dipoles from frozen native multipole/
frame helpers. Candidate total dipole in Debye is permanent plus
`DEBYE*sum_sites(u+v)/2`. No three-decimal total-M estimate is substituted.
For Wz repeating the laboratory z unit vector, the full susceptibility column is
`Hartree/(K*Bohr²)*sum_sites(solve(At,Wz))`, in a0³. Retain exact conversion constants
from the qualified original field19 manifest. Drive has zero susceptibility
derivative in this fixed linear family; stiffness affects both moment and curvature.

The native dual convention is preserved even though these intact-water Ed/Ep
are equal. Changing a driving field or a response tensor alone is not license to
use a nonconjugate reporting convention. Charged absolute moments retain the same
laboratory origin; the literal pair/C/P/B/E3 cycles are charge balanced.

## Admission and numerical policy

Require the actual full54 result and independent review `passed_baseline_matrix_inventory_review`,
with exact runner/config/source/manifest joins. Read both raw export pins for each
fixture and join them to admitted arrays and actual site maps. Original normalized
SPD criterion remains `lambda_min>1e-6*lambda_max>0`. Report spectrum/condition
and solve residual for every point; never clip, symmetrize or use a fallback.
An unstable fixture makes its metal point and dependent paired compositions
unavailable, retained explicitly. Other fixed points are still reported.

Before nonzero points, replay all54 baseline total energies, unanchored polarization
energies and reconstructed full moments with existing CG/rounding bounds; replay
the original ten energy values,34 pair moments, four C/P/B sets and eight E3/mu3
cycles. Join all four baseline susceptibility columns to the independent review's
already-qualified16 old field outputs. Do not rerun broad raw native admission.
Static water moments retain the original six-decimal site bounds and exact frame
reconstruction. Baseline numerical allowances are not new physical tolerances.

## Complete reporting and interpretation

For all81 available paired compositions retain the original10 energy gates
unchanged (1kcal/mol), all34 pair dipoles, four C/P/B vector sets, eight E3 energies
and mu3 vectors, and all28 descriptive field rows. Preserve both field sizes,
estimators, reference bounds and exact derived dependencies; only four state
curvatures are independent. Pair/trimer/full-shell energy and dipole changes
remain together so susceptibility agreement cannot hide a large energy cost.
No aggregate loss, ranking, optimized coefficient, selected point or new physical
response threshold is produced. All15 candidate force checks are **unevaluated**.

The four-amplitude family already has an exact Dy directional obstruction:
`dχzz/dl=-Hartree/(K*Bohr²)*exp(lg)*g/alpha_M*||[At^-1 Wz]_M||²<=0` in both
dimer and trimer, whereas their observed offsets require opposite directions
for exact matching. These provisional-reference offsets are not physical error
bars. The finite grid measures leverage/tradeoffs; failure of its points rejects
only those points, not the full continuous domain or all mechanics. Fixed
nonpolarization terms may contain compensation from the old energy/force fit;
energy deterioration alone does not exclude a future jointly declared energy/
response repair. No such extension is implemented or selected here.

One local numerical thread. Root owns actual saved-screen execution after source
review. No candidate forces, donor/held-out use, GPU/periodic or production release.
