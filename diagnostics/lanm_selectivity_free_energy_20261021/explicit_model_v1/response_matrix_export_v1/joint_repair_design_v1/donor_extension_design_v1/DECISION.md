# Minimal donor extension: separate local drive and oxygen stiffness

**Decision:** retain the stationary dual-field h/H functional and add a separate,
opt-in carboxylate-O403 channel with **two amplitudes per metal**: one radial
driving-field amplitude and one isotropic oxygen-stiffness amplitude. Four new
amplitudes total is the smallest proposed family with separate drive and
stiffness directions. Independent adjustment of the projected moment and
susceptibility is conditional on their sensitivity rank. This is a source-design recommendation,
not a coefficient choice, fit release or claim that four amplitudes can repair
all observed forces. Wait for the complete Dy susceptibility comparison and the
canonical rotation method check before choosing coefficients or a fitting loss.
Oxygen placement is a prospective modeling choice, not a location identified by
the La result. One tangent susceptibility scalar cannot distinguish donor,
metal or distributed stiffness. Any later selection of this placement requires
discriminating constraints outside its own fitted projection, with root
declaring the comparison; no alternative family or held-out read is launched here.

## Evidence and what the field experiment identifies

The consumed 17621 donor model and PBE0-D4 each pass the original sampled
derivative checks, but their rotation slopes have opposite signs: La −20.390410
versus +28.357940 and Dy −7.316093 versus +26.601243 kcal/mol/Å. Rotation is the
declared q=(1 Å)θ motion and also changes metal–oxygen distances. Model component
splits are bookkeeping, not corresponding quantum energy components.

The independently admitted model40 field result (`71d83185`, review `f21ff157`)
adds a distinct electrical constraint. Along the fixed laboratory tangent,
central model/PBE0 moments are La 0.114375/0.580438 and Dy 0.100684/0.635663 ea0;
the q=.005 rotation derivatives are −1.150972/−.717914 and
−1.173915/−.788829 ea0/Å. A field-independent short-range energy correction cannot
alter these fixed-geometry electrical quantities. It may still be needed for
the remaining mechanical error after the electrical constraints are addressed.

La quantum five-cell susceptibility is now independently admitted (`2b396c26`,
collection `4085e7f4`, all six consistency checks). Its energy curvatures at
h=.001/.0005 are 47.882390/47.881020 a0³; dipole slopes are
47.865830/47.862563 a0³. The model half-step curvature is 54.125717 a0³, about
13% larger, despite its much smaller zero-field projected moment. Increasing
response softness alone is therefore not supported as a joint repair of these
two La observables. This is numerical comparison to provisional PBE0, not a
physical accuracy certification. The full ten-cell panel subsequently admitted
all native endpoints but **failed numerical qualification: 14/16 checks pass**.
Dy odd-energy slopes differ from its zero-field dipole by 0.000171507 at h and
0.000250652 ea0 at half, exceeding the frozen 1e−4 ea0 limit. Its approximate
55.1 a0³ curvature and the apparent paired sign reversal are descriptive failed-
panel outputs only. **Neither Dy nor the paired susceptibility is admitted as a
target or evidence for choosing coefficients.** Saved-output diagnosis must
retain the failure and distinguish numerical/stage/nonlinearity hypotheses.

Two admitted central projections can constrain the need for a change in
equilibrium electrical drive and/or stiffness in that projection. They do not
identify site-resolved polarizabilities, permanent versus induced errors, charge
transfer, Pauli repulsion, anisotropy, the complete response tensor, or the unique
cause of the rotational energy error. A uniform field probes a global mode,
whereas the metal's internal field is spatially varying. Both field steps and
estimators remain numerical controls, not extra independent fitting targets or
physical error bars. Frozen-f ECP, PBE0/basis/state and SCF/final-grid limitations
remain; the canonical rotation experiment tests a different uncertainty.

## One specified functional

Keep all original multipoles, protonation, topology, intrinsic α, pdamp,
POLPAIR, intraligand scaling/exclusions and short-range parameters fixed during
this extension's numerical qualification. Preserve the exact 17621 donor keys;
do not silently replace their La/Dy donor parameters with water-fit parameters.
Treat both chemically equivalent O403 atoms identically, including the farther
oxygen. No nearest-oxygen switch or outcome-dependent atom selection is used.

For each metal–O403 pair let R_i=r_i−r_M, r_i=|R_i|. Reuse the present smooth
radial envelope f(r), its fixed 2 Å⁻¹ range and 2.5 Å normalization, and the
native fixed damping function

```
κ_i = γ_Mi / (pdamp_M pdamp_i)^3
φ_i(r) = [1 − exp(−κ_i r³)] / r³
w_i = 2 f(r_i) / [1 + f(r_i)]
h_i = a_M Q_M f(r_i) φ_i(r_i) R_i          (only O403 RHS sites)
H_ii = β_i I,   β_i = [exp(c_M w_i) − 1] / α_i   (only O403 blocks)
```

Q_M remains the original +3 monopole; the whole dry complex remains +2.
The stable contact series already implemented for φ must be retained.
With no metal or no O403, this channel is identically zero. At a_M=c_M=0,
branch to the exact old path before new arithmetic. There is no new environment,
angular, water, mixed-ligand or frame-dependent fitted coefficient. Bounds and
the nonzero qualification values must be prospectively declared separately.

Let G=A+H, A=α⁻¹−T, K=electric/dielec and f_ext the unchanged native field RHS.
Retain **both** native direct fields and solve

```
G u = Ed + h + f_ext,       G v = Ep + h + f_ext
U = K/2 [uᵀGv − uᵀ(Ep+h+f_ext) − vᵀ(Ed+h+f_ext)]
  = −K/2 (Ed+h+f_ext)ᵀ G⁻¹ (Ep+h+f_ext)  at stationarity.
δU = −K/2 [(δEd+δh)ᵀv + uᵀ(δEp+δh) − uᵀ(δA+δH)v]
```

Donor Ed≠Ep is allowed and must be measured, not removed by changing exclusion
groups. The conjugate moment is the original permanent laboratory moment plus
½Σ(u+v), converted from eÅ. SAVE-UINDUCE contains u alone. With b_t the repeated
fixed tangent vector, χ_t=C b_tᵀG⁻¹b_t, where
C=Hartree_kcal_mol/[K Bohr_Å²]. The drive coefficient a_M does not change χ at
fixed geometry. Because ∂H/∂c_M is positive semidefinite on the donor sites,

```
∂χ_t/∂c_M = −C (G⁻¹b_t)ᵀ (∂H/∂c_M) (G⁻¹b_t) ≤ 0.
```

Thus c_M can stiffen or soften the donor mode separately from the field offset.
The sign and nonzero size of ∂μ_t/∂a_M still require actual response evidence;
positive a_M is not assumed to increase the measured projection. Two amplitudes
have local rank two only if both ∂μ_t/∂a_M and ∂χ_t/∂c_M are nonzero; H also
changes the equilibrium moment and may need a compensating drive adjustment.
Two amplitudes
per metal leave μ, dμ/dq, χ and mechanical slopes overdetermined. Failure would
be informative; it would not automatically authorize more parameters.

## Conservative native integration and precise source hooks

The frozen native tree is manifest `82886b76`; exact files/lines are pinned in
SOURCES.json. Current `localresponsesetup` accepts only 500/501–36 records and
intact 36/37 waters, and rejects donor types 401–404. Current geometry puts β
only on the metal. Merely loosening the atom whitelist would not implement the
proposed oxygen response.

| Hook | Necessary contained change in a future isolated copy |
|---|---|
| `LOCAL_RESPONSE.f:localresponsesetup` | Separate donor keyword/coefficients; exact acetate graph and both O403 atoms; keep native-mode guards, zero bypass and water branch. Permit only consumed eight-site M–acetate, isolated seven-site acetate and ion controls. |
| `localresponsegeometry` | Build the two donor h vectors, two β values and their exact M/O coordinate derivatives. Read actual resolved donor γ/pdamp; do not inherit water γ. |
| `induce.f:dfield0a`, line1253 | Existing shared addition to field and fieldp already has the right dual algebra. |
| `induce.f:ufield0a`, lines1633–1634 | Existing per-site −β_i u_i/−β_i v_i already supports donor oxygen diagonals. |
| `epolar3.f` h correction; `epolar1e` reduced energy | Retain the single −K/2 u·h addition where needed. Do not add H energy again after the stationary reduction. |
| `localresponsegradient`; `epolar1.f:2488` | Replace the hard-coded metal-only contraction with the sum over corrected sites: K/2[Σ_i (u_i·v_i)∂β_i − Σ_i(u_i+v_i)·∂h_i]. Include the equal/opposite metal reactions and virial. |
| `epolar1a:rotpole` and `torque`, lines196/2417 | Recompute original frame-dependent contractions using the **new u and v** at every geometry. Never cache old frame torques or freeze C/O/H frame atoms. |
| `moments.f:199`, full induced readback | Retain ½(u+v), plus permanent moments reconstructed at one fixed physical origin. Extend the exporter guard to this exact topology before using it. |

The full coordinate derivative also includes the original
K/2[uᵀA′v−uᵀEp′−vᵀEd′], evaluated with the corrected solutions. Original native
frame derivatives/torque distribution supply those terms. The new h uses
laboratory pair vectors and H is isotropic, so it adds no separate fitted local
frame rotation. Original non-oxygen atoms can nevertheless acquire changed
forces through their native multipoles and induced response. Holding their
torques fixed would violate energy conservation. An eventual anisotropic
H=R H_local Rᵀ would require R′ terms and is outside this proposal.

Local positive diagonal entries do not prove G positive definite. Retain the
actual normalized spectrum, min>1e−6 max>0, symmetry/eigen residual tests, native
PCG stopping/readback and failure stop; no eigenvalue clipping or repair.

## Required finite qualification before fitting

Use only the two consumed centers, their eight existing rotation endpoints,
and exact previously admitted unaffected controls. Choose one nonzero a and c
per metal prospectively, then test **drive-only, stiffness-only and combined**:
six central parameter states, not a scan. A separately declared execution matrix
must fix amplitudes, Cartesian displacement sizes, call counts and limits.

1. Zero/off must replay saved 17621 energies and available forces on all ten
   donor geometries; water-only, isolated-acetate and both ion controls retain
   exact bypass. Check total charge, permanent local multipoles and all original
   parameter/exclusion bytes unchanged.
2. On the six central states, compare full Cartesian analytic gradients to
   native energy differences, including C/H frame atoms, at two declared small
   steps. Check reaction forces, zero-field global rotation covariance and
   translational/torque closure. Check raw energy from the dual stationary
   expression, both RHS, residuals and the unrepaired spectrum.
3. For those states use the already established ±.001/±.0005 field convention
   and fixed tangent. Verify −dE/dF against the **dual** conjugate moment and
   curvature against the same G inverse. Drive-only χ must retain the baseline;
   c response must have the derived monotonic sign. Retain existing field
   numerical resolutions (1e−4 ea0, .5 a0³ individual/1 paired), printing bounds
   and all failed channels; these are not model–QM accuracy limits.
4. On the eight existing rotation endpoints, verify the combined correction's
   energy/force consistency and mixed derivative ∂F_q/∂F_t=∂μ_t/∂q, with the
   native multipole frames moving. Reuse laboratory origin0 for charge+2; a
   moving COM would contaminate the mixed derivative. Preserve both q steps;
   add no fitted physical threshold to dμ/dq.

No implementation, coefficient evaluation, molecular call, fit, held-out read
or submission was made for this note. The active water fit and all original
source/results remain untouched. Resolution of the failed Dy qualification and
canonical reference review are interpretation prerequisites, not excuses to change the
already frozen reference protocols.
