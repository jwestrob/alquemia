# Frozen MACE monopole / OBC-II solvent descriptor

Approved under the active discriminator goal; declared before solver outputs.
The H-corrected full-protein medium/large raw contrasts differ by 923.916846
kcal/mol. This tests whether adding the missing water reaction energy changes
that discrepancy. It does not change model weights, charges, geometry or labels.

## Energy and scope

For each saved endpoint use

    E_descriptor(M) = E_MACE,vac(M) + G_OBC2(q_M; X, radii, epsilon)
    R_descriptor = R_MACE,vac + G_OBC2(Ca) - G_OBC2(La).

q_M is the endpoint's saved learned monopole coefficient, with verified total
charge. OBC-II contributes only reaction-field charging energy: its zero-solvent
limit at epsilon_out=epsilon_in=1 is zero. MACE's vacuum energy already contains
its direct electrostatics and self terms; do not add point Coulomb energy again.
There is no CPCM term in these MACE endpoints. Nonpolar surface terms are
excluded and exactly cancel from R for identical paired geometry/radii.

This is a frozen-density perturbative descriptor. It omits the dipolar part of
MACE's density in the solvent model and does not feed solvent fields back into
the learned updates. It is neither self-consistent MACE/GB nor an exact solvent
functional for MACE's Gaussian density. Corrected gradients/relaxation remain
unavailable: fixed-charge GB forces alone omit the charge-response chain rule.
No aquo reference, inherited baseline band, calibrated S or class is available.

Use the maintained OpenMM native GBSAOBCForce implementation, with an independent
same-equation check against installed CustomGBForce OBC-II on two real cores.
Pin installed source, version, platform and constants. If a NonbondedForce is
included for API consistency, match its charges, disable its reaction field,
put it in a separate force group, and extract **only the GB force group**.

## Frozen physical parameters

No cutoff/PBC, epsilon_in=1, epsilon_out=78.5, salt=0, water dielectric interpreted
at 298.15 K. Surface area coefficient zero. Retain actual physical radius values
already recorded in the 1H4I state: H 1.2, C 1.7, N 1.55, O 1.52, S 1.8 A;
selected metal 1.8 A for both Ca/La. Thus the geometric cavity is common across
endpoints, with no element-specific sign adjustment. This radius policy is
explicitly distinct from automatic mbondi2 assignment.

Use documented installed OBC descreen factors H .85, C .72, N .79, O .85,
S .96; explicit Ca/La factor .8. The latter is an **unvalidated generic cavity
hypothesis**, not a published Ca/La force-field fit. PQQ uses its actual MACE
charges and these element radii; no guessed force-field PQQ charge is supplied.
Unsupported elements fail. No tuning of radii/dielectric/scales to these results.
Numerical success will not validate these physical approximation choices.

## Finite evaluation inventory and checks

Use both H-corrected full La/Ca pairs from jobs 1200681/1200682. Four primary
CUDA-double GB evaluations; four matched rotated evaluations (existing 37-degree
axis [1,2,3], around the original metal); four identity-limit evaluations
(epsilon_out=1); one medium-La repeat; one medium-La translation [10,-7,3] A;
one medium-La Reference-platform full evaluation. **Fifteen full solver calls**.

Additionally use the original analytic-medium qm33 La/Ca saved real densities
for native-versus-custom OBC-II verification on the Reference platform:
two endpoints times two implementations = **four core solver calls**.
Total **19 GB energy/force calculations, zero MACE calls, zero DFT calls**.
Keep every call/failed attempt and actual cost, including context initialization.

Acceptance before interpretation: charge sums 1e-5 e; identity energy <=1e-6
kcal/mol and force <=1e-6 eV/A; repeat/rigid-transform and CPU/GPU energy/contrast
agreement <=0.01 kcal/mol, fixed-charge force agreement <=0.001 eV/A. Native vs
custom core energy <=0.001 kcal/mol (allowing their documented slightly different
Coulomb constants), force <=0.001 eV/A. No grid or box extent exists in this
nonperiodic no-cutoff model. Fixed coordinates and radii remove solvent-volume
changes between paired endpoints. Report component changes and checkpoint
disagreement descriptively; there is no favorable-score stopping criterion.

Use the existing runner/allocation policy and pinned OpenMM environment, without
installing over anything. Start on one A5000/16 CPUs/64474 MiB. Expected scale is
seconds per solver call plus initialization; cost is unmeasured and must be
reported. No project time/CPU budget. Baseline and old experiments preserved.

## Primary implementation references

[OpenMM GBSAOBCForce API](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.GBSAOBCForce.html)
documents charge/radius/scale inputs and separate solvent energy handling.
Installed `openmm/app/internal/customgbforces.py` supplies OBC-II (Amber igb=5),
its descreen coefficients and explicit zero-salt reaction energy expression.
The installed code and actual matched outputs take precedence over an API name.
