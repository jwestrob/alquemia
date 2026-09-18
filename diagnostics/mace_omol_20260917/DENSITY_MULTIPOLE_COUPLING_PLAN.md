# Next: exact-density coupling to real protein permanent multipoles

Declared after native supplied-field replay passes. This is the next component
of an explicit polarizable-environment hybrid, not another source fit or a
new biological comparison. The baseline and all earlier failed models remain.
Active-goal autonomy applies; no further per-analysis approval is needed.

## Why this is needed

The native adapter can now solve response from supplied fields and reproduce
its original energy. Actual quantum potential and electric fields already
exist on every environmental atom for the eight normalized vacuum endpoints.
AMOEBA has permanent quadrupoles as well as charges/dipoles; coupling those
requires the spatial second derivative of the potential. Omitting it would
silently change the environment model. This derivative concerns observation
points around a frozen density, **not nuclear DFT gradients or Hessians**.

Pinned `kmpole.f` converts native quadrupoles by `bohr**2/3` before rotating
them. For those actual internal Cartesian tensors, the direct scalar coupling
has the convention

```
V_direct = sum_i [q_i phi_Q(r_i) − mu_i · E_Q(r_i)
                  + Q_i : Hessian(phi_Q)(r_i)]
```

Use the full symmetric tensor contraction (both off-diagonal entries) and
the internal tensor, with no second division by three. Convert native eÅ and
eÅ² into atomic units, then convert the final energy once. Verify this convention
against pinned native source and real exported moments before interpreting
the component. Do not use an unverified XML tensor convention or invent moments.

## Exact inventory and numerical choices

Use the same eight normalized vacuum r2SCAN-3c states: GGR extended/connected
and alpha1F6S/6IP9, each Ca/La. Same real wavefunctions/densities, atom mappings,
source support and all outside-support environmental atoms; no filtering after
inspection. Reuse the already qualified potential/field observations and native
receipts after checking all hashes and coordinate/atom orders.

1. Three native initialization/rotation-only exports of the actual global
   AMOEBA2018 permanent multipoles, one per existing physical framework.
   No energy, induced solve, force or geometry change. Export exact native
   local/global tensors and retain source/axis mapping. The metal remains
   explicitly outside these existing framework exports.
2. Eight `orca_vpot` calls on the saved actual densities. At each environmental
   point, query the center, ±h on each axis and the four ±h/±h offsets for
   each of xy/xz/yz. Use h=0.01 and0.005bohr:37points per environmental atom.
   Obtain diagonal and mixed central second derivatives independently at each
   spacing. No SCF, new population fit, MACE, minimization or FF energy.

Freeze these numerical checks before outputs: finite values; repeated center
potentials agree with qualified saved centers within1e−8au; native global/local
moment reconstruction and trace conventions agree within1e−10 native units;
charge/dipole/quadrupole scalar contractions are rigid-coordinate algebraically
invariant within1e−8kcal/mol; the coarse/fine change in the quadrupole coupling
is≤0.02kcal/mol for every endpoint and≤0.01kcal/mol for every Ca−La pair.
The latter is≤1% of the existing2kcal partition screen and is an accuracy target,
not chosen from a desired classification. Retain max/RMS tensor differences
and all unrounded component energies. Do not enforce a vacuum Laplace trace
at observation points that may overlap electron-density tails.

Pure algebra verification may compare analytic Hessians of the existing real
projected monopoles with the same finite-difference stencil at the real query
locations. Clearly distinguish this numerical formula test from the eight
actually executed quantum-density utility calls; no fabricated quantum output.

The unchanged native multipoles outside source support provide a **component
diagnostic** at this stage. They are not yet a charge-closed hybrid boundary.
Do not claim the diagnostic sum is an environmental correction or classify it.
Report field/monopole comparisons as representation diagnostics, not a new
biological result or a tuned selection of charge schemes.

The previous eight13point utilities took about800summed CPU-seconds; this
37point inventory is roughly three times the spatial workload. Use8one-thread
workers on standard-shared with16GB as before, measuring actual cost and any
technical recovery. No GPU or high-level quantum evaluation is involved. This
is a finite manifest for reproducibility, not a project time/compute budget.

## Target model to finish after this component

A concrete candidate is a mixed representation with **exact density for direct
QM–environment coupling/driving** and the existing frozen projected CHELPG
monopoles **only as a GK reaction-field proxy**. Its tentative energy is

```
E_trial = E_DFT,vac + V_direct,density + Delta G_GK,permanent,proxy
          + Delta I_environment + T_MACE,short(full) − T_MACE,short(core).
```

The permanent GK difference is evaluated in one common full physical cavity
with and without Q-proxy moments, induction disabled. It contains Q reaction
self energy and Q–environment solvent terms; cavity/environment-only terms
cancel. Native vacuum Q–Q and proxy Q–E Coulomb terms are not included. Native
solvation component differences must therefore be isolated, not whole totals.

For induction, retain native environment-only vacuum d/p fields; add the same
exact Q-density field to each. Add the Q-proxy *reaction-field difference* to
both solvent fields. Use the qualified native solver and subtract the compatible
environment-only polarization contraction. All Q-source induced variables are
frozen; no independent classical response on a quantum source. This changes the
QM–MM cross coupling to unscaled density coupling after boundary preparation;
it is an explicit new model, not native AMOEBA cross exclusions by assumption.

This is a coherent candidate expression, **not a qualified full implementation**.
Before its metal pilot, finish and freeze a source/boundary preparation with:

- All actual physical atoms, including the metal, in a common Ca/La cavity;
  use an explicit common source-metal cavity policy, no guessed La FF defaults.
  An identical Ca2018-derived cavity for both endpoints is a candidate requiring
  an explicit declaration and verification of every resulting radius/scale.
- Remove FF permanent moments on every Q-support/cap-anchor atom exactly once.
  Reconstruct the per-residue exterior formal-charge ledger using actual
  AMOEBA charges and fixed bonded recipients. No global neutralization or
  favorable-result-dependent redistribution. Keep exterior dipoles/quadrupoles
  and full physical axis/covalent graphs traceable.
- Document the source cap/anchor treatment, the residual approximate GK source
  density, and missing short-range penetration/response physics. The failed
  compact fits do not become qualified induction sources by renaming them;
  here exact driving fields replace them. GK proxy accuracy remains a separate
  limitation to test through partition/numerical/boundary controls.
- A complete identity limit and finite native execution manifest, convergence
  refinement and unchanged2kcal GGR partition gate. Compare the same four
  consumed alpha/GGR directions without refitting thresholds or claiming blind
  validation. No compatible absolute reference/bands exist for this candidate.

Do not compute a final scalar by dropping an unavailable component. No CPCM
solvation is added to this vacuum-DFT expression. MACE's local readout is not
uniquely non-electrostatic, so overlap remains an explicit hybrid approximation.
The next executed work is the component inventory above; complete and declare
the subsequent metal manifest before running the full candidate.
