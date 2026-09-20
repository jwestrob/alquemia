# Bounded nonlinear donor accommodation

2026-09-20. Approved under Jacob's active overnight accuracy goal and explicit
permission to pursue promising contained directions. Declared before these new
energies. Baseline, released fast scorer and existing experiments are unchanged.

## Question and why now

Can an inexpensive, solvent-aware potential find physically allowed nearby donor
arrangements for each metal, so a compressed source geometry does not dictate the
contrast? Native DFT has now confirmed the large differential work on the first
PLM extra-Asp path, as well as the small 4MAE control. It supports investigating
this mechanism, not a sub-kcal correction, a relaxed minimum or a biological label.
The solvent contribution worsens differential agreement on that PLM path; retain
its components explicitly. Existing profiles descend to their boundary, so a
quadratic extrapolation is not justified. We will evaluate actual energies.

## Fixed first experiment

Use exactly the four original contexts/mappings from
`workspaces/accommodation_torsion_20260920/prepared_v3/design.json`: 1H4I, 4MAE,
PQQSEQ_83440678cbbd658047c9 and PQQSEQ_07ab500e3df76b30d71c, sample0 for both PLM
cases. All source atoms, caps, charges, spin, PQQ, water and proton inventories
are unchanged. PLM labels are unknown. The existing independent MMOL1770 profile
and the reference-fold campaign are comparisons, not new starting structures for
this pilot.

For each of the eight metal endpoints, start at q=0. Move only anchor-Glu chi3
and, when present in the original role definition, extra-Asp chi2. Rotate whole
groups with the existing source/cap Jacobians; all other modes are zero. Use the
same coordinates and bounds for both metals. Each angle is bounded to [-0.8,0.8]
radian, approximately 46 degrees: a nearby terminal-group adjustment, not a full
rotamer search. Preserve bonds/fixed atoms; reject newly introduced severe clashes
using the existing 0.55 A H / 1.0 A heavy-pair screen. Record donor contacts and
maximum source-heavy displacement. Unsupported trial geometry is a visible
failure, not a baseline substitution or a favorable skipped state.

Energy and analytic gradient:

    E(q) = E_native_OMOL100M,vac(q)
         + E_native_GFN2,ALPBwater,Tight(q) - E_native_GFN2,vac,Tight(q).

Use the verified native gradient recipe in accommodation_response.py: unchanged
native GFN2 parameters, native mixer, 300 K electronic smearing, Tight convergence,
explicit analytic EnGrad. MACE uses the existing pinned native float64 checkpoint
and charge/spin embedding. Project the complete Cartesian gradient through the
physical Jacobian. No numerical DFT gradients, added springs or fitted coefficients.
The primary-to-Tight difference at q0 is retained, not confused with relaxation.

Minimize each endpoint independently with L-BFGS-B on the one/two active angles,
starting once at q0: maxiter=24, maxfun=80, maxls=10, gtol=0.2 kcal/mol/radian,
ftol=1e-12. Supply energies relative to that endpoint's q0 in kcal/mol, avoiding
an arbitrary absolute-energy effect in the stopping test. These finite algorithm
limits define the experiment, not a project compute/time budget. Cache identical
coordinates/method/state within an endpoint and retain all actual evaluations.
No alternative starting points or outcome-selected retries in this version.

Accept a stationary candidate only when the optimizer succeeds, all absolute
projected derivatives are <=0.2 kcal/mol/radian, every angle is at least0.02radian
inside its boundary, geometry checks pass, and energy does not rise above q0 by
more than0.10kcal/mol (the existing solvent-transfer numerical scale). Report an
energy decrease and a boundary-limited result even if stationarity is unavailable;
do not relabel a constrained boundary point as a relaxed minimum.

For eligible interior candidates only, check the small projected Hessian using
analytic composite gradients at +/-0.005 and +/-0.010radian along each active
coordinate. No whole-protein/dense DFT Hessian. Require symmetry and refinement
errors <=max(0.1kcal/mol/radian^2, 0.01*maximum matrix magnitude), and smallest
eigenvalue greater than the observed refinement spectral norm plus0.1. Report
failures without eigenvalue clipping. No entropy/log-determinant or populations.

Report R0, endpoint accommodation energies, R_candidate and their components.
Any old-band decision is explicitly a developmental transfer check, not a newly
calibrated classifier. Do not choose q by a desired class. No automatic promotion.

## Execution and next evidence

Reuse the existing Slurm/ORCA task renderer/receipts and isolated Python engines.
One H200,32 allocated CPUs,200000MiB host memory, as in the working standard
scanner. Keep MACE warm; run the two GFN2 media through the existing runner with
an allocation-compatible MPI layout. Preserve every failed attempt and resume
only exact compatible tasks. Eight bounded optimizer starts; actual low-level
and GPU calls/timing are measured, not assumed cheap. No project wall-time cap.

No new DFT is submitted by this first plan. After results, a separately recorded
validation manifest can use analytic native r2SCAN-3c/CPCM endpoints at the fixed
energy-selected candidates. Canonical PQQ source mappings can be prepared in
parallel without molecular calls. The same successful rule must be evaluated
across established reference classes before considering a scanner change. A
classification movement on unlabeled PLM targets is never accuracy evidence.

Record real-input derivative/mapping checks, unavailable states, all optimizer
traces and exact executable commands. Keep implementation focused on this small
nonlinear response; reuse the existing workflow rather than building another
general molecular engine.
