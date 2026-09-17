# MACE as a local curvature source: archived physical GGR displacements

Declared under Jacob's active MACE discriminator goal and blanket pilot
authorization, before these MACE outputs. Direct whole-protein scoring and its
short-range component failed the alpha/GGR development ordering. Test a different
role: can MACE supply local stiffness while DFT supplies the energy and gradient?
No fitted biological weights, threshold changes or affinity claim.

## Inputs and finite inventory

Reuse all 20 exact, executed Stage C GGR XYZs and their DFT outputs from
`workspaces/ggr_mechanism_20260915/stage_c_tasks_v1/manifest.json`, verified against
`workspaces/ggr_mechanism_20260915/report_v1/collection_c.json`. Two representations
(58-atom extended amide,111-atom connected peptide), La/Ca, center and the existing
metal displacement ±0.02 Angstrom and peptide rotation ±1 degree. Retain original
hydrogens, caps, charges, waters and all heavy coordinates for compatibility;
do not import the newer global-H geometry. The original physical atom/cap mapping,
axes and 0.05 Angstrom maximum heavy displacement remain authoritative.

Evaluate these 20 structures with each pinned medium/large analytic MACE model:
40 MACE calls, then 40 native CUDA OBC-II calls from their individual saved
monopoles, using the already validated solvent recipe unchanged. Medium primary,
large sensitivity. Zero new DFT, no optimization/Hessian/MD/training. Reuse the
existing finite runner, resource policy and receipts. Prior local batches took
roughly 150 seconds per 20 MACE calls and 27 seconds per 20 solvent calls; this
is a minutes-scale GPU development test, not routine per-site scoring cost.

## Energy and comparison

Candidate inexpensive energy C(q) = E_MACE,vac(q) + G_OBC2(q_MACE(q), X(q)).
At each displaced point, both learned charges and solvent are reevaluated.
The scalar finite differences therefore include charge response. Saved fixed-q
GB forces omit dq/dX and MUST NOT be called gradients of C. Vacuum MACE analytic
forces are compared only with its own vacuum finite energy differences.

For h from the existing physical coordinate, retain odd=[C(+h)-C(-h)]/2,
even=[C(+h)+C(-h)-2C(0)]/2, g_secant=odd/h, K_secant=2*even/h^2.
Report La,Ca and their difference. A peptide rotation follows a curved path;
its directional second derivative is not a Cartesian Hessian eigenvalue.

Compare candidate even terms against the archived r2SCAN-3c/CPCM even terms.
This assesses a DIFFERENT inexpensive Hamiltonian as an approximation to DFT
curvature; it is not an identical-Hamiltonian solver test. Anchored predictions
are Delta U_DFT(±h) = ±h*g_DFT(0) + even_candidate. Compare with actual archived
DFT changes. No rescaling or spring fitting. Also report vacuum and short-range
even/odd terms as decomposition diagnostics, without selecting whichever passes.

## Frozen gates and limitations

- Existing execution charge tolerance 1e-5 e; exact XYZ/input/receipt pins.
- MACE vacuum analytic-gradient check: odd−h*g within
  max(0.02 kcal/mol,5% of |h*g|), matching the prior directional consistency test.
- Curvature-approximation screen: error of each La,Ca,R even term no greater
  than max(0.005 kcal/mol,25% of |DFT even|), in both representations/directions.
  The 0.005 floor is below the 0.01–0.1 kcal/mol local deformation signals, while
  25% is a coarse mechanical-accuracy target, unrelated to biological labels.
- DFT-anchored ±h prediction: error <= max(0.02 kcal/mol,25% of |DFT even|).
  This looser energy bound includes the already measured DFT odd/gradient error.
- Report all signs/negative curvatures. No clipping, soft-mode removal, entropy,
  full-protein claim or numerical relaxation correction. A passing local screen
  still lacks coupled scaffold curvature, stable subspace/trust-region validation
  and independent biological usefulness: response_model_not_validated remains.

If this fails, retain the failure and assess whether it reflects curvature,
first-derivative mismatch, representation or solver behavior. Do not automatically
add DFT calls or widen motions until a new experiment is explicitly documented.
