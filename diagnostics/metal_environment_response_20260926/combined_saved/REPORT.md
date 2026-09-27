# Saved additive comparison: mechanical work added, electronic gate still failed

**The classical components change each metal's hydroxyl response by the same
+0.417937 kcal/mol, leaving the La-minus-Ca response at +0.352795817 kcal/mol.**
Mapped physical directional derivatives can now be compared for the shared
hydroxyl and boundary motions. Adding these components does not repair the
native rigid-transform failure or provide a solvent-consistent relaxation model.
No new molecular evaluations were performed.

Current machine-readable result: `RESULT_v2.json`. The first `RESULT.json` has
identical numerical values; v2 additionally verifies each actual classical
receipt's serialized model identity against the executed manifest. Its predecessor
is retained. Source scientific outputs are unchanged.

## Actual correspondence and force mapping

Twenty-four shared configurations are available: A, B, repeat, rigid transform,
and two-sided hydroxyl/boundary displacements at0.001 and0.0005 Å, for both metals.
The eight classical-only metal displacements have no corresponding quantum cells
and remain explicitly unavailable for combined finite-difference validation.

Both electronic manifests and the classical ledger point to the exact same
archived1H4I source preparation. Every real QM atom and external charge coordinate
was compared against the actual stored classical geometry. Charge lists, electronic
states, native receipt/output energies, classical receipt/model identities and
force artifacts were checked. Maximum coordinate discrepancy is1.42e−14 Å.
No A/B, boundary or rigid-transform geometry was substituted to fill a matrix.

The54-atom electronic region maps onto the9,141-particle physical model. Its three
link-H gradients are projected onto **both** real CB/CA anchors using the physical
chain rule; Jacobians are recomputed at each actual configuration. The archived
sub-microangstrom cap serialization offsets are preserved and rotate with the
system. Real QM atom gradients and external-charge gradients are then accumulated
on the same physical atom identities. Classical force arrays are negated before
adding to electronic gradients. Quantum-core caps are not additional physical
particles. A real-fixture check verifies that this mapping preserves the sum of
the electronic gradients.

The added energy is exactly the declared ledger's retained bonded, MM-MM Coulomb,
MM-MM LJ and real-QM/MM LJ components. No second QM/MM Coulomb energy or induction
term is introduced. This saved comparison does not newly qualify the cap-boundary
approximation, hydration-fit metal12-6 parameters or GAFF2 PQQ cross LJ.

## Work along the original A→B perturbation

| Contribution, kcal/mol | Ca | La | La minus Ca |
|---|---:|---:|---:|
| Embedded electronic | −3.328669 | −2.975873 | +0.352796 |
| Retained bonded | −0.139473 | −0.139473 | 0 |
| MM-MM Coulomb | +0.557409 | +0.557409 | 0 |
| MM-MM LJ | 0 | 0 | 0 |
| QM-MM LJ | 0 | 0 | 0 |
| Combined dry candidate | **−2.910732** | **−2.557936** | **+0.352796** |

These zeros are actual differences of available computed components, not missing
terms filled with zero. This particular hydroxyl-H rotation leaves both LJ works
zero in the declared parameterization. Equal classical work explains why this
experiment's differential response is unchanged; it does not imply all environment
or metal motions would behave this way.

## Physical directional derivatives

All entries below are gradients in kcal/mol/Å along the same unit-normalized
physical source displacement. The residual is finite difference minus analytic
gradient; forces have the opposite sign.

| Metal/motion | Electronic gradient | Classical gradient | Residual,0.001 Å | Residual,0.0005 Å |
|---|---:|---:|---:|---:|
| Ca hydroxyl | −17.226351 | +2.375421 | −0.000516 | −0.001830 |
| La hydroxyl | −15.246885 | +2.375421 | −0.001666 | +0.000809 |
| Ca boundary | −4.209438 | −62.799434 | +0.000071 | +0.000505 |
| La boundary | −4.288811 | −62.802144 | +0.000704 | −0.018082 |

The selected boundary response is dominated by existing classical mechanical
strain, not the small metal-dependent electronic change. This supports inspecting
source mechanics before using the candidate for relaxation. It does not by itself
identify whether the large derivative comes from hydrogen preparation, boundary
partitioning, local geometry or their combination; root's separate component
strain diagnosis addresses that question.

The smaller finite-difference step does not uniformly improve agreement,
particularly for La's boundary direction. These are reported residuals, not a
new loosened admission rule. Original electronic numerical qualification remains
failed; no stronger validation claim is inferred from adding a smooth classical
function.

## Repeats, rigid transformations and missing scope

Repeated calculations retain small combined residuals: Ca5.7e−9 and La−3.6e−7
kcal/mol energy shifts. Rigid transforms retain **+0.0127919 / +0.0113603** kcal/mol
combined energy shifts for Ca/La. Corresponding maximum mapped physical gradient
residuals are0.0607912 /0.1244544 kcal/mol/Å. Classical rigid energy shifts are below
10−9 kcal/mol, so addition leaves the original electronic failure essentially
unchanged. Both original failed gate records remain embedded in the result.

There are no quantum metal-displacement cells and no combined metal-direction
finite-difference result. Solvent is absent. No optimization, motion population,
thermal entropy, binding affinity, classifier or biological test was performed.
The additive dry candidate is **not qualified for full relaxation**.

## Reproduce

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/metal_environment_combine_saved.py \
  --classical workspaces/metal_environment_response_20260926/component_checks_v1/FINAL_COLLECTION_1219497.json \
  --scout workspaces/metal_environment_response_20260926/reference_scout_v1/FINAL_COLLECTION_1219207.json \
  --force-checks workspaces/metal_environment_response_20260926/force_checks_v1/COLLECTION_1219319.json \
  --output workspaces/metal_environment_response_20260926/combined_saved_reproduction.json
```

Three real-artifact tests pass: matched coordinates with rejection of an actual
mismatched B geometry, gradient-sum conservation through the cap chain rule, and
full saved comparison retaining failures/missing cells. No fabricated scientific
values or fresh molecular calls were used. This analysis incurred no Slurm/GPU
allocation; local read/parse/algebra process cost was not instrumented. Original
scientific receipts contain their already recorded costs, not a new execution.
