# Local-response implementation readiness

22 October 2026. Bounded source/capability inventory; zero builds, native calls,
MD, parameter changes or new physical validation. Paths below are repository-relative.

## Answer

**The fitted nonzero local-response correction cannot currently run qualified
periodic explicit-water MD on either CPU or Tinker9 GPU.** The baseline native
engines have periodic/dynamics machinery, but the correction is an isolated
small-system CPU implementation. Its absent/zero bypass preserves baseline
execution; using that bypass would omit the fitted correction.

The CPU source explicitly rejects periodic boundaries, replicas, Ewald,
neighbor lists, prediction and systems over 28 active polarizable sites
(`source/induce.f:7989–7996`). Its hooks currently affect vacuum `dfield0a` and
`ufield0a`. The 57-program/636-FD admission qualifies the declared static
energy/force/field panel, not an integrator or periodic system. The separate
qualified donor prototype extends selected chemistry; it does not supply a
periodic water/protein sampling implementation.

## What already passed

| Evidence | Actual scope |
| --- | --- |
| GPU POLPAIR result `0e0b3acd…` | Four upstream Ewald/non-Ewald cases, 182 assertions, plus four published La/Dy water/acetate CPU comparisons. Maximum energy discrepancy 3.91e-9 kcal/mol and force-component discrepancy 4.981e-5 kcal/mol/Å. These are original POLPAIR/VDW parameters, not the new h/H correction. |
| GPU dynamics result `eb666c6d…` | 127 assertions: Verlet/Monte/Bussi trajectory regressions (40 each), VRESPA initialization only (7). Upstream fixtures do not establish the selected metal/water/protein Hamiltonian. |
| New CPU local response `82886b76…` / actual review `7be6786d…` | Small vacuum MUTUAL route, dual RHS, energy, Cartesian/rigid forces and static electrical conjugacy. No periodic or GPU correction execution. |

Tinker9 commit 44bcd7c898e68f6827c112df79b8361e443730d2 was built in CUDA
double precision for arch90. Those build/regression successes remain reusable
baseline evidence; they do not imply the new keyword or Hamiltonian exists on GPU.

## Concrete engineering gaps

1. **Scalable CPU dispatch.** Apply h to both complete native direct fields and
   −H*x to both complete mutual-field actions, including Ewald/list routes.
   Preserve original alpha, damping, POLPAIR, reciprocal/self terms and native
   frame torques. Extend energy/force/virial dispatch once per route. Replace
   the dense coordinate Jacobian with the existing map's O(number of selected
   metal–water pairs) contraction and selected-site storage; retain zero bypass.
2. **Define the periodic Hamiltonian.** The frozen radial shape has a nonzero
   tail. Minimum-image hard truncation changes it and can create discontinuities.
   The existing design proposes lattice-image sums with a declared derivative
   tail bound, consistent cell/force derivatives, molecule imaging and boundary
   convention. No actual implementation or admitted tail bound exists here.
3. **Preserve stability guarantees.** The current per-solve dense spectrum is
   unsuitable for a solvent box. The bounded local inverse-response factor does
   not prove total operator positivity. A matrix-free strategy must preserve the
   existing normalized spectral criterion; ordinary PCG convergence or estimated
   extremal eigenvalues alone cannot replace its certificate.
4. **GPU port after CPU reference.** The pinned Tinker9 common field wrappers,
   two-channel PCG, PME and energy/force infrastructure are available. New
   coefficient/site data, geometry reductions, dual field actions, local
   gradient/virial kernels and the fused `emplar` route are missing. CUDA and
   OpenACC are separate implementations/qualifications. HIPPO EXCHPOL is a
   different Hamiltonian and is not a substitute.

## Smallest bounded next engineering step

Reuse `J/sampling_extension_design_v1/IMPLEMENTATION_MAP.md`: prepare an
isolated CPU common-hook/O(Q) refactor with unchanged vacuum physics and exact
zero bypass, then qualify it against the existing admitted 57-program panel.
Before enabling a periodic flag, freeze an image/tail convention and a scalable
stability certificate. The map already proposes subsequent periodic static,
reimaging/strain and short NVE gates; their future inputs/budgets must be bound
before execution. No MD stage is released by this note or by component fitting.

Unknowns remain periodic-tail/stability implementation cost, actual CPU/GPU
throughput and numerical dynamics behavior with nonzero fitted terms. Physical
transfer, held-reference acceptance and converged exchange sampling remain
separate unmet stages in PLAN.

## Exact inspected authorities

Here `J` in prose expands to `explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1`.

- `diagnostics/lanm_selectivity_free_energy_20261021/PLAN.md` — SHA256 `73698c7b4aedc0d6fe0cdfa929b5cc63abca6b33b63962cf0a04872037a52470`
- `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/sampling_extension_design_v1/IMPLEMENTATION_MAP.md` — SHA256 `1e74d50ca8a8f29dfc065f99377803980d2ead18687ed298535b8c06fcc65722`
- `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/sampling_extension_design_v1/SOURCE_EVIDENCE.json` — SHA256 `015e958f10e5cbc10f90803c360cea468139ae7add53c4bc1c9cd269919a32c0`
- `diagnostics/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/native_v1/IMPLEMENTATION.md` — SHA256 `c6bb53b649c0c772ad47ef45c0934c46c465db3df2f878508ff98e0e8bc350e5`
- `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/native_v1/source_v1/MANIFEST.json` — SHA256 `82886b766d192b1fb0b89012dcd2bebe6315a2dcc83979b982441da0d464a3b4`
- `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/response_matrix_export_v1/joint_repair_design_v1/qualification_v1/independent_actual_review_v1/REVIEW.json` — SHA256 `7be6786d9388753735885f32b4b7164de2397f1cf5c6bc9e5d128b12a0e4f16f`
- `workspaces/lanm_selectivity_free_energy_20261021/gpu_engine_qualification_v1/execution_v2/RESULT.json` — SHA256 `0e0b3acde13b537543aa87ef5a342c71fe3b0207dc955908e88250b779bf497a`
- `workspaces/lanm_selectivity_free_energy_20261021/gpu_dynamics_regression_v1/execution_v1/RESULT.json` — SHA256 `eb666c6d50107da7415e8dce25b7acda811773012717099c4f2e7bdae591609d`

All 18 CPU/Tinker9 source pins listed in SOURCE_EVIDENCE were verified against current bytes. Existing report and literal CPU guards were read; no native result was regenerated.
