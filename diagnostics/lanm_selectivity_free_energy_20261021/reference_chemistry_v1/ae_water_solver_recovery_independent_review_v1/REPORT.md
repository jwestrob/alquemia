# Independent AE Dy–water recovery source review

**Passed source and actual zero-call preparation; no unresolved blocker.** This review authorizes no submission itself.

- Only nested QuadRegionStart changes1e-4 default to1e-6; exact original geometry and scientific deck retained.
- Fresh native directory, NoAutostart and PModel prevent using old unconverged GBW. FourMPI and actual-grant75% MaxCore remain.
- Native original collect_one checks normal ending, SCF,73electrons/sextet/noECP, atom identities/coordinates, finite EnGrad and energy closure. Wrapper adds164AO/charge/UHF/version/MaxIter and exact executed input.
- TRAH macro norm is checked against printed active TRAH threshold2e-6 separately from unchanged TolG2e-6; micro residual never substitutes. Ordinary SCF retains ConvCheckMode/signals/last gradient.
- Raw native RESULT and outer recovery RESULT are separate and identify wrapper/input delta; absent spin/occupation evidence cannot qualify state.
- One endpoint only; native or admission failure stays failed with no retry. Success would fill one missing provisional reference, not prove determinant/global state, pure f-electron effect, force accuracy or exchange selectivity.

Final source `8cf62605`, configuration `fd17852b`, manifest `4ae5132c`, preparation `6bd089ff`; the15 saved author controls pass. Independent metadata/pin checks: 54. Zero native calls; no broad tests repeated.

The added review correction explicitly reads the native TRAH gradient threshold separately from TolG and retains ordinary-SCF convergence mode. Original17260 failure and input/output hashes remain intact.

[Exact review receipt](/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/ae_water_solver_recovery_independent_review_v1/REVIEW.json).
