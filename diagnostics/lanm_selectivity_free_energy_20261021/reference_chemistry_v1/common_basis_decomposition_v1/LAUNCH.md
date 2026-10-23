# Native launch contract

No native calls were made during preparation. Root binds the independent review and writes the actual release. The two scouts are part of the68-cell inventory; the remaining66 release requires their raw-admitted results and `passed_scout_review` receipt.

```bash
python /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/common_basis_decomposition_v1/REFERENCE.py run --panel /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/common_basis_decomposition_v1/prepared_v1/PANEL.json --release ROOT_RELEASE.json --cell CELL --out /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/common_basis_decomposition_v1/run_v1/tasks/CELL/execution
```

Execute through RUN.sbatch: four MPI ranks, scheduler-default memory/time, allocation-derived MaxCore. `stage` is `scouts` or `remaining`. The remaining release has `scout_results` as the exact two-ID-to-result-pin mapping; its scout review repeats that mapping, panel and independent-review pins.

Native energy-check convergence is retained without converting false printed density or DIIS criteria into new convergence gates. Each output additionally requires both native Aux-J counts to match its full parent basis. Physical state and reference accuracy remain unqualified.
