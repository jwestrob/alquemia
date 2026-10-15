# Four-point native GPU reuse scout: ready

Root submitted5066; completion wakePID4052621. Inspect actual scheduler and
workspace SUBMISSION.json/receipts before acting; do not duplicate the scout.

Independent static review cleared the unchanged worker and exact four saved CPU points. Metadata freeze is complete; no numerical imports or scientific execution ran on login.

CONFIG SHA: `3701029eab54a511d7693a5d7567716aae6843fdd54b6d5659a09f45eaf27c02`.

Root submission (combined actual preflight/runtime probe and four fixed GPU calls):

```bash
sbatch --parsable --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_gpu_reuse_v1/slurm-%j.out --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_gpu_reuse_v1/slurm-%j.err /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_gpu_reuse_v1/RUN.sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_gpu_reuse_v1/CONFIG.json 3701029eab54a511d7693a5d7567716aae6843fdd54b6d5659a09f45eaf27c02
```

One GPU, normal scheduler defaults; no CPU fallback, optimization, solvent or production change. Actual results are `execution/RESULT.json` and `execution/EXECUTION.json`, with raw native GPU request/force receipts preserved. Runtime probe and native worker summary are separate evidence. Missing cells remain explicit; all four points must pass the unchanged 0.01 kcal/mol energy and 0.001 eV/Å maximum force-component gates.
