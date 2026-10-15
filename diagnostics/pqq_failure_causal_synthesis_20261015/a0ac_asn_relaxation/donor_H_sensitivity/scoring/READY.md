# Twelve-point H sensitivity: ready

Executed as5039:12/12 complete. Read ../RESULT_REVIEW.md and the workspace's
ROOT_REVIEW.json. Commands below preserve the original launch; do not resubmit
the completed workspace.

Actual Slurm5029 preparation is complete for both sources and all twelve cells. Independent static review cleared the scorer and actual input schema. This metadata freeze makes no numerical or model calls.

CONFIG SHA: `f5532b200c53ed92b0fd43eda8330cfaffd7f794af344931046ccbf2b7708909`.

The `run` action performs the complete actual preflight before starting the qualified seven-CPU warm worker. Root may submit it directly; `preflight` is optional and makes zero calls. No optimizer, solvent or classifier is invoked.

```bash
sbatch --parsable --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/run-%j.out --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/run-%j.err /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/RUN.sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/CONFIG.json f5532b200c53ed92b0fd43eda8330cfaffd7f794af344931046ccbf2b7708909 run
```

Optional preflight only:

```bash
sbatch --parsable --output=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/preflight-%j.out --error=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/preflight-%j.err /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/RUN.sbatch /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_scoring_v1/CONFIG.json f5532b200c53ed92b0fd43eda8330cfaffd7f794af344931046ccbf2b7708909 preflight
```

The result is `execution/FINAL_COLLECTION.json`; execution/cost and partial-results receipts remain available on failure. All twelve originals and new statuses remain visible. Per-metal hydrogen corrections and origin-to-fixed-flip/relaxed-La Ca-minus-La works retain unrounded energies and the predeclared 1 kcal/mol interpretive resolution. Root owns submission.
