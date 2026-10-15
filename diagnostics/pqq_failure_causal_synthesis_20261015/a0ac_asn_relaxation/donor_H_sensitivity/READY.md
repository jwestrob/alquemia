# Asn two-H geometry sensitivity

Preparation only: Ca4 and La4 × origin/fixed flip/actual relaxed-La × Ca/La. The output retains twelve original MACE/XYZ cells beside the derived coordinates. Root owns submission and the later decision about twelve energy calls.

Metadata-only freeze (no numerical imports):

```bash
python diagnostics/pqq_failure_causal_synthesis_20261015/a0ac_asn_relaxation/donor_H_sensitivity/PREPARE.py freeze --output workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_sensitivity_v1
```

Independent static review cleared the final source. Metadata freeze succeeded. CONFIG SHA: `482bcae3a9c830ad08135c9005bba59e78f6299792e88c13635f8c18e640a77a`. Use the immutable workspace wrapper:

```bash
sbatch --parsable --output=workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_sensitivity_v1/slurm-%j.out --error=workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_sensitivity_v1/slurm-%j.err workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_sensitivity_v1/implementation/RUN.sbatch workspaces/pqq_failure_causal_synthesis_20261015/asn_donor_H_sensitivity_v1/CONFIG.json 482bcae3a9c830ad08135c9005bba59e78f6299792e88c13635f8c18e640a77a
```

One CPU, standard partition, no memory/time overrides. Actual graph/parameter/geometry checks execute only inside Slurm. `execution/PREPARATION.json` contains both source rows and all twelve endpoint statuses; `execution/EXECUTION.json` records completion/failure and cost. Parameter ambiguity, missing source mapping or a newly severe contact remains unavailable. No molecular or numerical execution has been performed by the preparation author.
