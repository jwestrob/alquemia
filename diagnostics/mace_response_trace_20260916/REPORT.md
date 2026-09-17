# Charge-update tracing isolates the next preparation experiment

Both instrumented pilots completed: medium job **1200676**, large **1200677**.
Twelve actual MACE calls, zero new DFT. Every replay passes the declared
equivalence gate. Energy changes are exactly zero at saved precision; maximum
force changes are 1.06e-12 (medium) and 4.90e-11 (large) eV/Angstrom. The observer
does not change model tensors, weights, energies or gradients.

## What the traces establish

The normalization denominator is **not nearly singular** on these full systems.
Minimum absolute-sum/sum-absolute weight ratios across endpoints and updates
are 0.929004 for medium and 0.773036 for large. A value near zero would indicate
severe cancellation; these values do not support that proposed mechanism.

The contrast between Ca and La develops during the field-dependent updates:

| Representation / stage | Medium sum(abs(delta atomic charge)), e | Large, e |
|---|---:|---:|
| Full protein, initial local prediction + charge restoration | 1.070347 | 1.166164 |
| Full protein, first field update + restoration | 6.674905 | 6.495219 |
| Full protein, second field update + restoration | 5.115218 | 11.940021 |
| qm33, final | 1.076528 | 1.065031 |
| qm36, final | 1.080098 | 1.067862 |

Net Ca-minus-La charge remains −1 e throughout. Large's final charge-restoration
contribution changes by 26.370181 e in L1 norm between endpoints, partly
cancelling the preceding learned change. This is a coupled calculation; these
terms alone do not identify a unique physical error or prove all global
polarization is undesirable. The initial local prediction gives much less
extensive contrast than the final full-system response.

The observed normalization and two field updates agree with the installed
MACE implementation and the authors' described architecture. The model uses
learned global charge restoration, not an independently solved variational QEq
minimum. [Primary model description, sections II.2.4 and II.2.7](https://arxiv.org/html/2602.19411v1).

## Independent preparation finding

The [read-only bond audit](PREPARATION_AUDIT_PLAN.md) compared the actual
source-mapped protein coordinates with its existing Amber ff19SB bonded
parameters. It used no model evaluation or minimization.

| Protein bond group | Number | Mean length minus ff19SB equilibrium, Angstrom |
|---|---:|---:|
| C-H | 3459 | +0.103309 |
| N-H | 916 | +0.176687 |
| O-H | 92 | +0.224716 |
| C-C | 2390 | −0.000262 |
| C-N | 1469 | +0.001076 |
| C-O | 881 | +0.003391 |

The nearly uniformly long protein H bonds are a concrete preparation concern.
Force-field equilibrium distances are an independent geometrical reference,
not an exact quantum optimum. Some heavy-bond outliers also exist; all are
retained in the output table. Metal and PQQ bond parameters were not invented.

Next: [uniform protein-H bond projection](../mace_hydrogen_20260916/PLAN.md)
preserves all heavy/cofactor atoms and chemical states and tests whether this
input correction changes full-model behavior. No score/label selected the H
targets. This follow-up is part of the active goal, not completion of it.

## Cost, artifacts and commands

Medium: **167 s allocation, 2,672 allocated core-seconds, 184.575 actual CPU s**.
Large: **298 s allocation, 4,768 allocated core-seconds, 317.120 actual CPU s**.
Each used one A5000/16 CPUs/64,474 MiB host request. Total 465 GPU-allocation
seconds and 7,440 allocated core-seconds; no failed inference calls or retries.
Three trace admission/recovery tests pass (6.967 s), alongside six existing
pilot tests (7.139 s). No prediction/calibration/relaxation result is claimed.

Workspaces: `workspaces/mace_response_trace_20260916/{medium_v1,large_v1}/`.
Collections `collection_job_1200676.json` and `collection_job_1200677.json`
pin every energy, force, density, trace and execution receipt. All stage arrays
are saved per task as `charge_trace_arrays.npz` with summary `charge_trace.json`.
`comparison_v1.json` contains exact paired stage quantities. Bond audit:
`preparation_audit_v1/result.json` and `protein_bonds.tsv`.

Reproduce the saved-data comparison without running a model (fresh filename):

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  scripts/mace_response_trace_assess.py \
  --medium workspaces/mace_response_trace_20260916/medium_v1/collection_job_1200676.json \
  --large workspaces/mace_response_trace_20260916/large_v1/collection_job_1200677.json \
  --output workspaces/mace_response_trace_20260916/comparison_reproduction_v1.json
```

The existing `mace_hybrid.py` runner supplies dry-run, execute, collect and
report for trace manifests. Prepare via `mace_response_trace.py` using explicit
source collection, agreement and output arguments. Baseline remains unchanged.
