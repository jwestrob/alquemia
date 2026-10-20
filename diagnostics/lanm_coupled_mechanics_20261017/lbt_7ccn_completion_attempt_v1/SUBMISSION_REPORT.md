# One bounded model-convergence attempt

Job **15458 is queued for resources** as of 20 October 2026 09:27 UTC.
It resumes the exact 7CCN–La step92 checkpoint with the unchanged whole256
MACE plus fixed-charge ALPB model.

## Why this run

The preceding twelve steps lowered the model energy by4.226978kcal/mol.
The independent PBE0 plus the same solvent comparison gives4.726291kcal/mol
lowering, a0.499313 difference. This supports that particular increment.
The model remains unconverged, and its full-force disagreement with the
reference remains3.579797kcal/mol/Å RMS. The question here is whether the
unchanged model can converge with a substantial, finite additional allowance.

## Exact allowance and admission

- At most80 further accepted steps or200 further distinct evaluations,
  with cumulative ceilings172/292.
- Original force criteria: maximum0.10 and atom-vector RMS0.03kcal/mol/Å.
- Same source, all256 movable atoms, model, solvent, state, physical guards,
  L-BFGS history and trust radius.
- Actual preflight **15457 passed** in13 allocated CPU-seconds with no
  molecular calls or model initialization. It restored113 saved records,
  109 unique cache keys,92 historical requests and the exact checkpoint.
- Independent static and actual-preflight reviews passed. Root decision
  `b6d3f115` and run configuration `8a78864e` bind the actual evidence.

The allocation requests12CPUs and oneH200 under `standard_expand`, on the
previously qualified host. No explicit memory or time override was added.
The completion watcher is recorded in the run workspace as PID1733550.
There is no run result yet.

The attempt stops at the original convergence/failure conditions or its
declared allowance. Budget exhaustion remains incomplete relaxation.
There is no automatic extension, intermediate DFT or other-metal campaign.
Model convergence would not establish force accuracy, stable electronic
state, affinity or metal preference.
