# Nikasha

Nikasha is the current name of the La/Ca coordination-chemistry discriminator,
previously called Alquemia. The integrated PLM manuscript is the current delivery
target. Naming does not change a scientific protocol, calibration or prediction.

The repository currently uses standalone scripts, without an installable Python
package. `scripts/nikasha` is a thin entrypoint to the existing
`scripts/affordable_workflow.py`: it keeps the selected interpreter, arguments,
working directory and environment. Existing commands continue to work.
Protocol IDs, historical Alquemia schema names, cache keys, data paths and the
repository remote retain their recorded identities.

## Use the current workflow

Run from the repository root with the existing CPU environment:

```bash
export NIKASHA_PY=/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

"$NIKASHA_PY" scripts/nikasha standard --help
"$NIKASHA_PY" scripts/nikasha baseline --help
"$NIKASHA_PY" scripts/nikasha standard ensemble --help
```

`standard` supports explicit source requests, preparation, validation, execution
and collection. Its released automatic route uses the fast MACE+GFN2-solvent PQQ
scorer only for compatible PQQ preparations. `--mode dft-reference` preserves the
native r2SCAN-3c/CPCM reference route. `baseline` retains the staged contextual
water-preparation/DFT workflow for supported prepared sites. The optional
three-fold `standard ensemble` descriptor remains developmental.

Read-only validation of the actual completed release plan:

```bash
"$NIKASHA_PY" scripts/nikasha standard validate \
  --plan workspaces/pqq_fast_release_20260920/standard_1H4I_v1/plan.json
```

This launches no molecular calculation. For fresh source requests and existing
allocation wrappers, follow the [release commands](../diagnostics/pqq_fast_release_20260920/COMMANDS.md),
substituting `scripts/nikasha` for `scripts/affordable_workflow.py` if desired.
Keep every plan/output path explicit and use a new workspace for a new run.
`nikasha --help` intentionally preserves the older workflow's general-development
help; use `standard --help` or `baseline --help` for the current scoring routes.

You may invoke the executable directly with the intended Python environment on
PATH, or add this repository's `scripts/` directory to your own PATH. No global
installation, environment upgrade or package publication accompanies this name.

## Interpretation and current research

Scores describe protocol-specific coordination chemistry. Functional association,
affinity and physiological occupancy remain distinct evidence. Unknown PLM
structures are predictions; source-domain and preparation failures stay visible.
The completed shared-pool and adaptive geometry tests did not improve the
reference decisions enough to earn promotion. Read the
[September22 delivery](../diagnostics/nikasha_recovery_20260922/DELIVERY.md) for
the recovered comparison with DFT, actual research outcomes and manuscript
figures/tables. The released default remains unchanged.

The [agent operating guide](AGENT_PIPELINE.md) links current releases and dated
research records. The [PLM source inventory](../diagnostics/nikasha_plm_sources_20260922/SOURCES.md)
locates existing protein, genome, tree, motif and transcript tables. It introduces
no new biological matrix, substrate assignment or cohort rescore.
