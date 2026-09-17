# Protein-H correction completed; global model disagreement persists

All four full-protein calls completed, jobs **1200681/1200682**. The uniform
preparation corrects 4,467 stretched protein H bonds using the existing topology
and ff19SB equilibrium lengths. Every heavy atom and PQQ atom is unchanged;
charge, protonation, assembly, explicit waters and atom ordering are unchanged.
No new DFT. Baseline/default and original scientific records remain untouched.

| Measured quantity | Medium, original | Medium, H-corrected | Large, original | Large, H-corrected |
|---|---:|---:|---:|---:|
| Raw R = E(Ca)−E(La), kcal/mol | −405333.457127 | −405338.274741 | −404662.146779 | −404414.357895 |
| Direct contrast gradient L2 norm, eV/Angstrom | 10.914341 | 9.654376 | 100.247352 | 152.938320 |
| Maximum atom contrast-gradient norm, eV/Angstrom | 4.746752 | 4.040485 | 32.368257 | 59.596153 |
| Sum(abs(Ca−La atomic charge)), e | 5.115218 | 4.834678 | 11.940021 | 12.878487 |

The medium raw contrast changes by −4.817614 kcal/mol; large changes by
+247.788884. Medium/large disagreement increases from 671.310348 to
**923.916846 kcal/mol**. The corrected geometry does **not** resolve global
response. The H rule was chosen independently of energies/labels and is retained
as a defined preparation version; a favorable score is not its validity test.

This does not isolate solvent, learned response or heavy-coordinate strain as
the cause. Both calculations remain vacuum evaluations. No calibrated S or
classification exists. The changed core-H geometry has no matching DFT
endpoints, so **hybrid results remain unavailable**, explicitly rather than
reusing incompatible old DFT energies. No mechanical correction is enabled.

## Verified implementation and costs

New protocols:

- `mace_polar_1m_analytic_vacuum_protein_H_v1`
- `mace_polar_1l_analytic_vacuum_protein_H_v1`

Four real-fixture geometry tests verify identity/chemical-state preservation,
fixed heavy/cofactor positions, every H's unique attachment and target length,
orientation, rigid-transform covariance, rejection of corrupted real bond
inventories, and unavailable results. The complete **46-test regression suite
passes in 58.579 s**. No dummy scientific output was substituted.

Primary pair evaluations: **117.323761 s medium**, **243.858296 s large**.
Peak GPU tensor allocations: 10,251,213,312 / 16,087,777,280 bytes; reserved
12,140,412,928 / 22,540,189,696 bytes. Peak worker host RSS 1,727,168 /
1,990,624 KiB. One A5000/16 CPUs/64,474 MiB host request per job.

Successful allocations: 140 + 267 = **407 GPU-allocation seconds**. Including
the failed preflight below: **410 GPU-allocation seconds, 6,560 allocated
core-seconds, 426.892 actual CPU seconds**. No inference failed or retried.
Queue waiting is excluded. Preparation and analysis timing are retained where
measured; an exact all-inclusive preparation CPU total was not captured.

Preparation v1's preflight (job 1200679) demanded bitwise equality when replaying
floating-point geometry on another CPU and stopped before any inference. Its
dependent job 1200680 was cancelled by the owner without allocation. v2 uses
the already tested 1e-12 Angstrom replay tolerance; **the stored XYZ bytes are
identical in v1 and v2**. The original failure, manifests and receipts remain.
Hashes still enforce exact stored input bytes; fixed heavy atoms require exact
coordinate equality. No scientific tolerance was relaxed to pass a score.

## Next research step

The charge tracing excluded a near-singular restoration denominator and located
the broad response in field-dependent updates. The H experiment removes a
verified input concern but leaves the main response problem. Continue the
active goal by testing a coherent affordable solvent contribution on the saved
densities, followed by matched benchmark development. This report is an
intermediate result, not completion of the discriminator goal.

## Artifacts and commands

[Plan](PLAN.md); [unrounded summary and pins](result.json).
`workspaces/mace_hydrogen_20260916/pilot_v2/preparation.json` records every
source-to-prepared H move. The `medium/` and `large/` subdirectories contain
finite manifests, copied implementations, inputs and collections named
`collection_job_1200681.json` and `collection_job_1200682.json`.
`output_audit_v1/` exports all 9,141 mapped atoms, charges and direct gradients.
`sacct_final.tsv` includes the failed preflight and cancelled dependency.

Inspect the completed medium collection without another model evaluation:

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python \
  workspaces/mace_hydrogen_20260916/pilot_v2/medium/implementation/mace_hybrid.py \
  collect --manifest workspaces/mace_hydrogen_20260916/pilot_v2/medium/manifest.json
```

The same runner supplies dry-run, execute, collect and report. Reproduce
preparation into a fresh directory (submits no calculation):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python scripts/mace_hydrogen.py \
  --medium-collection workspaces/mace_analytic_20260916/pilot_v1/collection_job_1200470.json \
  --large-collection workspaces/mace_large_20260916/pilot_v2/collection_job_1200525.json \
  --bond-audit workspaces/mace_response_trace_20260916/preparation_audit_v1/result.json \
  --agreement diagnostics/mace_hydrogen_20260916/PLAN.md \
  --output workspaces/mace_hydrogen_20260916/preparation_reproduction_v1
```

No H pilot remains running.
