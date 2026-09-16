# MACE memory-work execution checkpoint

**In progress,2026-09-16.** Current job1200381 on RTX A5000,16 CPUs,64474MiB
requested host RAM. Current campaign `workspaces/mace_hybrid_20260916/blocked_v6`.
No full-protein result is claimed at this checkpoint. Read live receipts before
relaunching. Baseline/default, model weights and physical inputs are unchanged.

The first pair-blocked implementation passed original-autograd kernel checks
and all four complete GPU core checks (1200365). Further full-system attempts
identified large local neural edge and per-atom product tensors. Those are now
blocked/checkpointed too, including independent field tensor products. Full execution is gated on reproducing all four
original core energies/forces/densities at the predeclared tight tolerances.

Scientific scope and technical details: [agreement](MEMORY_AGREEMENT.md),
[local recovery](MEMORY_LOCAL_RECOVERY.md), [product recovery](MEMORY_PRODUCT_RECOVERY.md),
[implementation](BLOCKED_KERNEL.md). The current scientific protocol is still
`mace_polar_1m_vacuum_r2scan3c_subtractive_pilot_v1`.

## Read-only commands

```bash
cd /groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs
MACE_WORK="$PWD/workspaces/mace_hybrid_20260916"
MACE_PY="$MACE_WORK/software_v1/venv/bin/python"
MACE_MANIFEST="$MACE_WORK/blocked_v6/manifest.json"
squeue -j 1200381
"$MACE_PY" scripts/mace_hybrid.py dry-run --manifest "$MACE_MANIFEST"
"$MACE_PY" scripts/mace_hybrid.py compare-cores --manifest "$MACE_MANIFEST"
"$MACE_PY" scripts/mace_hybrid.py collect --manifest "$MACE_MANIFEST"
```

Job exits collect `blocked_v6/collection_job_1200381.json`; accounting watcher
PID4159484 writes `blocked_v6/job_1200381_accounting.json`. There is currently
no automatic offload retry: an earlier offload failure exceeded the requested
host-memory share and requires implementation inspection. Own pendingH200
1200309 was cancelled without compute after the initial core validation.

Preserved attempts: blocked_v1 jobs1200365,1200368,1200369; blocked_v2 1200370;
blocked_v3 1200371,1200372; blocked_v4 1200379; blocked_v5 1200380. Their failures/results are not
overwritten. 1200372 had already failed when cancellation was requested after
observing excessive live RSS; the correction receipt records this precisely.
