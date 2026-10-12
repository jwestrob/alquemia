# Corrected Nikasha parameter/context work — 12 October 2026

## Result and current decision

The live Linux native-ORCA parameter correction is reproduced. Sixty actual
controls are terminal:52 admitted, eight original/default-roundtrip density
failures preserved. Corrected La and14 corrected transition-metal aquo cells
match the archived independent controls. Zn/water population parity remains
unqualified. Read runtime/REPORT.md and the independent runtime_review report.

The strict208 corrected calibration remains unsupported (gap−4.556602kcal/mol).
The complete prior2796-point correction was imported, not rerun; see
recovered_strict208/IMPORT.json. Corrected origin-only203/3/2 and frozen geometry
208/208 remain separately scoped results over correlated folds from25proteins.
No threshold tuning or production promotion. Current operating guide and PLM SOP
now identify the earlier native energetic performance claims as historical.

## Active finite experiment

Two consumed boundary proteins now have exact original/expanded representations
of the same five saved donor geometries. Context expansion had also changed the
fourth selected torsion; this matched factorial separates those effects.
No new optimizer, physical inventory, DFT or folds.

MACE is complete:24 reused endpoint energies plus16 new cross-context calls;
four additional CPU/GPU origin replays agree exactly. Actual20 calls:113.884s
wall,590.696 process CPU-s,4.25GiB peak RSS, zeroGPU. See context_MACE/REPORT.md.

GFN2 component-qualified v2 is RUNNING under Curie, two serial numerical workers.
Actual workspace: workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/.
Read FINAL_COLLECTION.json if present before any rerun. The original four-input
xTB/tblite scout failed absolute-energy parity; a separate four-call diagnostic
isolated nonsc three-body dispersion. That density-independent term cancels
within each backend's same-geometry ALPB−vacuum contribution. Version2 requires
1e-6Eh per-metal transfer parity and1e-5e per-input populations, preserving all
original absolute mismatches and every failed gate. Primary energies stay actual
standalone xTB6.7.1, with tblite0.7.0 as component check.

The frozen v2 matrix is80inputs/backend, including four scout reuses. No rescue
or repeated2796-point panel. See context_GFN2/COMPONENT_V2_PLAN.md.
WorkerPID1406461; durable watcherPID1418506. Exact rearm/receipt command in
workspace WATCHER_REARM.json; event context-factorial-solvent-v2 queues root
01a0a63a-ee36-7483-b726-1f7f5b7f75cc. Inspect receipt before rearming.

## Next exact command

After the solvent final collection exists, run once into this new directory:

```bash
python scripts/compare_context_factorial.py \
  --mace diagnostics/orca_parameter_audit_20261012/context_MACE/RESULT.json \
  --solvent workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json \
  --output diagnostics/orca_parameter_audit_20261012/context_comparison
```

Then independently audit raw components, common-pool selections, origin/adapted
work, context interactions and two-anchor ordering. A positive two-anchor gap
is not full canonical recalibration. Decide any next same-search-policy pilot
from those results; no new optimization is automatically queued.

## Ownership and safeguards

Jacob cleared the pause and explicitly authorizes<=16 TOTAL numerical login
threads. Recorded verbatim atop AGENTS.md. Runtime4/rootMACE8 allocations are
released after completion; Curie owns2 numerical threads. Dalton owns file-only
guard/runtime review; Heisenberg owns completed runtime records and new guarded
dispatch wrapper. Root owns scoring/reporting, current-facing docs and integration.

Parameter guard commitcd97a66 implements exact task-local full-file imports,
shell-aware audits, runtime/profile/state admission and pre/postflight hashes.
Historical original/default controls remain diagnostic-only. The guard's Ca/La
registry uses actual current receipts; it does not globally patch old production
executors or qualify new-metal mixed profiles automatically. The new guarded
dispatch wrapper is separately tested; inspect its report before use.

Two headless Claude attempts produced zero science: first OAuth refresh failed;
one technical retry refused the explicitly authorized login scope. Their receipts
remain under workspaces/claude_parameter_audit_20261012/. Codex completed runtime
work without changing authentication or asking Jacob to repeat authorization.

No Slurm science jobs were submitted. The30minute readiness monitor remains
activePID669874; latest probes still fail account/partition admission. On recovery
inspect actual live jobs; do not duplicate login work or automatically restart
old LanM campaigns. Historical outputs, existing production and concurrent files
are preserved. Scoped commits/pushes remain authorized.
