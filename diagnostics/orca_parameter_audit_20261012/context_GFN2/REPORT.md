# Solvent factorial complete: 40/40 components qualified

All 80 declared Ca/La vacuum/water inputs finished in both independent backends.
Every one of the40 solvent-transfer components passes the separately frozen v2
gate. These components are ready for root's matched MACE/context/accommodation
comparison; this report does not itself claim improved biological discrimination.

The primary Linux xTB6.7.1 results reproduce all16 available exact archived
corrected standalone energies. Maximum independent solvent-component error is
**1.24519e-7Eh**, below the declared1e-6Eh tolerance; maximum all-atom charge error
is9.19400e-8e, below1e-5e. There are no failed, inconverged or missing inputs.

## What changed and what did not

The20 prepared context/geometry rows cover q9z4j7 and c5axv8, original/expanded
pockets, and a common origin plus four archived donor proposals. Both metals
receive all five geometries within each context. Exact source-atom/cap mapping,
paired coordinates, charge/spin and preparation gates are recorded in
`../context/REPORT.md` and the pinned PREPARATION.json. Geometry transfer changes
no water/proton/PQQ inventory within a context; context-only atoms stay fixed.
The original versus expanded search had also changed one selected donor mode,
which is why cross-evaluating the actual proposals matters.

All80 absolute xTB–tblite energy comparisons fail v1's1e-6Eh gate. This evidence
is retained. The causal s9 diagnostic identifies a difference in the
non-self-consistent three-body dispersion term. That geometry-dependent term
cancels in the same-geometry ALPB-minus-vacuum component actually used by
Nikasha. Root explicitly approved a new component qualification protocol after
the scout; this is an inspected development revision, not a claim that the
original test passed. [SCOUT_REPORT.md](SCOUT_REPORT.md) retains that chronology;
[COMPONENT_V2_PLAN.md](COMPONENT_V2_PLAN.md) freezes the derivation and gates.

Primary xTB parameters, both solvers, geometries, finite candidate pools and
class bands were unchanged. No empirical energy offset was fitted or subtracted.
All qualified components are differences of actual full primary endpoints.
tblite0.7.0 only checks the independent solvent component and atomic populations.
No production default, historical result or calibration was changed here.

## Outputs and exact costs

- [COMPONENTS.tsv](COMPONENTS.tsv):40 metal/context/geometry rows, primary
  vacuum/water energies in Hartree, accepted solvent components, independent
  errors and input hashes.
- [INPUTS.tsv](INPUTS.tsv):all80 inputs, raw energies, charge residuals, original
  absolute gate failures, reference-availability flags and result paths.
- [RESULT.json](RESULT.json):counts, software/result pins and actual per-process
  receipts. Complete endpoint charges, SCF output and runtime dependency hashes
  remain in the owned workspaces.
- Machine-readable integration contract:
  `workspaces/orca_parameter_audit_20261012/context_gfn2_component_v2/FINAL_COLLECTION.json`.

The factorial used160 singlepoints (80 primary,80 independent), including the
eight scout calls; the causal diagnostic adds four calls. Total measured process
CPU is2490.296s, summed worker wall2495.408s, peak single-process RSS111104KiB.
The factorial alone used2446.437 measured CPU-s. These are actual process
measurements on two concurrent one-thread workers under Jacob's login-compute
authorization, not Slurm allocations or measured GPU use. Preparation/model
loading by other agents is outside these totals. No ORCA, DFT, MACE, optimization,
new geometry or Slurm job was executed by this agent in the solvent task.

The80 primary calls used486.154 CPU-s; the80 independent checks used1960.283
CPU-s. The measured elapsed span for the76-input post-scout batch was1195.658s.
Independent qualification overhead is reported separately from primary inference;
these measurements are not a new whole-pipeline throughput benchmark.

Four actual-scout parser/algebra tests pass. They verify preservation of the
failed original gate, exact same-coordinate subtraction, missing-state handling,
and the causal dispersion diagnostic. These are distinct from the164 scientific
process calls above. No missing scientific executable was simulated.

## Next operation

The bounded execution is finished; do not rerun it. Root's existing
`scripts/compare_context_factorial.py` consumes the completed solvent collection
and MACE results. Its selection rule and conclusions belong in the combined
report; no new optimization or broader calibration is implied by completion.

```bash
python diagnostics/orca_parameter_audit_20261012/context_GFN2/TEST_REAL_SCOUT.py
```

Completion event `context-factorial-solvent-v2` was durably queued to root,
receipt `context_gfn2_component_v2/completion_event.json`; queue command returned0.
No owned molecular worker remains active. This branch recommends using these
qualified components for the finite comparison, while retaining the absolute
energy mismatch as an explicit limitation of cross-backend total-energy reuse.
