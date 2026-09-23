# Nikasha: approved parallel search and basin pilots

## Agreement and scope

Jacob asked what additional useful experiments could run in parallel. Root
proposed (1) two structure-informed alternative starts, (2) archived solvent
ranking followed where justified by finite solvent-guided proposals, and (3) a
small direct local basin-width feasibility experiment. Jacob replied: **“I
approve.”** This is approval to implement and execute those contained pilots.
No new per-command approval is required. The earlier scaffold/chemical-state
agents' designs remain separate; the redox DFT proposal is not part of this
execution plan.

All comparisons use real already-consumed reference inputs, native float64 OMOL,
and the existing native GFN2 ALPB-minus-vacuum recipe. Preserve source chemical
states, water/proton inventory, paired source mapping and all production defaults.
No new DFT, folds, MD, flexible label fitting or full-population rescore is in
this pilot. A changed protocol's old-band transitions are explicitly development
transfer results; no new threshold is fitted to the pilot or unknown PLM cases.

## Common eight-source pilot

Use the completed adaptive pools for these exact sources, retaining original
and adaptive scores side by side:

1. 1H4I, consumed Ca-associated crystal.
2. 4MAE, consumed La-associated crystal.
3. q9z4j7-pqq-la_model, canonical Ca-associated reference.
4. q88jh5-pqq-la_model, canonical Ca-associated reference.
5. a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-1.
6. a0a3f2yly8-pqq-la_model__conditioned_Ca__seed-1_sample-3.
7. a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4.
8. a0acd6b9f2-pqq-la_model__conditioned_La__seed-1_sample-4.

The latter four are known La-associated development structures spanning a
remaining error, a repaired inconclusive and conditioning-dependent behavior.
They are selected to diagnose mechanisms, not estimate independent biological
accuracy. The crystal/canonical controls preserve both reference classes.
Exact source manifests, completed collections and references are in INPUTS.json.

## Ownership and finite work

Second_shell owns structure-informed starts: at most two source-supported starts
per metal/source, 32 searches and 128 GFN2 singlepoints before reuse/deduplication.
Actual transfer uses the target's atoms/connectivity; both metals receive the
same admitted geometry pool. Freeze template choice, source mapping and physical
domain before evaluation. Unsupported starts remain explicit, not random noise
or copied foreign atoms. Exclude target/homology groups for transfer claims.

Khoury owns solvent guidance: first inspect all completed original30/full225
matrix rankings with no new molecular calculations. If discrepancies exceed the
existing 0.1 model-kcal selection scale, define a small finite energy-probe rule
on the common8 panel. Pin point counts and selections before new evaluations;
use exact native solvent endpoints for final selection. Do not restart repeated
composite-gradient optimization or build a surrogate merely to fill a branch.

Water_basins owns basin breadth: begin with one actual physical torsion on
1H4I,4MAE,Q9Z4J7,and A0ACD6B9F2 Ca-sample4, at common 300 K. Freeze coordinate
measure, domain, nested quadrature and boundary/convergence criteria before
evaluations. Retain nonlinear physical energy and state accounting. The result
is a conditional local integral, not binding free energy, whole-pocket entropy
or water occupancy. Failed/domain-dominated results remain unavailable.

Root owns common inputs, cross-branch comparison and execution coordination.
Each agent writes its own versioned scripts/diagnostic/workspace only; shared
executors remain unchanged unless root integrates a necessary focused extension.
No agent may run another branch's chemistry or resubmit completed work. Record
finite task manifests and actual resource receipts, including failures and
cache reuse. There is no arbitrary project-total time/CPU cap. Scheduler requests
must match current allocation policies and real concurrency.

## Outcomes that matter

Report exact Ca/La works, selected candidates, native/solvent contributions,
paired contrasts and old-band decisions against released and adaptive methods.
Keep failures and source coverage visible. Compare both folding-metal contexts
where paired examples exist. A lower endpoint energy, a wider training gap or
moving an unknown protein toward La is insufficient. Useful progress means a
specific known-reference/structural-robustness gain at practical incremental cost,
or a decisive negative result that stops the tested branch.

The pilot comparison can report changed raw contrasts before a full new canonical
calibration exists. Production promotion and a larger scoring campaign are
separate decisions. Agents must save actual results, runnable commands and vault
notes. Root coordinates the shared index before scoped commits.
