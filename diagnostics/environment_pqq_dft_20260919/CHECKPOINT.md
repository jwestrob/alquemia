# Live checkpoint — native charged PQQ check

Completed compact context results are committed by root as39d9383. Native OMOL
retains25/25PQQ and3/3consumed transfers; its gap contracts79.030608→2.315459model-kcal.
Separate static alpha/GGR native OMOL improves2/6→6/6 and native DFT strengthens
all six margins. This context scorer is research only. Root's water-preparation
promotion is separate (2fa14b0); no scorer/default change here.

## Current scientific work

Job1202478 is RUNNING on64CPUs, four16-rank endpoints concurrently. Exactly eight
native r2SCAN-3c/CPCM points on all four added-charge−1 PQQ contexts; fixed actual
OMOL XYZ, original protonation/water/graph and native SP settings. At the most
recent detailed read,0/8endpoints were complete; first four were SCF13–16, last
four waiting in the finite manifest. Recent iterations took46–54seconds. These
are actual progress readings, not scheduler estimates or converged energies.
Do not infer final classification from unconverged totals.

Scientific manifest:
`workspaces/environment_pqq_dft_20260919/prepared_v2/manifest.json`.
Original partial preparationv1 is preserved. Its input-comment parsing rejection
occurred before any submission/energy and is documented in PREPARATION_NOTE.md.
Two actual-artifact preflight tests passed, zero skips,1.104s.

The submitted sbatch already runs the frozen collector after successful execution.
It will produce `prepared_v2/collection_1202478.json`. Every failed endpoint keeps
its output/receipt, and missing components never become zero or baseline scores.

## Independent finite finalizer

A read-only local observer additionally waits for terminal scheduler state, audits
all eight endpoint receipts and produces:

- `diagnostics/environment_pqq_dft_20260919/COMPLETION.json`, including failures,
  actual allocation cost, and its own small measured observer CPU cost;
- `.../collected_1202478.json` only if all eight endpoints pass, retaining each
  original/context Ca/La energy change, raw contrasts and included CPCM diagnostic
  terms. No threshold/calibration/promotion or new solve.

Actual live process1638834, persistent execution session22209. Exact code/manifest
pins and behavior are in FINALIZER_RECEIPT_v2.json. The first shell-background
launch1635358 did not persist; it made no science call/output. Its attempted receipt
is retained as history; v2 is current. Do not duplicate the observer. Standard
partition is exclusive, so the observer uses local lightweight accounting instead
of reserving an entire compute node for a one-core report.

## Interpretation and next decision

All four added-Asp contexts shift native OMOL upward while all neutral expansions
shift it downward. This couples real anionic chemistry with total-charge/composition
changes; it does not uniquely establish a model defect. Original PQQ core total
charge separated classes; three expanded Ca cases now share total charge with ten
La cases and still separate. The reduced gap exposes a harder comparison as well
as a smaller robustness margin. Do not drop those cases or label-route around them.

Archived neutral-context native DFT ΔR is−2.033/−4.425kcal while its included CPCM
term changes by+32.476/+30.835kcal; native OMOL vacuum shifts−22.850/−35.247model-kcal.
Different densities/methods prevent a unique decomposition. The pending eight
DFT points test whether the charged shifts also occur in the established method.

No GFN2 calculation is authorized for automatic execution by this finalizer or
has been launched. SOLVATION_CANDIDATE.md records a possible compact matched
GFN2/ALPB−vacuum transfer, with primary software evidence and an explicit additive
energy expression. Parent directed waiting for these DFT results before selecting
that experiment. Continue from actual outcomes; no automatic correction rescue.
