# Saved comparison of scaffold model and subsequent response

This is a saved-data calculation, with zero molecular calls, new geometries or
optimizer steps. Root binds and runs it in a one-CPU Slurm allocation only after
the actual complete hybrid parent and baseline response pool exist. A missing
future pool remains `unbound`; it is not replaced by the earlier whole model.

## Fixed question and denominator

Keep the admitted Mex broad common7 pool fixed. Retain both LBT sources, including
the Lu-conditioned provenance of 7CCN. For each source report four ordered stages:

1. The completed Module 2 whole-LBT/Mex diagnostic, unchanged.
2. Baseline hybrid LBT at the **same La and Dy geometries selected by Module 2**.
3. Baseline hybrid selections from the complete 26-geometry parent pool.
4. Baseline hybrid selections after adding all four admitted response terminals.

The 104-cell parent (both regions, both metals, all 26 geometries) must replay.
Only baseline residues 3–14 enter the primary conditional comparison. Expanded
residues 2–15 remain the parent's separately recorded fixed-geometry sensitivity;
no absolute energies are subtracted between partitions and no favorable region
or source is selected. The new pool must preserve all 26 inherited geometries,
both metals and all four terminal slots, with 26–30 geometries / 52–60 cells.
No rejected search trial is added to its candidate menu.

Define the conditional balanced exchange
`R = (E_Mex,Dy - E_Mex,La) - (E_LBT,Dy - E_LBT,La)`.
Each system retains its own complete chemical inventory through this one-target
exchange. Mex retains three Nd spectators and its represented waters; LBT has
256 physical atoms and no represented waters. This is not an affinity, population
or experimental preference prediction. Physical MACE La1/Dy6 and effective
closed-shell frozen-f GFN2 remain explicit and unchanged.

The three consecutive changes separate: Hamiltonian change at the same selected
coordinates; reselection within the already saved parent pool; and additional
hybrid response plus its paired cross scores. These are an ordered accounting
identity, not a unique causal energy decomposition. The fixed Mex contribution
cancels from each successive change. Preserve per-component values and signed
`R(7CCN)-R(7CCO)` plus its absolute magnitude at every stage. The existing 1 kcal/mol
interpretation resolution is retained as an `abs(R) <= 1` flag, not uncertainty.

## Admission and arithmetic

Reuse the completed Module 2 comparator's `run()` and native-difference helper;
replay its search/pool/conditional outputs. Reconstruct the fixed Mex common7
from the frozen pool's saved-only admission/summary functions, including every
raw paired point and both terminal guards. Never invoke a pool execution path.

Reuse the response wrapper's zero-call `prepare()` for all prerequisite, full
104-cell, source, state, derivative and component-gradient admission. Reuse
`saved_search()` for all four response results/checkpoints. Restore every final
paired cell through the qualified hybrid helper, replay exact coordinate aliases,
both source-local menus and all selected native component works. All four full
force residuals, nonconvergence flags, termination reasons and bounds are retained.

Native differences use the existing eV and Hartree conversion constants and
component order. Compare recomputed saved scalars with the existing comparator's
`1e-7 kcal/mol` arithmetic tolerance; this is not a new physical acceptance gate.
Native point/state/coordinate/configuration and full-gradient admissions retain
their existing exact or previously qualified gates. Missing/failed evidence makes
the final comparison unavailable. A real missing-cell test removes one admitted
La cell from a copied final menu and requires unavailable status.

Record exact parent/source/config pins and application costs without summing
cumulative search counters with allocation counters. Scheduler receipts stay
root-owned. No backend, numerical retry, optimization or source selection is
permitted. The adapter has no molecular evaluation call path.

## Interface

Bind `CONFIG.template.json` with the actual hybrid and final response pool pins.
Run inside Slurm with the qualified Python environment:

```text
python COMPARE.py --config CONFIG.json --output NEW_OUTPUT_DIRECTORY
```

Outputs are `RESULT.json` and a compact `REPORT.md`. Before any final result,
source inspection and AST parsing are the only local checks.
