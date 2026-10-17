# Module 2: saved broad-versus-restricted comparison

## Question and denominator

How much additional accommodation is available under each system's unchanged
energy model when all represented coordinates can move? Retain all six initial
search slots: LBT 7CCO/7CCN × La/Dy and Mex 8FNS × La/Dy. 7CCN is
Lu-conditioned. Mex retains the EF3 target, three Nd spectators and represented
waters. Missing or state-incompatible terminals remain unavailable.

Only saved results are read. This analysis performs zero molecular calls, zero
optimization steps and no reselection outside the completed, source-local paired
pools. Root binds exact result/config/source hashes and immutable snapshots of
actual scheduler receipts. No directory scan chooses the newest or best result.
Versioned failures and the successful original scout remain in the cost ledger.

## Comparisons

1. For each declared search, compare its exact old restricted terminal to its
   new terminal. Recompute component work from saved native units, and full-atom
   force metrics from saved gradients. The initial force is from the original
   start point, never the first trace entry of a resumed allocation. Preserve
   termination, finite budget, physical admission and raw convergence.
2. Compare each new paired pool to its exact old pool: LBT 12767 separately by
   source, and Mex 12731 common5. Keep selected candidates, both per-metal works,
   their component changes, ties and unavailable cells. The pools themselves
   perform admission and selection; this adapter verifies lineage and arithmetic.
3. Retain the conditional exchange from comparison 12769 for continuity:

   `R_new(source) = R_old(source) + (W_Mex,Dy - W_Mex,La) - (W_LBT,Dy - W_LBT,La)`.

   W is additional selected work from the corresponding old per-metal minimum.
   Compute by components and stable per-metal differences, not subtraction of
   large cross-metal absolute totals. Report both source-local R values and the
   signed/absolute source spread. Do not take a minimum across LBT source pools.
   This is an **unmatched scaffold Hamiltonian diagnostic**: whole-peptide MACE
   plus ALPB for LBT versus local MACE plus embedding/MM/ALPB for Mex. It is not a
   qualified common-Hamiltonian preference prediction or binding free energy.
   Module 3 separately tests matched local-electronic/MM representation.
4. Preserve the earlier 12823 electronic substitution result as a separate
   model-sensitivity reference. It used an earlier LBT pool; it does not supply
   substitute energies for this comparison and does not establish GFN accuracy.

## Frozen checks and interpretation

- Exact SHA256 source/config/point pins; exact old selected native components
  and old selected candidates; exact declared source/state denominators.
- Replay native component work using the existing EV=23.06054783061903 and
  Hartree=627.509474 kcal/mol conversion constants. Derived scalar closure may
  differ by at most 1e-7 kcal/mol from different summation order. This is only
  reporting arithmetic, never a relaxed molecular or force-admission gate.
- Raw full-force convergence is unchanged: maximum atom-vector norm ≤0.10 and
  RMS atom-vector norm ≤0.03 kcal/mol/Å. Budget exhaustion remains incomplete
  relaxation even if energies improve or an expected exchange sign appears.
- Recompute final force metrics from the saved terminal gradient and check the
  search report to 1e-10 absolute plus 1e-12 relative. The force is dE/dR;
  gradient/force sign has no effect on these norms.
- A complete six-slot report may contain unavailable scientific cells. No
  sign-based admission, new threshold fitting, affinity or accuracy claim.
- The inherited 1 kcal/mol resolution is descriptive, not uncertainty or a
  classifier threshold. Report whether changes exceed it without treating that
  as a validation gate.

## Costs and saved-data tests

Each actual job/result is counted once. Search counts include the original
LBT La/7CCO scout exactly once; cumulative optimizer counters are not added to
per-allocation native counters. Rejected evaluated trials remain included.
Original failures remain explicit with their measured costs. Missing native
counters remain missing rather than being inferred as zero.

Report molecular/component counters, model/call wall time and observed CPU/RSS
where present separately from scheduler elapsed time × allocated CPUs. GPU
allocation is reported only from measured/recorded allocation data; no utilization
or speedup is inferred. Root should snapshot watcher receipts after terminal
state because delivery metadata can otherwise change while a pin is in use.

The single Slurm analysis also performs real saved-data checks: original selected
native states replay; six-slot preservation; component and force closure;
search/pool lineage; old R=component sum; unavailable-pool suppression using a
removed actual menu; and one-copy scout/job accounting. Test-only masked rows
are labelled and never counted as observations. Output RESULT.json, SEARCHES.tsv,
POOLS.tsv, CONDITIONAL.tsv, COSTS.json and a short REPORT.md.

## Command contract

`COMPARE.py --stage template --repo REPO --output INPUTS.json` creates a pinned
old-reference template with six explicit missing new slots. This only reads
source/data metadata. Root fills actual paths or pins and prior attempts.

`COMPARE.py --stage bind --config INPUTS.json --output CONFIG.json` freezes
source/spec and all supplied result/receipt pins without numerical analysis.

`COMPARE.py --stage compare --config CONFIG.json --output FRESH_DIRECTORY`
runs the saved-data analysis, requiring an actual Slurm job. One CPU is enough.
No waiting collector or automatic resubmission is implemented.
