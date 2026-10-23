# Explicit continuation of the four-parameter radial water solver

## Question and immutable start

Can additional bounded joint steps bring the same four-parameter, eight-cell development problem inside its existing controls? Job17589 stopped at its declared12-step limit while still improving; it did not establish either an optimum or family impossibility. This new, explicitly declared solver-budget test changes only the starting point, allowed additional work and exact prior-cache admission.

Start at the frozen17589 selected vector `[3.800640492160888, 0.27995568058029474, 3.490554961303394, 0.2649340807140388]`, ordered La water Rmin/damping, Dy water Rmin/damping. Its RESULT SHA is `31cb5caa7f96a3f312cc3367964149923665faa16b41d8993e664ad22eeb259c`; objective is0.4970127555469721. The result must be complete, error-free and exactly bound to original radial source `e808ae45f65380255cb02dbace2865932972873989ab3ee8782ccb5f050b22e0` and evaluator `1cf6025e27291523ac7c7242a07245da22c8927d67f4323c87392974133eda56`. The four remaining failures are individual RMS on water minus/plus; paired, energy and original-center checks pass. No prior failed check is removed.

The original eight consumed cells, reference E/forces, coordinates, body definitions and objective remain exact. Training RESULT `3aefeb74f40aaf51e327fdc5c0ffd85f7f00e6e769b07c080f3ece16c345192d`, objective `5b83cf5f3e77aeb96ec0f4d0a0c3f765883f406e8f4110648972e5bac0996f6b`, body helper and original four-center diagnostics are reused without modification. Water48/CC sensitivity remains separate from the uncorrected PBE0 objective. No donor, reserved24, native-reference, geometric or energy-surface expansion occurs.

## Unchanged physics, search and stopping

Keep the original four radius/damping bounds, fixed metal-water epsilons, charges, polarization, water-H, bonded and permanent-multipole terms. Original SI radius values remain the bound origins; the17576 baseline values remain the data/cache baseline descriptor. Neither is the new search start. Actual initial search/model/score must replay the exact17589 selection with zero calls.

Retain the four equal energy/La/Dy/paired loss blocks, all individual/paired fixture RMS≤2 and maximum component≤5 kcal/mol/Å, all nine energy controls≤1 kcal/mol, and every previous center check. No tolerance or weight changes.

At most **32 additional outer steps**, **384 new full parameter vectors** and **4096 new native program invocations**. Each step retains eight central derivative slots at±0.005 normalized coordinates, clipped to the original box with actual-width differences; the unchanged saved joint direction helper applies the0.025 normalized trust box. Evaluate scaled-GN half/full then bounded-LS half/full, in that order, completing the full derivative and trial batches before selecting the lowest original objective. Current wins ties. Repeated full vectors and identical effective cell inputs reuse prior work.

Stop immediately after a complete selected update if all unchanged controls pass, if a complete trial batch makes no actual objective improvement, or at32 additional steps. First native/admission failure stops new starts, retaining sibling receipts and the eight-cell denominator. No fallback, epsilon change, further continuation or transfer is automatically released. A controls-pass stop is not a global optimum or physical qualification; a budget/no-improvement stop is not proof of family impossibility.

## Cache, cost and execution contract

Use a fresh workspace, ledger and native cache. Preserve17589 and all older results, native records and ledgers byte-for-byte. Admit only the parent's complete native cells with exact source/helper/executable/geometry/rendered-key/parameter identities and raw E/force records. Index those immutable records by the new request's effective physical input; retain their original native authority and source identity. Never relabel an old record as a newly executed record. An unchanged opposite-metal subvector must reuse its cells. The selected eight starting cells must all reuse with zero native calls.

New native calculations retain one serial analyze plus analytic testgrad per new physical cell. The conservative ceiling is128 programs/step×32=4096:64 across the eight one-coordinate derivative vectors and64 across four full trials. The full-vector ledger permits385 including the initial point. Existing completed cells reduce actual cost; no repeated chemistry is required to initialize the continuation. Execution uses one8-CPU Slurm task, at most8 native processes and one numerical thread per process, standard_expand with default memory/time. Aggregate child CPU is recorded once over the whole driver; overlapping per-receipt CPU deltas are not summed.

Before root release, run source controls and actual saved zero-call preparation: exact17589 selected model/score replay, original family/criteria, read-only inherited-cache and opposite-metal reuse, and new independent limits. Root owns release/submission. Old source, results and ledgers remain protected. This is a finite solver-budget test, with no claim of protein transfer, bulk-water calibration, affinity or electronic-state qualification.

## Runtime-only checkpoint write reduction

The growing full driver RESULT is saved at each batch start and end, and immediately on a failure. It is no longer rewritten after each successful future within a batch. Every cell native receipt, native result and evaluator ledger keeps the existing durability. All slots are declared before starting; complete-batch selection, fixed ordering and failure semantics are unchanged. Incomplete running snapshots are not completed scientific results.
