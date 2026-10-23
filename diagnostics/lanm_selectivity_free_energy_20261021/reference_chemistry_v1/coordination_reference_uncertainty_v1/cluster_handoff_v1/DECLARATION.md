# Eight consumed cluster energies: versioned execution handoff

Source/preparation only. Root owns release and submission. Preserve the original
26-cell design, its failed Tight18 controls (32/35), all tighter-control attempts,
and every original input and adapter. No new donor, CP fragment, optimization,
force, basis, geometry, field, held-out outcome or fit is included.

## Scientific inventory and prospective staging change

Reuse the eight full-real complexes already in original prepared_v2/PANEL.json:
`cluster__{tight,tighter}__aquo{8,9}_central_{La,Dy}`. Coordinates, order, charge+3,
effective singlet and physical-state metadata, ECP46/55, orbital/AuxC assets,
NoFrozenCore, RHF, VeryTightSCF/MaxIter300, local-MP2/localization, CC and iterative
triples settings remain exact. The original deck and recovered native parser
are reused. CC contains no D4. All native admission criteria are unchanged.

The original first cluster was **Tight W9 Dy**. This handoff prospectively changes
the first/only initially releasable cell to **tighter W9 Dy**, the largest cluster
at the recipe that now passes all18 consumed canonical-overlap controls and35
unchanged checks. Original Tight remains a lower-accuracy sensitivity comparator;
its three failed controls are never relabeled. This change is declared before
any cluster call. No existing cluster native START/output/result was found.

After the actual scout passes native admission, an independent actual review
and root acceptance of measured cost, memory and scratch are required before a
separate release of the other seven exact cells. Each cell has one immutable
execution directory and at most one native invocation. A failed/unavailable
scout cannot release the remainder. No automatic retry or alternate method.

## Resource change only

Original resources were four MPI ranks, one thread each, default memory/time,
with template MaxCore3000 MiB/rank and actual MaxCore=floor(.75*grant_MiB/4).
New resources are **24 MPI ranks on exactly one node**, one thread/rank,
standard partition/standard_expand QoS. No memory or wall-time override.
MaxCore=floor(.75*actual granted MiB/24), whether the default grant is reported
per node or per CPU. Require a positive grant and matching Slurm rank/node
counts; retain the grant, MaxCore, input, command and complete artifacts.
Template decks change only nprocs4→24; the runtime changes only MaxCore to the
actual allocation-derived value. No claim of linear 8→24-rank speedup is made.

## Complete comparison and limitations

For each tier use all four energies:
`C = E(DyW9) - E(LaW9) - E(DyW8) + E(LaW8)`.
The original cutoff-sensitivity criterion remains
`abs(C_tighter-C_tight) <= 0.25 kcal/mol`, evaluated only after all eight native
admissions. Retain all four endpoint shifts. Report both contrasts against
the same saved PBE0 electronic energies and separately PBE0+D4; never add D4
to CC or mix new CC totals with old PBE0 fragments into a CP/B cycle.
The existing PBE0 electronic C is -0.640844592 kcal/mol; the original final
PBE0+D4 C is -0.534206667 with D4 contribution +0.106638333 (small printed-stage
closure differences remain numerical evidence, not a new target adjustment).

Small-fragment local/canonical agreement qualifies this recipe only on those
consumed controls. Even passing cluster cutoff sensitivity does not establish
basis convergence, physical frozen-f/SOC/multiplet accuracy, forces, dipoles,
collective B, sampling or affinity. No physical-accuracy threshold or target
replacement is introduced. All printed SCF criteria, T1/localization warnings,
native residuals, failed attempts and costs remain explicit.

## Existing cost evidence

The actual tighter40-electron/201-orbital donor complexes took4242–5575 s on
eight ranks (9.43–12.39 allocated CPU-h), with reported batch MaxRSS about
32,767,000 KiB. Water-dimer complexes took152/161 s on eight ranks. W8 has88
correlated electrons/391 orbitals; W9 has98/434. There is no measured local-CC
cluster runtime or memory bound, so those small timings cannot honestly predict
the scout. The old canonical engineering projection was104–111 four-rank hours
per W9 and1295 CPU-h for the four canonical clusters; it was neither measured
nor a local-CC budget. Actual largest-scout resource evidence governs continuation.
Parent scheduler allocations must be counted once; application sums cover only
admitted rows and do not replace accounting of failed/unavailable allocations.

## Release/review contracts

Source review: `passed_cluster_handoff_source_review`, exact implementation,
panel and preparation pins. Scout actual review:
`passed_cluster_handoff_scout_review`, exact implementation, panel, collection
and result pins, `native_checks_passed=true`, `cost_reviewed=true`.
Root continuation review: `accepted_cluster_handoff_scout_feasibility`, the same
collection/result, exact independent scout review, and
`cost_and_memory_accepted=true`, `scratch_reviewed=true`.
These later authorities remain null until actual results exist. The prepared
bundle and tighter18 root admission alone authorize no cluster execution.
