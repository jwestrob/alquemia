# Independent scorer review before molecular launch

Read-only review of the root-owned `scripts/lanm_global_occupancy.py`, the frozen
whole-chain plan and actual nine-state preparation. No molecular calls, scorer
edits, alternate geometries scored or new scientific analyses were performed.

## Consequential integration defects identified and corrected by root

1. The initial state loader expected pin records, whereas preparation supplies
   inline state dictionaries. It now resolves the actual sibling `state.json`
   from the pinned mapping and requires exact inline equality.
2. A failed required search originally disappeared from the candidate set and
   could be published as successful origin/reduced-pool accommodation. The runner
   now records `failed_search`; collection leaves accommodated scores unavailable.
   Valid finite geometry rejections remain separately explicit.
3. Collection originally iterated only returned states, dropped unexecuted states,
   and permitted an empty required-set check. It now validates named state/cell
   coverage and iterates all declared states, retaining unexecuted rows.
4. Reused native receipts were recorded by preparation but omitted by collection.
   Collection now replays the actual original task/receipt with parent and
   scientific-key identity checks.

Root additionally corrected text-based SCF parsing, native TolE/mixer verification,
17-digit XYZ serialization and exact-coordinate scientific reuse keys. These
were reviewed in the current code. A reporting-only caution was sent: a
capability-only origin calculation should retain that status in collected output
and must not be described as evidence of accommodation.

## Scientific checks

The composite sign and units are correct: native MACE eV converted to kcal/mol,
plus native GFN2(ALPB−vacuum). Dy−La followed by Hans−Mex at equal occupancy has
the declared positive-relative-La-preference interpretation. Optimizer gradients
are negative MACE forces in eV/Å. Origin-relative objective subtraction leaves
the gradient unchanged. Source water inventories remain separate; there is no
cross-source geometry pooling or two/four occupancy-population calculation.

The Cartesian component box, covalent bond-ratio, source C-alpha oriented-volume
and new heavy-clash gates implement the declared geometry checks. Every one of
the nine real q0 systems passes unchanged and has no initial severe heavy clash.
Physical Dy multiplicities 11/21 fit the pinned model's categorical spin range;
that is implementation support, not a magnetic-ground-state validation.

Seven real-fixture tests pass in 2.689 s (`SCORING_TESTS_v2.txt`): actual inline
schema/all origins; malformed actual state/native spin; coordinate box/nonfinite
coordinates; water bond stretch and source chirality; newly overlapping heavy
atoms; selection/unit algebra using already completed local-core values; staged
real-XYZ missing-cell/state/search rejection. No fabricated successful molecular
outputs are used. The local-core values are numeric regression fixtures only,
not whole-chain candidate results or a new scientific comparison.

Reviewed scorer SHA256 at this checkpoint:
`5933550397a71d4653d7fba896e1b7343b4ae47ccbb2385271bce38fc5726f72`.

No remaining issue identified that blocks the declared feasibility calls. Actual
GPU/native feasibility, solved-state receipts, cross-cell geometry/state and
complete-pool results still require execution; this review does not substitute
for those checks. Root owns all scorer edits, immutable snapshots and execution.
