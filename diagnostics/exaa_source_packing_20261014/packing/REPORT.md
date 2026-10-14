# Saved-motion packing preparation: complete

Both declared ExaA two-Ca sources, sample0 control and sample3 target, now have
complete **249-atom origin/Ca-proposal/La-proposal pools for both metals**.
Actual Slurm job4089 prepared all12 inputs successfully. No new model energy,
solvent calculation, optimization or fold was performed.

## Physical change

The fixed structural rule in [PLAN.md](PLAN.md) retains the old192-atom context
and adds complete source fragments contacted within4.3Å by actual moved donor
heavy atoms across either saved proposal. The fragment-identity union is shared
by both sources and all candidates. It adds:

- Glu213 and Trp349 peptide groups, including their source-mapped bonded partners;
- Phe442, Leu485 and Trp591 sidechains.

Leu485 was not selected by a special exception. Its selection follows the same
all-heavy contact rule as the other additions. Phe442's backbone was already
represented; its sidechain is newly included. Newly included source atoms remain
at their real source coordinates unless they already belong to the archived
physical motion. Source hydrogen positions are reused, including the preserved
old motion; no additional hydrogen preparation was introduced. New sigma caps
follow their real retained/omitted source atoms without independent movement.

This is a context test that can change local electronic, dispersion/repulsion
and solvent contributions together. It cannot by itself attribute a contrast
change specifically to Leu485 or to a classical steric penalty.

## Actual preparation checks

Both original contexts replay exactly. The saved kinematic proposals replay
exactly; origin replay error is at most2.22e-16Å. Every shared real/cofactor atom
and surviving cap is copied exactly into its corresponding expanded candidate.
Retained source covalent-bond lengths change by at most2.89e-15Å across the saved
motions, below the predeclared1e-8Å gate. The existing0.45Å overlap failure gate
passes; minimum all-atom distances are1.01Å.

Both sources retain Ca charge−3/La charge−2, multiplicity1 and the same expanded
state identity. PQQ chemistry and water inventory remain unchanged. The active
Ca is resolved from the exact source selector, including sample3's chainD
assignment; the other Ca remains in the full source. Its nearest expanded
retained-heavy distances are23.7867Å and23.9880Å, beyond the frozen4.3Å boundary.
The newly supported neutral Phe/Leu sidechains pass the pinned ff19SB complete
heavy-membership and per-heavy-atom source hydrogen-parent checks.

## Files, costs and next step

Authoritative result:
`workspaces/exaa_source_packing_20261014/packing_v1/execution/PREPARATION.json`
(SHA256 `2e443ee5ca87b490b9959409d535e0efeec8f95af847f5bd0fd0e348702e9108`).
It contains both declared rows,12 endpoint XYZ pins, charges/electron counts,
per-candidate maps, source/preparation/state pins and the immutable old matrices.
`UNION.json` retains all actual selection contacts and the common membership.

The measured preparation used3.366s elapsed,3.463CPU-s and198740KiB peak RSS on
one CPU (`node-64-768g-5`). Scheduler allocation accounting belongs to root's
job receipt; these are process measurements. Exact command and pins are in
[READY.json](READY.json).

Root can now evaluate the12 new endpoints using the declared warm MACE and
qualified native fixed-charge solvent path, reusing the12 old endpoint cells.
Compare each metal's work and the change in paired contrast within each source.
No classification improvement, calibration transfer or production promotion is
claimed by this preparation result.
