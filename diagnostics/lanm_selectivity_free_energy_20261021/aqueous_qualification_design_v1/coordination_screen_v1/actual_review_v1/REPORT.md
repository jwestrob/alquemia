# Independent actual review: frozen water coordination screen

Job **17670** completed normally, and independent replay confirms the scientific
screen **fails three of its four unchanged checks**. No new molecular calls,
parameter changes, or control reruns were made in this review.

## Actual result

The balanced contrast is
`[E(Dy,W9)−E(La,W9)]−[E(Dy,W8)−E(La,W8)]`.

| Quantity | Result |
|---|---:|
| Frozen model contrast | −5.485437060 kcal/mol |
| Provisional ECP-PBE0-D4 contrast | −0.534206667 kcal/mol |
| Model minus reference | **−4.951230393 kcal/mol: fails ±1** |

Independent net-force and centroid-torque arithmetic reproduces all 57 body
components per force comparison. The ion and water categories retain equal
weight, and torque is divided by the declared 1 Å length.

| W9 force comparison | RMS | Maximum component | Original criteria |
|---|---:|---:|---|
| La | 3.811021 | 7.961232 | Fail |
| Dy | 4.685648 | 9.949310 | Fail |
| Dy minus La | 1.054138 | 2.267946 | Pass |

Units are kcal/mol/Å; limits remain RMS ≤2 and maximum ≤5. The passing paired
force check shows cancellation of part of the individual errors on this fixture.
It does not repair the failed finite energy contrast or either individual force
comparison.

W9 also changes the first eight water placements. This tests **coordination and
geometry transfer together**; it cannot isolate ninth-water insertion, predict
coordination populations, or identify a unique missing physical term. The
frozen-f ECP reference remains provisional: effective singlets, La/Dy physical
multiplicities 1/6, cores 46/55, and 88/98 explicit electrons for W8/W9 are
preserved. No reference accuracy or physical-state qualification is inferred.

## Replay and cost

All four new native outputs, both saved W8 model anchors, and all four saved W8/W9
reference energy/gradient results replay exactly. Input bytes, actual producer
authorities, identical per-metal W8/W9 model keys and reference decks, full atom
counts, and finite gradients pass **139 saved checks**. Independent energy and
body arithmetic reproduce the saved comparisons. All four comparison slots remain
present.

Native cost: **3 scheduler seconds ×1 CPU =3 allocated CPU-seconds**; four native
programs used 1.970710 measured aggregate child CPU-seconds and 2.075364 summed
native wall seconds. The completed saved review used one thread, 0.768265 wall
seconds and 0.689720 CPU-seconds; an earlier strict provenance check stopped after
0.120742 wall seconds before scientific replay.

## Historical provenance defect

The release pins pre-submit `ROOT_REVIEW.json` at SHA `62605349…`, but that path
was overwritten by root's terminal review. Root confirms no exact original
snapshot is known. The failed strict attempt and observed terminal bytes are
preserved under the actual-review workspace. **Full historical release-review
integrity is unavailable**; this review does not reconstruct or relabel it.
The release, adapter, bound inputs, native outputs, and terminal result pins
remain intact and support the independent scientific replay above. No rerun is
needed to restate or conceal this documentation defect.

## Exact authorities

- Actual result: `d8de1f0ba01852762f6a012ca241b7072603140eecf68ba07947110d3f04f34d`.
- Native adapter: `d295ab62e125ccdb478a379b8f4d784a2c80c424398399ff6df0f8a1822df88a`.
- Bound inputs: `9310cdd2e6ad278c20aa77461ecb3d470a629c270e04602080267a886200aa9d`.
- Root release: `52d9b1c18b5f07276b8178f76a8bf1f4f73f5627194607dd0a7dd1db5480c9dd`.
- Independent receipt: `6e8e2b4e5f9bb761ad8bbb80e68283724addbe7e471c6351b82dd074f37e0a4f`,
  at the matching workspace `actual_review_v1/REVIEW.json`; it contains full
  source/input/output pins, force components, cost, and the provenance exception.

Preserve the failed screen and frozen candidate. This result authorizes no
automatic refit, bulk/protein sampling, or claim about exchange free energy.
