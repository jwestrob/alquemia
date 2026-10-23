# Frozen water candidate: balanced eight-to-nine-water development screen

Declared before new candidate evaluation. Root owns release and submission.

Use the exact frozen six-water-parameter candidate from job17605, the SI water
Hamiltonian, and the original unoptimized `aquo9_central_La/Dy` geometries.
This is consumed development evidence against provisional frozen-f ECP-PBE0-D4,
not an independent validation set or a hydration-population prediction.

## Question and four fixed checks

Does the candidate transfer from its eight-water training anchor to the existing
nine-water configurations, with acceptable metal-dependent work and rigid forces?

The sole energy check is model-minus-reference error in

`C = [E(Dy,W9) - E(La,W9)] - [E(Dy,W8) - E(La,W8)]`.

Require `abs(C_model - C_reference) <= 1 kcal/mol`. The free-water term cancels
between metals. Do not compare raw different-element energies or infer a
coordination equilibrium from this contrast. Original W9 also changes the first
eight water placements (oxygen z=±0.875 Å versus W8 ±1.25 Å); this is a balanced
coordination-and-geometry transfer test, not an isolated ninth-water insertion.

At W9, separately require La, Dy, and Dy-minus-La rigid-response RMS <=2 and
maximum body-component error <=5 kcal/mol/Å, using the unchanged category-balanced
net-force and centroid-torque/1Å helper `matched_dry_v1/OBJECTIVE.py` (bf65a20a).
All four checks are required; missing/failed cells remain unknown and cannot pass.
No scalar fitting objective, threshold change, optimization or parameter fitting.

## Exact reuse and preparation

- Re-admit both saved W9 and W8 ECP reference E/F results from
  `reference_chemistry_v1/ecp_remaining_v1/tasks/ecp_pbe0__aquo{8,9}_central_{La,Dy}/execution/`.
  Their original REFERENCE.py/panel, native input recipe, states and coordinates
  must match. For each metal the W8/W9 native input decks are byte-identical;
  the separately pinned coordinate files differ. Charge remains+3, effective
  multiplicity1, physical La1/Dy6, ECP cores46/55. Explicit electrons88/98.
- Reuse exact frozen17605/17612 W8 model anchors and their raw-output admission
  through the already reviewed water-only model helper. Do not rerun W8.
- Copy the original typed W9 XYZ and bonds from explicit_model_v1/prepared;
  replace only key/prm with the exact corresponding frozen17612 W8 key/prm bytes.
  Only types500/501 and36/37 occur, so no donor override is active. Preserve all
  original W9 coordinates/order and SI water intramolecular terms.
- Check the17605 entire vector ledger and17612 row inventory contain no W9;
  preserve original pre-repair W9 AMOEBA outputs as a different historical model.

## Finite native execution

One root-released1-CPU job, standard_expand/default memory/time, one native
thread. Exactly two cells, each `analyze system.xyz E` followed by
`testgrad system.xyz Y N`: maximum4 native programs, no new Cartesian FD or DFT.
The existing small-system native invocation, normal-completion and gradient
parsers are reused. Verify executable/source/input hashes at runtime. First
native or admission failure stops subsequent starts; no automatic retry.
Record all receipts and finite arrays, retain both cell slots, and do not sum
overlapping per-program child-CPU receipts. The saved scorer runs only after
native completion and replays available raw E/F before evaluating the checks.

A pass advances this frozen SI-water candidate's consumed coordination screen.
It does not qualify BIO flexible water, salt association, periodic dynamics,
equilibrium hydration, protein interactions or exchange free energy.
