# q9 saved-motion structural diagnosis

Question: does the saved large accommodation response coincide with changed
metal coordination or with contacts/environment constraints still absent from
both finite contexts? Geometry can locate candidate interactions; it cannot
assign a unique energetic cause.

Inputs are the pinned q9 source parent PDB, original/expanded context atom maps,
origin and the two saved La proposals in the completed factorial PREPARATION.json.
No energy calls, optimization, new structure, label access or thresholds fitted
to scores. Coordinate arithmetic is limited to one numerical thread.

Definitions fixed before numerical inspection:

- Use `scripts/coordination_policy.py` chemically typed direct donors and its
  **3.1 Å** coordination threshold; retain its **4.0 Å** diagnostic shell.
  Report distances as well as the binary count and denticity by source residue.
- Report every changed inter-residue real heavy-atom pair involving a moved
  atom within **4.0 Å** in at least one saved geometry. Descriptive distance bins
  are **≤3.5 Å**, **3.5–4.0 Å**, and **>4.0 Å**. These are contact descriptors,
  not energy or clash admission criteria. Caps are excluded from real contacts.
- Separately flag N/O donor–H···N/O acceptor geometry using explicit parent H:
  D···A≤3.5 Å, H···A≤2.5 Å and D–H···A≥120°. Assign parent H only to an N/O
  within1.3 Å and the same residue, preserving ambiguities. Report raw distances
  and angles, and avoid asserting a hydrogen-bond energy. Protein carbonyl and
  carboxylate O, alcohol O and unprotonated histidine N are eligible acceptors;
  amide/backbone/Arg nitrogens are not acceptors.
- Identify all source atoms omitted from both contexts within4.0 Å of any
  moved atom in any saved geometry, including donor, backbone and explicit H.
  Parent atoms outside the moved set remain at their real source coordinates.
- Retain the existing0.8 Å heavy-displacement admission boundary and recorded
  angular/physical feasibility flags; this is not an unconstrained minimum.
- Inventory only the known donor-response DFT campaign and same-source frozen
  fold DFT manifests, following their explicit pins. Exact geometry, atom/state
  inventory and method compatibility must be demonstrated before calling an
  archived DFT endpoint a validation of these motions.
