# LBT electronic/MM partition preparation

17 October 2026 UTC. **Both real LBT sources now have a complete baseline and expanded partition.** This is preparation and geometric mapping qualification; no molecular energy or force evaluation has been performed for these new partitions.

| Region | Full residues | Real electronic atoms including metal | Link H atoms | External charges |
| --- | --- | ---: | ---: | ---: |
| Baseline | 3–14 | 170 | 2 | 84 |
| Expanded | 2–15 | 208 | 2 | 46 |

All four source/region preparations pass, defining eight La/Dy states. The electronic and whole-system charges are both−2; the shifted field sums to zero within recorded roundoff. Full donor chemistry, hydrogen-bond selection, cap-anchor Jacobians, charge redistribution, native bonded/cross terms, exclusions and pure12–6 cross interactions are explicit in PREPARATION/MAPPING/LEDGER artifacts. No extra induction term was added.

Original job12950 failed four common-MM serialization comparisons before energies. The sole difference was target sigma in an MM-only LJ force where target charge and epsilon are both exactly zero. V2 compares that inactive attribute separately after proving the target has no relevant exceptions/offsets. Actual generated La/Dy MM XML must remain byte-identical to12950, including the physically active, metal-dependent QM–MM interactions. All four exact-XML checks pass. The original failure is preserved.

Recovery12953 completed4scheduler seconds/1CPU (application2.854s); failed12950 cost5seconds/1CPU. No molecular calls, optimization or quantum calculation. Complete evidence: `workspaces/lanm_coupled_mechanics_20261017/lbt_partition_v2/execution/RESULT.json`; exact recovery proof: `RECOVERY_V2_EVIDENCE.json` in this diagnostic directory. Native parameters, software and source hashes are pinned in CONFIG_V2.json and deployed configuration.

Next: bind these partitions to the reused electronic/MM/solvent executors and qualify the complete derivatives. Constructed maps alone do not qualify the coupled model.
