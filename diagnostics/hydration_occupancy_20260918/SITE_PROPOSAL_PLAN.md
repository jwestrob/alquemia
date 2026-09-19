# Aligned experimental water sites — preparation only

Part of the approved occupancy plan: propose missing waters from aligned
experimental structures, using identical proposal rules for both metals.
Use only the two already consumed 1F6S/6IP9 structures. No new biological control,
energy evaluation, inserted production water, occupancy label or result selection.

Fit the common, source-mapped nonwater protein heavy atoms retained in the fixed
expanded contexts, by equal-weight least-squares proper rotation/translation.
Transfer each experimentally observed variable-water O in both directions.
Report all positions, fitting RMSD, closest actual target water, closest actual
nonwater/nonmetal heavy atom, and distance to the target metal. Scan the target
source structure, not just the small core, when identifying clashes/waters.

Before inspecting these transfers, use a1.0A nearest-water distance to label an
already represented water site; flag a proposed missing site if any nonmetal
protein/cofactor heavy atom is within2.0A or if its target-metal distance exceeds
3.2A. These are conservative proposal-screening rules, not free-energy/occupancy
criteria or tuned classification parameters. Keep rejected positions in output.
Do not insert any transferred position into the20-endpoint occupancy experiment.
Further computations need separate explicit input inventories; no jobs here.
