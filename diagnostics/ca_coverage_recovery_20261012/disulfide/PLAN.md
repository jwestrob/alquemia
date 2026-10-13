# MMOL_2048 disulfide preparation recovery

Question: does the rejected Ca-conditioned seed1/sample0 contain a preparation
mapping defect, or does its actual cysteine chemistry fail the preserved model?
The old `unsupported disulfide partner` message is not a diagnosis by itself.

The fixed union includes CYS127. Existing code obtains SG–SG connectivity from
OpenMM, then requires two CB hydrogens and no SG hydrogen on both cysteines.
The partner-H inventory check runs before the existing 1.8–2.3 Å S–S distance
check. SourceGraph assigns hydrogen parents by nearest heavy atom within each
source residue. These are concrete hypotheses to inspect, not reasons to relax
the gate or assume the original fold is defective.

## First finite diagnostic

Three pinned real sources from the same protein:

- Failing Ca-conditioned seed1/sample0.
- Accepted Ca-conditioned seed1/sample1.
- The original fixed-state La-conditioned seed1/sample0 anchor.

One Slurm CPU runs `DIAGNOSE.py`, reports actual atom coordinates, hydrogen
parents/nearest-heavy distances, SG neighbors and source topology bonds, and
replays the unchanged cysteine and complete-expansion gates. No energies,
hydrogen generation, new state selection or structural relaxation occurs.
The exact source/preparation/runtime pins and command are in
`workspaces/ca_coverage_recovery_20261012/disulfide/diagnosis_v1/PLAN.json` and
`PREPARATION.json`. Root owns submission and completion notification.

## Decision after the real diagnostic

If a mapping defect is established, implement a separately versioned repair in
this owned directory. Preserve source heavy coordinates, physical membership,
proton inventory, formal charge, the fixed union, paired Ca/La coordinates and
the original state anchor. Require existing bond, overlap and physical mapping
checks. Compare the real accepted controls before handing preparation to scoring.

If the chemistry or source geometry is genuinely unsupported, report the exact
atoms and defect and the specific source requiring a new fold. Do not manufacture
a missing partner, accept an unphysical bond, silently change thiol/disulfide
state, broaden the context, or tune a repair to a score. Original failed
preparation and all frozen canonical bands remain unchanged.
