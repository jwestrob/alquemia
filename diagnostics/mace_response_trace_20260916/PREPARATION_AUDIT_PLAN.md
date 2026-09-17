# Protein bond geometry integrity audit

Read-only follow-up under the active goal. Use the same pinned 1H4I physical
state and the already used Amber ff19SB protein topology, without generating
new atoms or moving coordinates. Compare every protein bond length to its
existing force-field equilibrium length. Report hydrogen/heavy bond groups,
all deviations and the largest absolute deviations, identifying any stretched
bonds at atoms with previously observed large MACE force changes. This is a
preparation diagnostic, not validation of metal/cofactor parameters or evidence
that classical bond lengths uniquely define the correct quantum geometry.

No metal/PQQ force field is assigned. Their atoms remain in the source physical
inventory but are outside this standard-protein bond comparison. Zero new
MACE/DFT calls, no minimization and no score modification. Write all bond rows
and exact source/software hashes to a new workspace directory.
