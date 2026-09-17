# Audit the vacuum counterpart of the failed reaction-field term

Declared after full-boundary GB job1200975 completed and failed both partition
and all-four ordering gates. This audit adds no score, calibrated class or
scientific execution; it analyzes the same frozen real charge vectors.

For all four representations and both endpoints, compute the point-charge
vacuum cross interaction C_M=k sum_ij q_QM,M,i q_env,j/r_ij from the actual
projected QM distribution and fixed exterior. There is no charge overlap by
construction; reject a zero distance. Use the existing project constant
332.063713299 kcal A/(mol e^2), retain its source pin and all input hashes.
Report C_Ca-C_La alongside the actual GB cross contrast, their sum, and the
archived MACE full-minus-core contrast. Do NOT label that sum a corrected hybrid:
MACE's implicit electrostatic contribution is not uniquely separable, so adding
bare C to its context could double count interactions.

Bin the direct contrast by exterior atom distance from the metal at
[0,6,12,18,24,30,36,infinity]A; retain per-residue and per-atom terms. Boundaries
come from the pinned6A graph cutoff and3interactionblocks, not labels/results.
18A bounds changed readout centers;36A is a conservative bound for an external
input to influence the metal-substitution contrast through a shared readout.
These two bounds are not interchangeable. The shared charge-feature mask and
fixed atom inventory remove endpoint-dependent global charge conditioning.
Do not assume all electrostatics inside36A are represented accurately.

Checks: component/bin/residue/paired algebra closes within1e-7kcal, repeat the
same direct sums after the existing joint rotation/translation within1e-7.
No geometry, atom, charge, cavity or solvent-state changes. This quantifies a
possible cancellation problem; it does not identify a unique cause for failure
or validate a forcefield/ML energy split. No new DFT,MACE,GB,training or cluster
allocation. Record actual local wall/CPU/memory and preserve the failed parent.
