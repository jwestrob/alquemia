# First whole-protein LanM result — 24 September 2026

Job1213018 ran458s on1H200/32CPUs (14656 allocated core-seconds,458 GPU-allocation
seconds). It returned a failed overall status because all four native ORCA
launches failed before computing energies. The dependent1213040 was cancelled
without running; do not reinstate its automatic eight-system continuation.

## MACE executed, with a physical failure in both proposals

All129 attempted MACE calls completed in451.420s, including native/adapter force
comparison, both bounded searches and the admitted origin pair. Peak framework
GPU allocation was73,678,407,680bytes; native and adapted energy were identical,
with maximum force-component discrepancy6.84e-14eV/Å. This establishes working
full-source energies/forces, not successful accommodation or selectivity.

Both searches reached60iterations and were rejected by the frozen covalent-bond
gate. La elongated Asp85CB–CG and Glu91CG–CD to2.705/2.725Å; Dy elongated
Glu91CG–CD and Asp107CB–CG to2.864/2.596Å, from approximately1.52Å. They retain
C-alpha chirality and have no new severe heavy clashes, but lose acceptable
carboxylate connectivity. Both have hundreds of Cartesian components at the
movement boundary. Their large vacuum-energy drops are not admitted as physical
accommodation. A shared pool containing only the origin is not an improvement.

This failure concerns the unconstrained Cartesian proposal policy on the tested
source. It does not establish that a covalent-geometry-preserving whole-protein
search or within-lanthanide discrimination is impossible. No revised search
model is launched here. All proposed coordinates and actual forces are retained.

## Native startup failure and contained technical recovery

The batch requested one Slurm task with32CPUs, while each ORCA call requested
eight MPI ranks. OpenMPI saw insufficient task slots and aborted startup. This
was our launch-layout error. The native Hamiltonian and its feasibility have
not been tested by those failed starts. ORCA returned exit0 despite error
termination; the existing receipt checks correctly marked every task failed.

The recovery prepares those same four origin inputs in a fresh directory,
byte-identical coordinates and templates, unchanged native eight-rank recipe,
charge, effective multiplicity and numerical settings. A CPU-only allocation
requests32 Slurm tasks×1CPU, with four concurrent eight-rank calls. No
oversubscription or partition-rule exception is used. No MACE computation is
repeated. Successful native outputs will be collected separately; this recovery
does not restart the remaining eight molecular systems or claim a valid
accommodated score.

Recovery inputs: `workspaces/lanm_global_occupancy_20260923/native_feasibility_retry_v1/manifest.json`.
Original MACE output: `workspaces/lanm_global_occupancy_20260923/mace_feasibility_v2/result.json`.
Original failed native outputs and all cost records remain preserved.
# Technical recovery submission

Job1216461 repeats only the four unchanged native origin cells, using32 Slurm
tasks ×1CPU and four concurrent eight-rank calculations on standard/memory.
No GPU, new MACE calls, DFT or continuation to other systems. Input and XYZ hashes,
scientific keys, chemical states and pinned implementation match the original
tasks; the four-task dry-run, collector syntax and shell syntax checks pass.
`native_feasibility_retry_v1/manifest.json` SHA256:
`81b3a6eba18a7efc728f84e10dbc2fdb3b78b885d9fe10e17736071226ee0f65`.
Submission is not successful molecular execution; consult its actual receipts.
