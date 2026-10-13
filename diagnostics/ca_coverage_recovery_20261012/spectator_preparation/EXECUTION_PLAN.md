# Spectator-aware preparation implementation v1

The geometry-only Q9 audit admitted all five explicitly two-Ca sources. This
finite first execution freshly prepares those five and replays the previously
scored one-Ca P12293 sample0. No electronic endpoints, model forwards, new
optimization, calibration or classifier decision occurs in this stage.

The original PLAN.md/AUDIT.json remain the pre-implementation record. This
implementation follows that plan with the first execution order explicitly
set by root: one-Ca regression first, then all five Q9 sources. The independent
1FLG structural-mapping audit has already run in the refold work and is reused
as motivation; no duplicate crystal protonation is required here.

`freeze` performs only path/hash/JSON binding. `execute` requires a verified
RUNNING Slurm allocation before scientific imports. Use five one-thread source
workers after the sequential one-Ca regression (six allocated CPU slots).
Workers use explicit Python spawn, avoiding inherited initialized OpenMM/native
thread state. The unchanged protonator internally positions hydrogens; no new
accommodation or electronic optimization is introduced.
There are no resource/time/memory flags beyond these shared CPU slots. Every
source ends with a separate ADMISSION and an explicit unavailable status on
failure. The execution directory is exclusive-created; no automatic rerun.

`PROTONATE.py` retains the original seeded, CPU-one-thread, forcefield=None,
standard-residue-only PDBFixer policy and exact canonical protonator. Its sole
source-occupancy extension admits one monatomic Ca spectator aliased CA in
addition to the active native Ca aliased LA and the one PQQ. No metal changes
native element in the full source. No H is admitted on those nonstandard
residues. Missing heavy-atom repair is rejected. The original helper is intact.

The original source-to-core and fixed-envelope preparation functions are reused.
Only the local call's helper binding supplies the new narrow protonator; no
shared module or production function is changed. Active selectors come from
actual per-sample audits, including Q9 sample3 active D/spectator B. The core
uses its original PQQ(3−), dry state, homologous roles, typed contacts and
Ca/La charges. Original canonical La-source fragment unions and their exact
state/protein signatures are mandatory. No adaptive expansion or water changes.

Before scoring eligibility, both real Ca must remain in the full normalized
and protonated sources, all source heavy atoms must retain identity and
coordinates within the existing 0.001 Å PDB serialization precision, the
spectator must be absent from all local source/cap/bond maps, and its distance
to **every retained physical heavy atom** must exceed 4.3 Å. Record the actual
minimum, active/PQQ/anchor distances, and full-source pins. A distant spectator
is excluded explicitly from the local Hamiltonian: this does not model its
whole-protein electrostatics or mechanics. Any near spectator is unsupported.

The one-Ca replay must match the actual archived core and initial context at
1e−12 Å, plus the original expanded state, paired mapping and coordinates under
the unchanged replay check. Failure blocks the new preparation. A later batch
may reuse its exact passed regression receipt with identical code/config;
all other source chemistry is prepared freshly. It cannot reuse a failed or
unrun regression.

`FINAL_PREPARATION.json` retains every declared source. New two-Ca IDs remain
separate from original one-Ca IDs. Admitted rows carry `case_id`,
`original_source_id`, `selection_id`, `source`, `representations.context`,
`maps`, `state_signature`, `state_key`, `frozen_state_key`, `state_anchor`,
`anchor_case_id`, `mapping_eligible`, `paired_maps_equal`, and the
`spectator_admission` pin. Its receipt uses `status: passed, passed: true` only
on success. The top-level `canonical_reference` is the unchanged corrected
canonical collection; no new threshold is fit or inherited silently.

After actual preparation root/another agent may bind these rows to the existing
warm CPU proposal and corrected solvent executors. Preparation itself makes no
prediction. The original 125-source benchmark and its unavailable rows remain
immutable; these are separately named source-occupancy recovery results.
