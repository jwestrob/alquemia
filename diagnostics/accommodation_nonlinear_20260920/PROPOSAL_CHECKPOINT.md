# Proposal experiment checkpoint

2026-09-20, 12:53 UTC. Root reviewed the immutable manifest and implementation
before submission. The fixed [proposal plan](PROPOSAL_PLAN.md) is unchanged.

**MACE proposal job 1204169 completed all 60 starts successfully.** It made
283 actual native energy/force calls in 86.678916 executor wall seconds on one
H200: 2,773.725317 allocated core-seconds and 86.678916 allocated GPU-seconds
for the measured executor interval. Full scheduler allocation costs remain
separate. All origins passed replay and retain 25/25 canonical and 3/3 crystal
calls before accommodation. No composite-selected result is available yet.

The generated primary GFN2 manifest contains exactly 120 new proposal endpoints,
zero missing proposals and zero q0 reruns. Its actual runner dry-run passed.
**Solvent job 1204171 submitted**, 64 CPUs, eight concurrent eight-rank tasks,
128 GiB requested memory. Its wrapper collects even if the runner exits nonzero;
failed members remain unavailable, with no substitution of the origin result.

Run directory: `workspaces/accommodation_nonlinear_20260920/proposals_v1/`.

- `manifest.json`, SHA256
  `edaf320f6da3725c9270f8e5c33d1046c6558623a630d428d929eb771f1b0b28`.
- `GPU_SUBMISSION.json` and `GFN2_SUBMISSION.json`: actual commands/receipts.
- `proposal_execution_1204169.json`, `gpu_summary_1204169.json`: measured work.
- `proposals/<endpoint>/`: actual optimizer traces, energies, Cartesian forces,
  active gradients, boundary flags, geometry and failed-attempt fields.
- `after_proposals_1204169.json/.md`: honest pre-solvent collection; selected
  scores remain unavailable.
- `proposal_GFN2/manifest.json`, SHA256
  `4acb0ef8f0ff5abefc08e651bb3af99885fb835c843bdfb0bb860c1c0eac4f26`.
- `proposal_GFN2/PREFLIGHT.json`: 120-task dry-run.
- `final_1204171.json/.md`: expected actual automatic collection after the
  solvent job; this checkpoint does not imply that those energies exist yet.

Six real-artifact preflight tests passed, zero skips. No DFT, extra optimizer
starts, threshold changes, production changes or new biological labels.
The large PLM vacuum proposal work is a mechanistic observation; await actual
composite selection before interpreting the score. Commands are in
[PROPOSAL_COMMANDS.md](PROPOSAL_COMMANDS.md).
