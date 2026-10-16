# Added donor/backbone freedom improves structural consistency

The matched ten-coordinate search reduces the source-local La/Dy contrast
discrepancy to **29.509 kcal/mol**, compared with **40.452** after equal extra
search effort in the original seven coordinates. The starting range pool was
45.805. Thus the extra physical freedoms help beyond continued optimization.
This is useful evidence for modeling scaffold response; the remaining source
dependence is still far too large for a reliable preference prediction.

## What physically changed

Both arms start at the same four saved range endpoints, with identical energy,
bounds, numerical gates and finite search allowance. The candidate adds intact
Asp3/Glu11 side-chain rotations and Ile10 backbone phi, selected from earlier
force evidence. Both arms use the same deterministic feasible-direction fallback.
All original256 atoms remain in the energy. Every new candidate passes the
existing covalent, mapping, contact and state checks. No source was discarded.

The strongest extra response occurs for La in the heavy-source7CCN geometry:
the candidate adds−19.685 kcal/mol from its starting point, while the matched
control adds−5.830. The response is primarily local MACE energy with a solvent
contribution. The final paired pools cross-score every admitted geometry with
both metals, preserving20 previous cells and adding8 own plus8 cross cells.
All36 matrix cells and all8 arm/state slots are available.

The global common menu still selects a7CCO-derived La geometry and a7CCN-derived
Dy geometry. Its contrast moves−4.318 kcal/mol from control to candidate. This
is a conditional energy response, not an affinity or validation against ITC.
The common global menu removes source selection by construction; the separate
source-local comparison exposes the remaining structural sensitivity.

## Limits and next decision

All three remaining candidate searches reach the12-step allowance; other
searches stop after finite geometric/energy rejection. These are admitted finite
proposals, not stationary minima. Bond lengths, protonation, water inventory,
physical model and source-dependent NMR refinement remain conditional.
Neither a smaller source spread nor favorable response completes the LanM goal.

Preserve this finite mechanics comparison. The next chemical test examines
explicit water near Asn5 using equal water inventories. In parallel, qualify
a local electronic field correction and solvent forces for a whole-protein
mechanics route on Mex-LanM. No automatic further torsion/bound expansion,
production promotion, repeated quantum campaign or PQQ change follows.

## Execution and reproduction

Jobs12556/12558/12560 complete:128/580/99seconds, each12CPUs and oneH200.
This round makes66MACE and132GFN2 calls;9684 allocatedCPU-seconds and807
allocatedGPU-seconds. Cumulative new monomer branch through this round is
19634CPU-seconds and1574GPU-seconds, including its preparation/support failures
and excluding older Mex calculations. Allocation is not measured utilization;
see ACCOUNTING.tsv for scheduler TotalCPU/memory records.

Exact runnable commands, frozen adapters/configs and terminal receipts are in
`workspaces/lanm_coupled_response_20261016/lbt3_basis_{scout,remaining,pool}_v1/`
SUBMISSION.json. They are completed commands, not instructions to rerun.
RESULT.json pins the full workspace matrix and searches; SUMMARY.json supplies
individual menus, works, selections, source spreads and the declared denominator.
