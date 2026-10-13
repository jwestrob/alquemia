# MMOL_2048 preparation and native pool recovered; solvent score pending

**The excluded source contains an intact disulfide. A malformed existing hydrogen
coordinate caused the preparation failure; a new protein fold was unnecessary.**
The versioned one-H repair now passes the unchanged199-atom envelope, physical
mapping and frozen chemical-state checks. Native accommodation now completes
the common six-cell geometry pool; no new composite classification is claimed here.

## What actually failed

Slurm2864 inspected the rejected MMOL_2048 Ca seed1/sample0 source, accepted
Ca sample1 and the original La sample0 state anchor. The source topology contains
the real Cys127–Cys128 bond at2.117557 Å. Cys128 HB2 was present, but its CB–H
distance1.376989 Å exceeded the source graph's actual1.340000 Å parent-assignment
cutoff. The resulting missing H-parent entry triggered `unsupported disulfide
partner` before the S–S distance test. Both controls passed without alteration.

## What changed and passed

Slurm2870 restored only that existing HB2 along its original CB–H ray to the
pinned ff19SB CYX equilibrium1.09 Å. The explicit XML parameter is
`protein-H1`/`protein-2C`, length0.10900000000000001nm; topology SHA
`086bfdb1c05e7d5cb1d330e6c39fab5056a95e89492cdbd3044c746db13e9d09`.
After PDB serialization, CB–H is1.089972477 Å and H displacement0.287017421 Å.
Every other coordinate and atom identity is unchanged. This is radial
preparation repair, not an angular optimization or a different protonation state.

The original source/protonation records remain immutable. New records point to
the repaired source, retaining the original core XYZ and fixed fragment union.
No builder cutoff or admission rule changed. Actual results:

- Target:199atoms, paired maps equal, exact original small-context replay.
- Frozen source-state key matches the existing canonical group:
  `f8eeff480ce4ee77f17fbccba0dd5f458922d734abd2320c94b3c997ee2cf2eb`.
- Two untouched controls reproduce their archived states and physical maps.
  Maximum coordinate-copy residual is2.22e−16 Å; the La control's XYZ hashes
  are identical, while Ca control serialization differs only at this rounding scale.
- Core coordinates, heavy geometry, proton/water inventory, charges, selection,
  canonical bands and original failed benchmark row remain unchanged.

Protocol: `Nikasha_disulfide_unmapped_source_CH_radial_repair_v1`.
Real receipt: `workspaces/ca_coverage_recovery_20261012/disulfide/repair_v1/RECEIPT.json`.
Result SHA: `147aedd094e2a43837c519d7a6a53fa7bbd3a7a7f5ac6260b43793b741d2b367`.

## Work and next step

Diagnosis:4.335s measured script wall time. Preparation:7.875s wall,
6.664processCPU-s, peak RSS177300KiB on one Slurm CPU. No energies, forces,
optimizations, new protons or waters were calculated.

Slurm2875 completed the one-source native stage on seven CPUs at
node-64-768g-5: both origins, both unchanged four-mode proposals and both cross
evaluations form a complete shared pool. The two existing runtime force fixtures
passed with exact energies and maximum force residual9.55e−14 eV/Å. This stage
used33 actual forwards,463.764s elapsed,2549.236 measured modelCPU-s and
3246.348 allocatedCPU-s. These are this CPU recovery's measured costs, not
GPU throughput or a full new-site latency estimate. All started model calls
have terminal receipts. The same-protein preparation controls were not rescored.

`SOLVENT_PLAN.md` defines corrected solvent Slurm2880, submitted after completed
2875 on node-64-768g-5 with12CPUs and completion watcher PID3530160. It uses unchanged
component gates and frozen bands, at most12 fresh inputs per backend. It retains
the125 primary rows and appends preparation recovery, with separately pinned
MxaF technical recovery when available. Solvent/composite outcomes remain
unavailable until actual execution. A repaired input or native pool is not
evidence of better class discrimination.

Original125-source results, frozen calibration and every original failure remain
visible. Later results must be appended as named preparation recovery.
