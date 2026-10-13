# MMOL_2048 existing-hydrogen preparation recovery

Job2864 found an intact Cys127–Cys128 disulfide (S–S 2.117557 Å).
The partner exists. Its existing HB2 lies 1.376989 Å from CB and therefore
fails the original nearest-parent covalent distance gate. Both same-protein
controls retain two correctly mapped CB hydrogens. No new fold is warranted
by this defect alone.

## Frozen repair

Use the exact three sources from diagnosis_v1: affected Ca sample0, successful
Ca sample1, and frozen La sample0 state anchor. Consider only unmapped existing
HB2/HB3 atoms of actual source disulfide CYS residues. Require the ff19SB CYX
template to assign CB, CB to be the nearest heavy atom, no thiol HG, the original
S–S gate 1.8–2.3 Å, and the original CB–H distance above the unchanged parent
cutoff but at most1.5 Å. Unsupported cases remain unsupported.

Restore only the overlong bond radially along its existing direction using the
pinned ff19SB harmonic-bond equilibrium: protein-2C/protein-H1,
0.10900000000000001 nm (1.09 Å). This is source-coordinate preparation, not an
energy minimization or protonation change. PDB serialization allows at most
0.001 Å bond-length roundoff; displacement must remain within0.5 Å. All other
coordinates and every atom identity must be exact. No hydrogen is added/deleted.
Do not change parent-distance cutoffs, source graph code, forcefield chemistry,
protein identity, core coordinates, charge, water/PQQ state, or selection.

Clone the protonation receipt and core parent solely to point to the versioned
source. Retain archived core endpoint XYZ. Use existing motion_envelope.prepare_task
unchanged, including old small-context exact replay, frozen union, paired-coordinate,
charge/parity, physical cap/motion mapping and canonical state-anchor gates.
The target must require exactly one H repair. Both controls must require none
and reproduce their actual archived coordinates and physical maps using the
existing compare_probe checker. Failure remains explicit; no automatic retry.

## Execution and interpretation

One CPU Slurm preparation job; no model/DFT/xTB evaluation, no optimizer, no new
calibration. PLAN.json pins source code, actual source artifacts, topology,
diagnosis, existing builders and frozen reference. A successful preparation
would recover input coverage only; score and classification remain uncomputed.
The original unavailable row stays unchanged alongside this named recovery.

```bash
sbatch --parsable diagnostics/ca_coverage_recovery_20261012/disulfide/run_repair.sbatch "$PWD/workspaces/ca_coverage_recovery_20261012/disulfide/repair_v1/PLAN.json"
```

Root owns submission and completion wake. No memory/time/account/QOS flags.
