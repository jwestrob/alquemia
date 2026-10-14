# 1FLG A/B: unchanged accommodation, then frozen-charge endpoints

Declared before new energies. Bind only actual terminal two-chain preparation
records from the separate crystal preparation branch. Its source and chemical
limitations remain explicit, including the conditional oxidized PQQ state and
any separately documented remote standard-terminal completion. This adapter
changes no source geometry, protonation, region, parameter, optimizer or gate.

One seven-CPU Slurm job keeps the existing float64 OMOL worker warm for both
independent crystal chains. Reuse the exact original CPU transport, four-mode
force-selected SLSQP settings (one start per metal, existing maxiter200), actual
two qualification fixtures, source maps and common origin/Ca-proposal/La-proposal
pool. Both metals are cross-scored on each admitted geometry. No second search,
new start, threshold adjustment or automatic restart is introduced.

At most four origin calls, four original optimizer starts and four cross calls,
plus two actual runtime qualification calls, are planned. Optimization may make
multiple actual calls under its unchanged finite policy; report those calls,
not this endpoint count as an observed cost. No GFN2, ALPB or DFT calls occur in
this adapter. The root owns their separate finite solvent stage.

Recorded source admission and a Slurm physical preflight are both required.
The preflight verifies the paired coordinates/states, original physical
kinematics at q=0, unmodified kernel/runtime/settings and the exact historical
fixture energy/force gates. A source-specific failure keeps that chain
unavailable while a separately valid chain may proceed. A failed shared model
qualification blocks all new sources. Both chains remain in every output.

Each actual common-pool cell exports XYZ, state/counts and its existing MACE
receipt/energy in `MATRIX_INPUT.json`. Root's frozen-charge endpoint scorer must
reuse these MACE values. Missing required cells stay explicit; no static-score
fallback or fabricated control fills a gap. No new calibration or production
promotion is part of this preparation/transport adapter.

Preparation is metadata-only and can bind a completed real artifact:

```bash
python diagnostics/exaa_source_packing_20261014/crystal_scoring/NATIVE.py prepare \
 --preparation ACTUAL_CRYSTAL_FINAL_PREPARATION.json --output NEW_NATIVE_WORKSPACE
```

After root review, the exact execution entry is:

```bash
sbatch --output=NEW_NATIVE_WORKSPACE/run_%j.out \
 --error=NEW_NATIVE_WORKSPACE/run_%j.err \
 diagnostics/exaa_source_packing_20261014/crystal_scoring/RUN.sbatch \
 NEW_NATIVE_WORKSPACE/CONFIG.json
```

Root supplies real absolute paths and arms the usual completion watcher. No
test or molecular call has run when this file was written. Scientific imports,
preflights, model calls and saved scientific exports execute only on Slurm.
Started partial calls and original failures are preserved; the fixed first
attempt directory prevents accidental execution restart. The separate `export`
command can collect existing source results to a fresh named file on Slurm
after an interrupted job without evaluating the model.

## Frozen-charge solvent continuation

`SOLVENT.py` reuses the pinned root `MATRIX.py` preflight, full xTB/tblite
vacuum/water component gates, and native ALPB water-diagonal/vacuum-charge replay.
It verifies and reuses every actual MACE receipt: **zero new MACE calls**.
It does not edit the running packing implementation. At most twelve new
geometry/state cells give 48 SCC calls and 24 native ALPB evaluations. Their
costs are separate from the earlier model search. No ORCA default parameters
are used; the exact published-parameter standalone runtime stays pinned.

An incomplete native pool remains unavailable and does not dispatch solvent
for that source. A complete independent chain can proceed. Both chain rows
remain in the final denominator; failed solvent cells prevent pooled selection.
Where both origin solvent cells succeed, the existing scorer retains its
independent origin policy even when another solvent cell fails.

Only exact unchanged-Q9 local-state preparations receive the separately frozen
operational, mathematical and origin reference bands from canonical collection
SHA256 `ee00fa03aba64c23e6a737fb2b3b4a9f40dc139e58bdb9225e8ed10c9af2b58d`.
There is no refit. Report these as conditional source comparisons, with crystal
redox unverified and both chains kept separate. They are not fresh blind tests
or evidence of physiological occupancy. The larger packing contexts continue
to receive no transferred class bands.

The root can queue the continuation behind the actual native job:

```bash
sbatch --dependency=afterany:NATIVE_JOB_ID \
 --output=SOLVENT_PARENT/crystal_solvent_%j.out \
 --error=SOLVENT_PARENT/crystal_solvent_%j.err \
 diagnostics/exaa_source_packing_20261014/crystal_scoring/SOLVENT.sbatch \
 NEW_NATIVE_WORKSPACE/MATRIX_INPUT.json NEW_SOLVENT_WORKSPACE
```

The Slurm wrapper binds actual terminal native outputs before running solvent.
A missing terminal artifact fails visibly; it never reruns the native search.
Root must create only the log parent in advance and leave the new solvent
workspace absent for immutable creation. All scientific work remains on Slurm.
