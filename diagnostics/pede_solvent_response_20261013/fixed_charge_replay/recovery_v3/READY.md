# XYZ reader convention recovery v3

Slurm4001 built the isolated native library and evaluated all four water-charge
diagonals. The gate rejected all four hydrogen-bond and SASA matches; no vacuum
charge evaluation ran. Actual Ghb residuals were4.83–8.44e-9Eh; SASA residuals
were−1.35 to−1.99e-9Eh. Gelec and reference-shift checks passed. The raw results,
rounding envelopes and failed qualification are retained in v2.

Source inspection establishes an input-unit mismatch in the wrapper:

- Legacy `xtb_mctc_convert` defines `autoaa=0.52917726`.
- The actual executable's `xtb_io_reader::readMolecule` delegates to
  `mctc_io::read_structure`.
- `mctc_io_read_xyz` uses `mctc_io_convert::aatoau`, derived from its CODATA
  constants. The pinned mctc-lib version is0.4.2, also the installed dependency.

V3 imports the actual reader's conversion module for the unchanged saved XYZ
coordinates. This is the only executable-line change. Native solvation
parameter constants retain their original definitions; there is no global
constant replacement or parameter adjustment. The actual source chain is
pinned in `technical_recovery`. Whether this resolves every replay residual
remains an execution result, not a claim from static inspection.

Manifest:
`workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3/manifest.json`

SHA256 `d556ea7a7e60a8b503f44281f40404e93cb7fd11ede26b91aa946bb253a9ec1c`.

Root submitted the unchanged finite build/gated replay as Slurm4002 on the
healthy standard node,12 task CPUs, with actual completion watcher armed.
`SUBMISSION.json` records root's command and watcher. No compilation or
scientific execution occurred on login. No SCC/DFT call is introduced.

The reproducible command (already submitted; do not duplicate) is:

```bash
REPLAY_WORK=/groups/banfield/projects/environmental/sr/srvp2020/Jacob/lanthanide_binding/on_density_scanner/alchemical_bvs/workspaces/pede_solvent_response_20261013/fixed_charge_replay_v3
REPLAY_SHA=d556ea7a7e60a8b503f44281f40404e93cb7fd11ede26b91aa946bb253a9ec1c
sbatch --nodelist=node-64-768g-5 "$REPLAY_WORK/implementation/RUN.sbatch" \
  build-test-execute "$REPLAY_WORK/manifest.json" "$REPLAY_SHA" \
  "$REPLAY_WORK/execution_v1"
```

Every charge, source-state, component expectation, rounding bound, informative
ceiling,1e-9Eh allowance and four-water-before-four-vacuum gate is unchanged.
The previous four inadmissible diagonals are rerun under the corrected input
convention; no accepted SCC, DFT or vacuum-charge component result is repeated.
Root owns actual result inspection and interpretation.
