# Native capped-union scout

Prepared for root submission. Physical preparation4539 passed all four endpoints; wrapper build4543 succeeded by linking the unchanged qualified xTB6.7.1 static libraries. The library is serial. No compilation, scientific imports, tests, or native evaluations ran on the login host during this implementation.

## Execute

Metadata-only binding:

```bash
/usr/bin/python3 diagnostics/native_cavity_shielding_20261014/native/RUN.py prepare \
  --preparation workspaces/native_cavity_shielding_20261014/preparation_v1/execution/PREPARATION.json \
  --build workspaces/native_cavity_shielding_20261014/build_v1/execution/BUILD.json \
  --output workspaces/native_cavity_shielding_20261014/native_v1
```

Root submits the copied `RUN.sbatch` with the absolute CONFIG path and printed SHA256. The wrapper requests one shared node and four CPUs, with no memory/time overrides; native and BLAS threads remain one. It records actual scheduler memory and rejects insufficient RAM before dispatch. Native internal arrays are estimated from their source declarations, with25% plus1GiB headroom. The wrapper does not copy or print a second dense matrix.

The first six updates are sample0Ca local, union, cap-removed, translated union, rotated union and permuted union. Repetition and cap-zero vectors reuse each matrix. Only a passing scout permits the remaining three endpoints' union/cap-removed updates: twelve total, with zero SCC/MACE/DFT calls. No separate first-two scout was submitted.

The actual saved charge/state/component-v2/native qualification pins remain required. Local replay, scalar/block closure and electrostatic invariance retain1e−9Eh allowance. Eight-decimal uncertainty applies only to saved core charges; exact added zeros and zeroed caps have none. The unchanged1e−6Eh informative ceiling gates individual and propagated endpoint differences. Paired and double-contrast uncertainty flags are reported explicitly. Native HB/SASA/reference and shape are retained separately; full-protein solvation is not added to a score.

## Pins and review

- CAVITY.f90: `0a61d74b43bfb2bf9f8a3c0afb1e12416900a847ad55796ec6f9db832c5e5380`
- RUN.py: `c2b8201a655dea878ebf99dbc93abbfc52a5c90e55a970a4401c0b6b605b3c81`
- RUN.sbatch: `40a241f1a7c4d61bca24ad2e75716324d32eb2411148b870191bf46e2a21fa4f`

Gauss independently cleared the source and final runner without numerical execution. Actual native integration is pending root's Slurm job. Failed or blocked cells remain explicit; execution and collection receipts retain partial evidence, actual child CPU/RSS, scheduler cost and dispatch-intent versus execution counts.
