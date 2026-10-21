# Independent actual directional comparison review

**The declared numerical consistency gate fails. Both native inputs are available and the comparison was evaluated; this is a scientific test failure, not a technical unavailable result.** Preserve the frozen threshold and no-automatic-retry stop.

## Exact result and admission

- Collector 16582 RESULT: `1f2e53ffb590ae766f3230c037519220218ccb6bb17241d8168523619145eb81`.
- Executed DIRECTION.py: `cd4396c02cf50905307acbf90772d062ae4a6a12aaf9914e01562167dca1ab32`.
- CONFIG_COLLECT: `3b958db04d7def0f2d1f270b57f6089a4a21a51cd0a1ee82497bd95ce9d65d4e`.
- Frozen DIRECTION: `f20be9142797fbcbe6f6a468380dc617880f534de1bbf2d64644fc6d9a7165c4`.
- Minus parent: `8e73f55a229e89bd2cec3f7b9e1ec4553ff82d46417b9936e18cc2454294c5db`.
- Plus parent: `9cedffa6dcd0ce1a53dd8cea5b29a78b11477f8c4b8a3053e33d5dcd5f37ae32`.

Seven selected collector/input pins were independently verified. Both collected native rows exactly equal their admitted singleton parent rows, including state, receipt, energy string and SCF observations. Scientific binding matches the actual preparation; the direction hash is unchanged. Stored central projections match the frozen direction. Both declared rows are available, no exception or traceback is recorded, and every molecular/runner/optimizer counter is zero.

## Reported comparison

The following values are copied from the executed decimal-string comparison; no energy difference or projection was recomputed on login.

| Quantity | Reported value |
| --- | ---: |
| Minus final energy, Eh | −6765.683736715161 |
| Plus final energy, Eh | −6765.681619165231 |
| Plus-minus energy, Eh | 0.002117549930 |
| Centered FD force, kcal/mol/Å | −66.4391321371518410 |
| Central PBE0 analytic projection, kcal/mol/Å | −66.45319991036233 |
| Signed FD-minus-analytic residual, kcal/mol/Å | +0.0140677732104890 |
| Frozen tolerance, kcal/mol/Å | 0.011645319991036233 |
| Printed-energy rounding bound, kcal/mol/Å | 3.13754737e-8 |
| Saved MACE vacuum projection, kcal/mol/Å | −9.027042525850963 |

Top-level and comparison statuses are both `failed_directional_consistency`, with `passed=false`. The source intentionally exits 1 for this evaluated gate failure. The scheduler event consequently records FAILED/exit 1:0, seven seconds on one CPU. It must not be relabeled as an unavailable native calculation or treated as authorization to rerun the collector or endpoints.

## Consequential interpretation

The strict one-direction, one-step qualification was not earned. Printed scalar rounding is too small to account for the reported residual, but this test alone does not separate finite-step, SCF stationarity, grid/integration or electronic-state effects. The retained native energy/density/DIIS flags remain false despite gradient-based convergence; no all-criteria convergence claim follows.

At the same time, the sampled PBE0 energy slope is near its analytic projection and far from the saved MACE projection. This sampled numerical discrepancy therefore does not account for the large model/reference disagreement along the selected direction. The result supports a real difference in the sampled energy responses while preserving the failed numerical gate. It does not identify whether learned-model error or the different electronic Hamiltonian, basis and state conventions dominate that difference.

Only one direction and one displacement magnitude were sampled. No full-gradient accuracy, physical minimum, solvent, metal preference or natural-transfer qualification follows. MACE selected the direction; a PBE0 self-consistency failure is not itself a MACE defect. No automatic half-step, new grid, relaxed tolerance, model repair or optimization follows from this result.

Review used saved fields, literal equality, source text and hashes only. No scientific imports, numerical replay, jobs or source/config changes were made.
