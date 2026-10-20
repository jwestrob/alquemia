# Independent static review

**PASS for the declared zero-call Slurm preparation.** No material source or scope blocker remains. The two native singlepoints still require the actual successful preparation and its bound root execution decision; this review is not a numerical result.

## Frozen evidence

| Artifact | SHA256 |
| --- | --- |
| DIRECTION.py | `cd4396c02cf50905307acbf90772d062ae4a6a12aaf9914e01562167dca1ab32` |
| CONFIG_PREPARE.json | `cec2f6e393f191c656f7dd6fd57441895c14f2bb5fd6a75014608f2b5428f43f` |
| SPEC.md | `c79d67370983f6100bc9e3ce3fea45a5641af26a4f924fdc3af41efd84f8c609` |
| READY.json | `f495a9811205b4b5decd2f7707f5addab4ea3f31bb97ab566a59bd0769769dd0` |
| Preparation launcher | `f76554920719a2d1118dab4454e0d06a3242ce2ecf8f4f2288d641f0172b4890` |

All eight direct CONFIG artifact pins and all four READY artifact pins were hash-verified. Actual central native/force metadata and the relevant original runner, collector and reference helper interfaces were read. No scientific modules were imported or numerical values recomputed.

## Source and admission checks

- The original `f251ca1e` inputs and saved-run admission bind native 15545, collector 15838 and the exact step 725 force report. Saved physical IDs, coordinates, charge/state, raw analytic gradient and MACE vacuum force are checked before constructing a direction. Both arrays retain all 256 atoms.
- The direction is the normalized full Cartesian MACE-minus-PBE0 vacuum force difference, without alignment or projection. Displacements use the declared ±0.01 Å full-vector norm. Fresh physical metadata, positions and XYZ are validated against that construction; singleton signatures must match the preparation. The only native deck change is removal of EnGrad.
- `single_point` and `engrad_path=None` match the original runner interface. `ggr_sensitivity.executed` validates input/output/receipt provenance and treats absent expected engrad as appropriate for an SP. New admission requires normal termination, native SCF convergence, one finite final energy, ORCA 6.1.1, whole charge −2, singlet/RHF, 1030 explicit electrons, La ECP46, D4, PModel and vacuum without external charges. It preserves the original SCF stopping semantics and diagnostics rather than imposing an all-residual pass.
- Native execution remains one worker with 12 MPI ranks and allocation-derived memory. The preparation launcher uses one CPU, the pinned original Python environment and one thread per library. The preparation and collection branches contain no molecular runner call; each released singleton contains one runner invocation. Root owns the global two-call limit and fresh job scheduling.
- Both declared slots initialize unavailable. The final source records the requested task and preserves failed pre-selection, zero-call parents; complete parents still require exactly the selected task, one call and exact saved native admission. Collection requires both available endpoints to produce a derivative verdict. There is no retry path.
- Decimal arithmetic subtracts the two final native energy strings before conversion. The negative centered energy slope is compared with the central PBE0 force projection using the frozen `0.005 + 0.0001*max(abs(F_FD), abs(F_analytic))` rule. Printed-energy rounding is reported separately. Geometric bookkeeping tolerance is not substituted for the scientific consistency gate.

## Interpretation

The test addresses the existing reference's energy/force consistency in one deliberately selected direction at one step size. MACE selects that direction; a failed PBE0 self-consistency check does not itself establish a MACE defect. Failure cannot uniquely distinguish finite-step, grid, SCF or electronic-state effects. A pass does not qualify every force component, choose between the different training/reference Hamiltonians, establish a physical minimum or validate solvent, preference or transfer.

The central endpoint's failed printed energy/density/DIIS residual flags remain preserved despite native gradient-based SCF termination. No central rerun, changed native settings, new solvent/model call, optimization or automatic follow-up is admitted by this review. Actual Slurm preparation and subsequent native receipts remain decisive.
