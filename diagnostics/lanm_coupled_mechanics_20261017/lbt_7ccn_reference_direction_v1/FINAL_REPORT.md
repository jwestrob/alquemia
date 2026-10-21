# Small reference inconsistency; large model disagreement remains

21 October 2026 UTC. Both declared whole-peptide singlepoints completed and
passed unchanged native admission. Saved comparison 16582 **fails the frozen
directional consistency gate**. That failure is retained without a retry or
tolerance change.

The sampled PBE0 energy slope is nevertheless close to its analytic force and
far from MACE's force. This small inconsistency does not account for the large
model/reference disagreement in this direction. The remaining attribution
question includes the different electronic methods used for MACE's training
target and this PBE0 comparator.

## What was measured

The exact model-converged 7CCN La structure at step 725 contains 256 physical
atoms, whole charge −2 and physical/effective singlet. The two new positions
are `x ± 0.01*u` Å, where `u` is the normalized full Cartesian MACE-minus-PBE0
vacuum force difference. The 0.01 Å is the norm over the entire coordinate
vector; there is no alignment, atom selection or optimization.

Both new ORCA 6.1.1 calls retain the original PBE0-D4, def2-TZVP/def2-J,
RIJCOSX, DefGrid3, VeryTightSCF, PModel and special La basis/ECP46 recipe.
Removing EnGrad requests only the two singlepoint energies. The central
analytic gradient from 15545 is reused. No new solvent or MACE calculation
occurred.

| Quantity | Value |
| --- | ---: |
| Minus energy, Eh | −6765.683736715161 |
| Plus energy, Eh | −6765.681619165231 |
| Force from centered energy difference, kcal/mol/Å | −66.4391321372 |
| Central PBE0 analytic force projection, kcal/mol/Å | −66.4531999104 |
| Central MACE force projection, kcal/mol/Å | −9.02704252585 |
| FD minus analytic residual, kcal/mol/Å | +0.0140677732 |
| Frozen allowed absolute residual, kcal/mol/Å | 0.0116453200 |
| **Directional gate** | **Failed** |

The collector subtracts the printed energy strings using Decimal arithmetic
before unit conversion. Its reported energy-printing rounding bound is
3.13754737e−8 kcal/mol/Å, well below the observed residual. This does not bound
SCF, integration-grid, state or finite-step errors.

The predeclared rule remains
`abs(F_FD − F_analytic) <= 0.005 + 0.0001*max(abs(F_FD), abs(F_analytic))`.
The observed residual exceeds it. The collector's Slurm `FAILED`, exit `1:0`,
records this scientific outcome: both endpoint rows are available, all joins
reproduce, and there is no collection exception or missing endpoint to repair.

## Interpretation and limits

- The strict force-consistency qualification remains false. One direction
  and one step do not qualify all components, an electronic state or a minimum.
- The energy-derived force remains close to the PBE0 analytic value, while
  MACE's projection is much smaller in magnitude. The measured small
  inconsistency therefore does not explain this large directional discrepancy.
- This does not prove that PBE0 is physically correct or uniquely identify a
  learned-model defect. MACE-OMOL targets ωB97M-V/def2-TZVPD with different
  lanthanide spin conventions; exact method equivalence has not been shown.
  See the [method evidence](HAMILTONIAN_CONTEXT.md) and
  [OMol25 calculation details](https://arxiv.org/html/2505.08762v1#S2.SS6).
- Neither electronic comparison independently validates the shared fixed-charge
  solvent contribution that drives the recent relaxation benefit. There is no
  new metal-preference, natural-protein transfer or production claim.

Both native runs converge by the original gradient-based SCF check after 18
reported cycles. Orbital-gradient and orbital-rotation thresholds pass; printed
energy, maximum/RMS density and DIIS thresholds do not all pass. Both retain the
early small-gap warning. These observations were preserved under the original
admission, not converted into an all-residual convergence or stability claim.

## Decision

Close this finite PBE0 directional test with its failed gate preserved. No
automatic grid sweep, half-step retry, central rerun or further MACE relaxation
follows.

The next contained step is [verification of the actual training-reference
recipe](../lbt_omol_target_reference_v1/SOURCE_REVIEW_SCOPE.md). A comparison
against that target could help distinguish comparator-method differences from
the remaining model/target discrepancy. Exact basis/ECP, spin, guess, numerical settings and executable
lineage must be established before separately declaring any native scout.
Source verification alone releases no calculation.

## Cost, review and preservation

| Job | Work | Scheduler outcome | Allocated CPU-seconds |
| --- | --- | --- | ---: |
| 15845 | Saved-input preparation | Completed, 6 s × 1 CPU | 6 |
| 15847 | Minus native singlepoint | Completed, 21,623 s × 12 CPUs | 259,476 |
| 15848 | Plus native singlepoint | Completed, 21,616 s × 12 CPUs | 259,392 |
| 16582 | Saved comparison | Scientific gate failed, 7 s × 1 CPU | 7 |

This finite phase used **518,881 allocated CPU-seconds**, no GPUs and exactly
two new native energy calls. The central gradient and earlier model/solvent
costs are recorded in their original phases. Allocated cost is separate from
measured CPU use.

- Minus RESULT: `8e73f55a229e89bd2cec3f7b9e1ec4553ff82d46417b9936e18cc2454294c5db`.
- Plus RESULT: `9cedffa6dcd0ce1a53dd8cea5b29a78b11477f8c4b8a3053e33d5dcd5f37ae32`.
- Collector RESULT: `1f2e53ffb590ae766f3230c037519220218ccb6bb17241d8168523619145eb81`.
- Original source: `cd4396c0`; collection config: `3b958db0`.
- [Independent native review](ACTUAL_NATIVE_REVIEW.md): `d87d880f`.
- [Independent collection review](ACTUAL_COLLECTION_REVIEW.md): `3c8bc24f`.
- [Native evidence archive](native_results_v1/ARCHIVE.json): `e86624f0`.
- [Saved comparison archive](collection_result_v1/ARCHIVE.json): `7c67a8c2`.

Native evidence, complete endpoint results, inputs, runtime receipts, SCF
diagnostics and actual accounting are preserved. Large density/restart files
remain remotely pinned; compressed text copies are lossless. Original failures,
frozen protocols and production are unchanged.
