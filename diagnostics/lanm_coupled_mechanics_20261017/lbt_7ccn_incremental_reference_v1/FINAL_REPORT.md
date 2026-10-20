# 7CCN-La incremental reference: energy lowering supported, force discrepancy remains

The admitted whole-PBE0 reference supports lowering the energy along the exact
step80→92 continuation. With the same saved fixed-FF ALPB contribution, PBE0 gives
**−4.726291 kcal/mol**, versus **−4.226978 kcal/mol** for whole-MACE. The
MACE-minus-PBE0 work discrepancy is **+0.499313 kcal/mol**, within the declared
1 kcal/mol reporting resolution. This resolution is neither an uncertainty
estimate nor a native admission criterion.

Native endpoint15321 and saved collector15456 completed. The old step80 La
reference was reused; the single new native endpoint is the already selected
step92 La geometry. Both endpoints remain nonconverged geometries. This is a
consumed, Lu-conditioned 7CCN peptide fixture, not an unseen-protein benchmark.

## Incremental work

All values below are reported by the saved collector, in kcal/mol. Differences
were taken in native energy units before conversion.

| Contribution | Step92 minus step80 |
|---|---:|
| Whole-PBE0 vacuum | −2.205301 |
| Whole-MACE vacuum | −1.705988 |
| Shared fixed-FF ALPB | −2.520990 |
| PBE0 + fixed-FF ALPB | −4.726291 |
| MACE + fixed-FF ALPB | −4.226978 |
| MACE minus PBE0 work discrepancy | +0.499313 |

## Full256 force evidence

The saved force tables retain all 256 physical atoms. Forces are negative
gradients; RMS below is the RMS atomic-vector magnitude. Units are
kcal/mol/Å. These are endpoint summaries, not evidence of monotonic force
improvement during the continuation.

| Quantity, including fixed-FF ALPB | Step80 | Step92 |
|---|---:|---:|
| PBE0 maximum atomic force | 12.826316 | 7.904442 |
| PBE0 RMS atomic force | 4.204686 | 3.639505 |
| MACE maximum atomic force | 9.739698 | 5.674648 |
| MACE RMS atomic force | 1.987798 | 1.086248 |
| MACE−PBE0 maximum atomic force error | 7.660828 | 7.718445 |
| MACE−PBE0 RMS atomic force error | 3.574975 | 3.579797 |

Both methods have lower endpoint force magnitudes, while their force discrepancy
remains about **3.58 kcal/mol/Å RMS**. The energy-work agreement therefore does
not establish improved force accuracy or a minimum.

At step92, the saved force projections onto the normalized step80→92 chord are
**+0.110539** for PBE0+ALPB and **−0.009955** for MACE+ALPB. The chord uses the
full256×3 Euclidean norm without alignment. These descriptive projections do
not replace full-force convergence checks or define another search.

## Native admission and limitations

The new endpoint passed the unchanged native analytic admission. SOSCF was
initialized and the output reports 15 cycles. Its energy-change residual,
−5.5297×10⁻¹⁰ Eh, meets the printed 10⁻⁹ Eh criterion. The printed maximum/RMS
density, orbital-gradient, orbital-rotation and DIIS residual criteria remain
false; `all_printed_residuals_pass` is false. The small HOMO/LUMO-gap warning
(−0.019 in the native warning) is retained. The reused old endpoint likewise
has `all_printed_residuals_pass=false`. Full residuals, thresholds and warning
text are preserved in the [archive review](collection_result_v1/ARCHIVE_REVIEW.json)
and raw native output.

Physical-state and force-consistency qualification remain false. The reference
uses the conditional free NH3+/COO− termini and five deprotonated acidic side
chains, whole charge −2, La singlet/ECP46 and no spin–orbit treatment. Post-added
ALPB does not polarize the vacuum density; native convergence does not establish
electron binding or a stable ground state. No minimum, affinity, metal preference
or general accuracy conclusion follows.

## Calls, costs and archived evidence

This phase made **one new native DFT call** with 12 MPI ranks. It made no new
MACE, ALPB, GFN, MM or optimization calls. Preflight and collection made zero
molecular calls.

| Allocation | Elapsed seconds | CPUs | Allocated CPU-seconds |
|---|---:|---:|---:|
| Preflight15319 | 4 | 1 | 4 |
| Native endpoint15321 | 24,866 | 12 | 298,392 |
| Saved collector15456 | 6 | 1 | 6 |
| **Reference-phase total** | | | **298,402** |

GPU allocation was zero. These are scheduler allocation costs, not measured
utilization; they exclude the previously recorded optimizer scout and reused
old reference. [Exact accounting](collection_result_v1/COST.json) retains the
raw scheduler records.

- [Endpoint archive and pins](endpoint_result_v1/EVIDENCE_PINS.json): result,
  manifest, directions, both full-force reports, native output/engrad, inputs,
  runtime/execution receipts and terminal health snapshot. Three binary
  wavefunction/density files are excluded explicitly and retained in the workspace.
- [Collector archive and pins](collection_result_v1/EVIDENCE_PINS.json): exact
  source/configuration/launch records and all saved results, directions and forces.
- [Collected result](collection_result_v1/execution/RESULT.json.gz):
  SHA256 `5908cb2d3b3fcb07b8ac764e3ddfb78fbc56beb29196790388cf79213af4ad72`.
- [Native result](endpoint_result_v1/execution/RESULT.json.gz):
  SHA256 `3a173d81072ab269922d17db1596c130e8132bd2f803614b00d2fe59455c10de`.

All archived execution files are lossless gzip copies with original and
decompressed SHA verification. This report transcribes saved results without
scientific recomputation. Completion receipts were captured with delivery
unconfirmed; no delivery acknowledgement was written here. Root owns event
acknowledgement, integration and any subsequent mechanics decision.
