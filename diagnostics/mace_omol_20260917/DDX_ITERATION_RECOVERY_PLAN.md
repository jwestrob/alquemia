# Complete the declared PCM checks with independent solver states

The unchanged 2FW0 La coarse replay completed: the dielectric system needed
337 iterations, followed by 79 single-layer iterations. It satisfies the same
1e-10 iterate-change criterion; energy is -40.85096540267961 kcal/mol. Its
264.183-second solve establishes slow convergence for this state, not a general
convergence or accuracy claim. Original 300-iteration failures remain failures.

Complete the original eight groups / twenty state roles with one explicit
numerical recovery version. Keep all physical sources, grids, cavity parameters,
FMM orders, charge representation, source evaluation and acceptance tolerances
unchanged. Use maxiter=1200, native logging, and a fresh Model for each newly
executed state. A failed Model retains an error flag in the Python interface;
it must not prevent the next endpoint/control from being independently tried.

Reuse every computed state from the completed direct-source job1201279, and
the matching logged La result1201280. Verify input/software identity and retain
their original execution parameters, coefficients, receipts and failed attempts.
Prepare the exact missing-state inventory only after the original job finishes.
No duplicate new solve for a reusable state. A prepared group may execute
independently using the existing 64-CPU/128-GiB allocation policy. At most eight
group jobs, one serial solver at a time within each; no automatic retries.
No project time/compute budget. This remains development cost, not production.

Collection must distinguish reused and newly computed states, actual forward
starts, failures and unattempted work. Combine coefficients only after verifying
identical cavity points. Cross-endpoint reciprocity, refinement, rotation,
repeat, zero-source and contraction checks remain exactly as frozen in
DDX_SOURCE_SELF_PLAN.md. A numerical iteration setting cannot change a label,
physical state or scientific acceptance threshold. No full hybrid score,
calibration, DFT/MACE call or baseline change is authorized by this recovery.

The separate matrix-storage experiment remains deferred. Its cache affects the
second linear system and was not needed to recover the first failed state.
