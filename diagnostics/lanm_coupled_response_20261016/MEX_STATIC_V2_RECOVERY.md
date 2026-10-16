# Mex static replay recovery

The original `MEX_STATIC.py` (SHA256
`9e458e4fe373db6b068b9c6d9b6753a63f391436eb4ac68abc4bad02344977d9`)
and failed job12673 are preserved. That job made zero molecular calls.
Saved-only diagnostic12680 completed all four12570 replays: energies,
populations, raw QM gradients and raw external-site gradients were exactly
equal. Differences were confined to derived cap-force scatter and its assembled
atom gradients/forces, at most `8.881784197001252e-16` kcal/mol/Å. The real-QM,
point-charge and translation-diagnostic fields were exact. Evidence is pinned in
the V2 configuration under `previous_failure` and `replay_diagnostic`.

`MEX_STATIC_V2.py` changes only this saved-derived-arithmetic admission. Raw
fields and their artifact hashes remain exact. Derived keys, dimensions, finite
values and force signs are checked explicitly. Real-QM and point-charge
components, all global translation diagnostics, and all entries outside actual
cap anchors still require exact equality.

At the four caps' retained/omitted anchors, the per-component bound is
`64 * eps64 * max(1, S)` in kcal/mol/Å. For each cap component, `S` is the sum
of absolute operands of its length-three Jacobian-transpose/gradient dot product,
converted once to output units and accumulated over any caps sharing an anchor.
For total gradient/force entries, `S` also includes both cap operand sums and the
absolute real-QM and external-charge components added at that atom. Thus a small
result caused by cancellation does not set the arithmetic scale. Two evaluations
of these three-term dots, at most four cap accumulations, unit conversion and
four-component assembly involve only a few dozen float64 rounding operations;
the factor64 is a conservative arithmetic allowance, not a chemical accuracy
tolerance or a fit to the observed residual. Global reductions are not covered
by this short-operation bound and deliberately retain their exact check.

The adapter writes all per-field residuals and bounds to `SAVED_REPLAY.json`.
It then uses the original accepted mapped gradients without replacing them by
the replayed values; the saved embedded-minus-isolated field reconstruction
remains exact. Exceeding the bound or any unchanged gate stops before new calls.

All original scientific configuration fields must equal12673's configuration.
The only configuration additions are pins `previous_failure`,
`replay_diagnostic` and `recovery_note`, plus the V2 `implementation` pin.
The original static plan, Fortran wrapper, source, models, parameters, states,
charge rules, ledger, solvent and proposed two-MACE/two-ALPB call matrix remain
unchanged. Root owns the new submission. No numerical execution was performed
while writing this recovery.
