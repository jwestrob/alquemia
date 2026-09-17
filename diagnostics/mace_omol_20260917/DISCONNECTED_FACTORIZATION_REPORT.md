# Two bound forwards reproduce the four-forward descriptor

The saved-output audit passes60checks. Largest difference2.3574e-8 model kcal,
against the predeclared0.01 tolerance. It covers all five primary development
scores, GGR sodium, and all detached repeat/rotation/farther variants. Zero new
model/DFT/solver/gradient/training calls were made.

With the charge feature masked, the unchanged protein readouts cancel between
the two disconnected states. The fixed1H4I detached metal terms are Ca
-18430.794927644074 and La -850.2720512362149 model eV. They are actual native
node plus embedding readouts, not quantum ion/aquo energies or a fitted offset.

R_mask = [T_bound,Ca - T_bound,La - (C_Ca-C_La)] *23.06054783061903.

The optional prepared-input `--factorization` path therefore needs only two
bound forwards, with an explicit verified reference. The four-state path remains
available and the running canonical job is unchanged. Both output types record
their evaluation ID; no four-state components are fabricated in the two-call
report. Classification remains unavailable until compatible calibration exists.

Real1H4I interface preparation now reuses exactly two actual bound receipts and
reports the equivalent scalar. Three factorization tests pass44.009s, including
corrupted-reference rejection and no invented detached components. No new
inference was launched to validate the interface.

Cross-environment replay initially refused the v1audit due to six auxiliary
NumPy sum differences, <=2.1477e-8 model kcal. All reference constants, scores,
pass/fail outcomes and the tolerance were unchanged. V2uses math.fsum for stable
component summation; exact verification now works in both pinned environments.
Original report, differences and failed interface timing remain archived.

Authoritative current proof: workspaces/mace_omol_20260917/factorization_report_v2/
result.json. Frozen source: factorization_source_v2/implementation/.
Two-call interface: prepared_interface_pqq_two_call_v2/; completed report:
prepared_interface_pqq_two_call_report_v2/. Larger panel equivalence remains
pending. Read FACTORIZATION_NUMERICAL_BANDS_ADDENDUM.md before publishing any
execution-specific reference. No new accuracy observation or default promotion.
