# Separate geometry-response utility from full entropy qualification

The first completed coupled check (1F6S11Ca) gives response-direction energy
errors0.00030/0.00131kcal/mol, whereas its soft even-energy error is approximately
-0.00975kcal/mol, exceeding the frozen0.005kcal/mol floor. That failed curvature
check remains failed; do not loosen its tolerance or release harmonic entropy.

A surrogate-guided geometry step can still be useful when its tested descent
path and actual native trial energy/gradient support the step. Recenter eligibility
therefore requires the originally frozen `response` direction checks and native
proposal energy check to pass. Keep the independent soft-mode result alongside it.
Each subsequent proposal is checked against new native DFT, and unsupported
steps stop; the two-round and0.20Å/0.35rad limits remain unchanged. Full-spectrum
thermochemical qualification still requires the soft/coupled checks and supported
basins. This separates relaxation from entropy, as required by the research brief;
it does not promote a failed curvature test or refit a biological threshold.

This decision is made during consumed method development, after the first four
native direction outputs, before any recentering run or occupancy estimate.
