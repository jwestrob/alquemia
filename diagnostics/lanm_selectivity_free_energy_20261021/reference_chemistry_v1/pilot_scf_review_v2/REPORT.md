# Native pilot: supplemental SCF review

Saved output only. The frozen native adapter and its original admission remain unchanged. Numeric values below supplement the adapter's empty YES/NO residual field.

| Cell | Native checks | Energy | Max density | RMS density | DIIS | Orbital gradient | Orbital rotation |
|---|---|---|---|---|---|---|---|
| ae_pbe0__acetate_bi_central_La | True | 2.001e-10 (pass) | 2.231e-05 (fail) | 1.195e-06 (fail) | 0.00199 (fail) | 1.516e-05 (fail) | 1.23e-05 (fail) |
| ae_pbe0__water_central_La | True | 2.199e-09 (fail) | 8.081e-07 (fail) | 7.472e-08 (fail) | 0.001271 (fail) | 4.914e-07 (pass) | 1.344e-06 (pass) |
| ecp_pbe0__acetate_bi_central_Dy | True | -1.795e-10 (pass) | 7.506e-06 (fail) | 5.628e-07 (fail) | 0.001535 (fail) | 8.617e-06 (fail) | 7.489e-06 (fail) |
| ecp_pbe0__acetate_bi_central_La | True | 3.058e-11 (pass) | 2.087e-06 (fail) | 1.753e-07 (fail) | 0.00272 (fail) | 2.992e-06 (fail) | 2.626e-06 (fail) |
| ecp_pbe0__water_central_Dy | True | 7.738e-10 (pass) | 6.101e-06 (fail) | 5.224e-07 (fail) | 0.003231 (fail) | 1.96e-06 (pass) | 7.375e-06 (fail) |
| ecp_pbe0__water_central_La | True | 3.158e-11 (pass) | 7.3e-07 (fail) | 7.49e-08 (fail) | 0.001627 (fail) | 1.955e-07 (pass) | 2.738e-07 (pass) |

Tolerances in order: 1e-9, 1e-8, 1e-9, 1e-8, 2e-6, 2e-6 (native atomic units).

Completed 6/8; no terminal result yet for ae_pbe0__acetate_bi_central_Dy, ae_pbe0__water_central_Dy.

The four frozen-f cells admit the remaining frozen-f panel as diagnostic coverage under the existing protocol. They do not qualify strict SCF/force accuracy. Explicit-f Dy admission remains separate; do not extrapolate the La results to Dy.

No physical-state, spin-orbit, basis-limit, correlated-reference, derivative-consistency or experimental selectivity qualification follows from this review.
