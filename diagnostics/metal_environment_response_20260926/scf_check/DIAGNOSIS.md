# A concrete convergence-policy hypothesis — 27 September 2026

All four refined-grid A/rigid calculations stop while their final SOSCF
iteration still changes the density much more than the printed TightSCF
thresholds. The final energy/gradient criteria have not forced small density
steps. This makes SCF stopping policy a justified isolated diagnostic. It does
not establish the cause of the remaining rigid energy difference.

| Endpoint | Stop reason | RMS density change | Maximum density change |
|---|---|---:|---:|
| Ca A | Energy check, cycle 25 | 1.8519e−4 | 3.9824e−3 |
| Ca rigid | Gradient check, cycle 29 | 1.8613e−4 | 5.7868e−3 |
| La A | Energy check, cycle 34 | 2.6123e−4 | 5.2923e−3 |
| La rigid | Energy check, cycle 34 | 1.8905e−4 | 5.1912e−3 |

Printed density tolerances are respectively 5e−9 and 1e−7. Ca rigid also stops
with energy change 4.4914e−8 Eh versus TolE 1e−8. The machine-readable companion
pins each actual output and retains all final residuals and printed tolerances.
DIIS errors are retained but can be stale after switching to SOSCF; they must
not be mistaken for the active terminal orbital gradient.

The [ORCA 6.1 SCF manual](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/scf.html)
distinguishes mode 2's energy checks from mode 0's broader criteria. Mode 0
still permits some overachievement exceptions, so changing the keyword alone
will not qualify the resulting wavefunction. Explicit residual checks remain
necessary. ConvForced prevents following calculations after failed convergence;
it does not independently tighten the convergence definition.

Proposed finite diagnostic: four Ca/La A/rigid endpoints, same refined grid,
Hamiltonian, electronic state, geometry, initial-guess policy and iteration
limit; change only SCF check mode to 0 and force convergence before derivatives.
Retain the existing TightSCF numerical thresholds rather than introducing
another accuracy sweep. Audit energy, density and active orbital criteria from
actual terminal output. Apply the unchanged rigid energy/gradient gates.
No B response, optimization, numerical gradient, new state or grid sweep belongs
in this diagnostic. Original failures remain recorded regardless of its result.
