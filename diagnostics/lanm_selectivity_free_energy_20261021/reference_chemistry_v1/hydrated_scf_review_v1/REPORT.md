# Hydrated native-output review

All six mixed-donor and both reused aquo8 references pass their original native
checks. Their saved numeric residuals are retained in RESULT.json. The native
SCF stopping logic is supported; no retroactive all-six-printed-threshold gate
is imposed. See ../scf_mode_interpretation_v1/REPORT.md for the actual SOSCF
trace interpretation and a finite observable-based precision test.

The hydrated comparison remains a conditional fixed-geometry electronic test.
It has not qualified method/basis/electronic-state accuracy or ensemble free
energies. Five of six mixed-donor outputs stop through the energy criterion
while printed gradient/rotation values exceed their separate listed thresholds;
mono-acetate La also meets those listed values. This motivates a bounded
precision comparison if required by the target accuracy, not an unsupported
claim that these calculations failed.
