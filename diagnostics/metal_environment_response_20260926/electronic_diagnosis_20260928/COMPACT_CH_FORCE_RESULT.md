# Compact metal-load diagnostic

Conditional fixed-core projected gradients, not affinity or relaxation work

Positive gradient load means moving the metal toward its preselected nearest oxygen raises energy. Force has the opposite sign. Units: kcal/mol/angstrom.

| Source | nearest O index / distance Å | La gradient | Dy gradient | Dy−La |
|---|---|---:|---:|---:|
| Hans8DQ2 | 8 / 2.452468 | 28.5549 | 19.0087 | -9.54614 |
| Hans8FNR | 8 / 2.292688 | 63.5997 | 49.3438 | -14.2559 |
| Mex8FNS | 5 / 2.351439 | -56.9408 | -56.9815 | -0.0406902 |

Full endpoint status, gradient norms, oxygen distances and capped-system translation residuals are retained in JSON. No force is subtracted to impose translational invariance. Directions are source-specific, so their load differences are not a common-coordinate finite difference between different sources. No physical protein/link-force projection or state/affinity qualification is claimed.
