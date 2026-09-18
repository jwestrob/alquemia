# Refine the same conductor-like boundary; retain failed coarse tests

The first conductor inventory completed20/20native solves in285allocated wall
seconds. Its numerical report passes420/430checks and fails ten, including
endpoint refinement, one endpoint rotation and one coarse reciprocity check.
The GGR2FW0-minus-2FVY source-self contrast is1.907,2.098,1.892kcal across the
three resolutions, compared with18.129kcal for the archived GK component.
This large model dependence is a useful lead, not a qualified new score or
proof of which approximation is more accurate.

Run a new refinement version of exactly the same conductor model. Preserve
all physical inputs, radii, charges, epsilon, energy prefactor, smoothing,
FMM orders, native solver settings and numerical acceptance thresholds.
Change only the prescribed angular basis and integration resolution:

| Role | lmax | Lebedev points |
|---|---:|---:|
| Coarse reference, reuse actual v1 refined outputs | 12 | 590 |
| New primary | 18 | 974 |
| New refined | 24 | 2030 |

Keep the same eight-group/twenty-role inventory. Reuse the four actual12/590
original-geometry endpoints. Execute sixteen new roles: four18/974original
endpoints, four24/2030endpoints, four18/974rigid endpoints, two18/974zero-source
states and two18/974repeats for2FVY. Source files remain byte-identical; every
new state gets a fresh native Model. No PCM cache substitution, scientific
retry, charge/geometry/water change, new DFT/MACE call or biological scoring.

One64CPU/128GiB allocation, groups serial. Record actual setup/source/solve,
CPU/memory and allocation costs. This higher-resolution inventory is one-time
development validation; ordinary primary-pair cost must be reported separately.
No project time/compute budget. Failure to converge or fit memory remains
explicit; do not fall back to a lower resolution or change the physical model.

Keep all v1 gates exactly: new primary-to-refined endpoint and contrast changes
<=0.1kcal; rigid changes<=0.05; reciprocity<=0.05; source potential<=1e-8au,
psi<=1e-12; contraction/zero/repeat<=1e-8kcal. Include between-structure checks.
The first inventory remains a failed frozen experiment even if this version
passes. No baseline reference, calibrated threshold, full-environment energy
or predictive-usefulness claim follows merely from numerical acceptance.
