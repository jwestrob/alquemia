# Finite electrostatic input derivatives — 26 September 2026

The reusable potential/field chain rule passes five algebra tests on archived
1H4I coordinates and charges. This qualifies the implementation's algebra,
**not MACEPOL-EF, an embedded Hamiltonian, physical boundary forces, or prediction**.
No molecular executable, GPU job, installation, or new scientific endpoint ran.

## Energy and force convention

`scripts/metal_environment_fields.py` uses atomic units throughout: positions
bohr, charge e, potential Eh/e, electric field Eh/(e bohr), force Eh/bohr.
For `d = x_QM - x_MM`, `phi = sum(q/r)` and `F = sum(q*d/r^3)`, zero at infinity.
For supplied derivatives `a = dE/dphi`, `b = dE/dF`, each pair contributes

```
dE/dx_QM = q * [-a*d/r^3 + b/r^3 - 3*d*(b.d)/r^5]
dE/dx_MM = -dE/dx_QM
```

The code returns the negative coordinate derivatives. Add intrinsic electronic
forces at fixed phi/F and separately accounted MM/cross interactions. It does
not add another Coulomb energy or assume learned charges/dipoles are energy
derivatives. MM charges are fixed, not responsive. Covalent cap projection must
use the actual source Jacobian; no unsupported boundary mapping is supplied.

## Real fixture and tests

Pinned core/environment PQR files and hashes are in `RESULT.json`. The archived
PB environment includes zero-charge cavity centers coinciding with all 47 core
atoms. Tests exclude exactly zero-charge rows (49 total), retaining 9,094 actual
nonzero MM charges; charged overlaps remain explicit errors. This fixture
adaptation is not a proposed preparation policy.

Both derivative channels are exercised with the deliberately mathematical
probe `q_core.phi + 0.5*sum(F^2)`. This probe is **not a physical energy model**
and its values are not reported as molecular evidence. Directional differences
use `h=1e-4 bohr`, tolerance `2e-8 Eh/bohr`, fixed before residual inspection.
QM and MM directional residuals are respectively `1.19e-12` and `1.18e-12
Eh/bohr`; total-force residual is `9.03e-17 Eh/bohr`. Joint rigid transformation
changes the probe by `1.02e-15` and transformed forces by at most `1.58e-16`.

Additional checks cover unchanged repeatability, chunk-size independence,
`Q*delta_phi` gauge accounting, corrupt real-fixture overlap rejection, and
fingerprint invalidation. The cache helper requires nonempty state, geometry,
environment, model, and protocol records. Actual adapters must include all
scientific settings and externally updated fields; hashing cannot discover
missing caller metadata or repair a calculator's internal stale cache.

Run:

```bash
python -m unittest discover -s tests -p test_metal_environment_fields.py -v
```

Five tests passed, zero skips in 2.06 seconds (test-run duration, not molecular
timing or allocated compute). Missing real fixtures explicitly skip; altered
fixture hashes fail. Actual candidate/interface checks remain outstanding,
including charge-inclusive gauge response, boundary projection, intrinsic plus
environmental force consistency, and electronic-state support.
