# Separate angular-basis error from surface-integration error

The second conductor inventory completed all roles and passed 444/450 checks.
Its GGR between-structure source-self contrast changes only 0.01879 kcal/mol
between 18/974 and 24/2030. However, all four endpoint energies still change
0.265–0.417 kcal/mol, and two rigid-transform endpoint changes slightly exceed
0.05 kcal/mol. The frozen full qualification remains failed.

The previous tests changed angular basis and quadrature together. Before
expanding either further across proteins, separate these numerical effects
on the already consumed 2FW0 Ca/La pair. This structure is selected for the
largest endpoint refinement error, not for its biological classification.

Execute eight new native conductor solves, two endpoints at each setting:

| Setting | lmax | Lebedev points | Comparison |
|---|---:|---:|---|
| integration18 | 18 | 2030 | archived 18/974: integration only |
| integration24 | 24 | 3470 | archived 24/2030: integration only |
| basis30 | 30 | 3470 | new 24/3470: angular basis only |
| integration30 | 30 | 5810 | new 30/3470: integration only |

The source/cavity, physical model, energy prefactor, native stock executable,
exact Coulomb potential, FMM12/12, eta=.1, shift=0, epsilon=78.3, solver
tolerance=1e-10 and maxiter=1200 remain identical. Fresh native Model per state.
No charge, radius, geometry or label adjustment. Record input and native
parameter hashes, source checks, raw/scaled energies and actual costs.

Use the existing 0.1 kcal/mol endpoint and contrast refinement screen for
each comparison; retain 0.05 kcal/mol reciprocity and the original source,
contraction/passivity checks. This diagnostic does not replace the complete
rotation/repeat/two-structure qualification. It may identify what needs finer
resolution; it cannot itself qualify a full score or reset previous failures.

One 64-CPU/128-GiB allocation, serial groups. The measured 24/2030 solves take
about 163 seconds per endpoint; larger grids/bases will cost more and are
one-time numerical development work. Record actual setup/solve/report and
allocation costs. No project compute/time budget, automatic retries, DFT,
MACE, force, response iteration or biological scoring calls.
