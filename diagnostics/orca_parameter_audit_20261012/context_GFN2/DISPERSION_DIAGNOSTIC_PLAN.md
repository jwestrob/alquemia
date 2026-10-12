# Narrow dispersion diagnosis, declared before extra evaluations

The four-cell primary Linux xTB scout reproduces the archived corrected Mac
energies within 1.1e-12 Hartree. All tblite charge distributions agree within
8.4e-8 e, but absolute energies differ by 6.47–6.73e-5 Hartree. The original
1e-6-Hartree parity gate remains failed; the remaining 76 matrix points are held.

Inspection of xTB6.7.1 `src/scf_module.F90` shows the non-self-consistent
three-body D4 contribution uses a 40-Bohr cutoff. tblite0.7.0
`src/tblite/disp/d4.f90` uses 25 Bohr. Pairwise cutoffs also differ (xTB's SCF
matrix uses 72 Bohr, tblite uses 50 Bohr). This is a concrete possible source of
the non-electronic residual, not an established atomic reference shift.

Execute exactly two real, already consumed scout geometries (Ca vacuum and
La vacuum) under each backend with **only D4 s9 changed from 5 to 0** in separate
exported parameter files/dictionaries. Four new molecular calls maximum, no
optimization, no coordinate/state changes, no retries. Preserve all original
parameters and outputs. This is a causal component diagnostic, not a proposed
scorer or a revision of the independent-check gate. All other settings and
two-total-thread restriction match PLAN.md.

If deleting this term removes the discrepancy to the original 1e-6-Hartree gate,
the discrepancy lies in the non-self-consistent three-body term to that precision;
the distinct cutoffs are then a source-supported explanation requiring separate
verification before a claimed exact correction. If not, retain the residual and
do not expand the matrix. Geometry-dependent dispersion must not be called an
element-specific constant or subtracted empirically to force agreement.

Source links:
- https://github.com/grimme-lab/xtb/blob/v6.7.1/src/scf_module.F90
- https://github.com/grimme-lab/xtb/blob/v6.7.1/src/disp/dftd4.F90
- https://github.com/tblite/tblite/blob/v0.7.0/src/tblite/disp/d4.f90

Runtime and source hashes, exact inputs, modified parameters, raw output and
actual timing are retained in `workspaces/orca_parameter_audit_20261012/context_gfn2_diagnosis_v1/`.
