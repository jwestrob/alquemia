# Source audit of cross-code discrepancies

This records source inspection after the actual La-aquo baseline 22327 failed
the fixed 1e-4 Eh archive comparison. Its ordinary SCF converged and passed
electron-count, occupied-orthonormality, density-idempotency and orbital-gradient
checks. Electronic difference: +0.000825580291 Eh / 0.518059454 kcal/mol.
No archived-reference or implementation-validation gate is relaxed here.

## ORCA versus PySCF

- The archived La ECP has a zero local g channel and s/p/d/f nonlocal projectors.
  The loader retains this operator, including radial powers. PySCF's named
  default La ECP is different and was not substituted. The independent NWChem
  parser agrees with the custom input conversion.
- Supplied H/C/O primitives and contracted shells agree with PySCF def2-TZVP.
  Positive amplitudes differing on single-primitive p shells disappear under
  contraction normalization. ORCA reports the same 391 AOs, 167 shells and 88
  explicit electrons. Its internal light-element primitive tables were not
  printed, so exact internal-library equivalence remains a provenance limit.
- ORCA used RI Coulomb with 919 auxiliary functions and COSX exchange; PySCF
  evaluates direct J/K. The last ORCA COSX correction is +0.000045391 Eh. That
  is smaller than the discrepancy, but it does not bound total COSX error.
- ORCA's DFT grid has 275,796 points and integrates 88.000143713891 electrons.
  PySCF uses 3,652,504 unpruned points. This establishes a quadrature difference,
  not a quantitative energy-error estimate.
- ORCA documents PBE mu=0.21952 and beta=0.066725; LibXC7 uses
  0.2195149727645171 and 0.06672455060314922. A fixed-density measurement is
  declared in PLAN.md and submitted as 22416, xc_fixed_density_diagnosis_v1.
  It does not modify a default or claim complete functional equivalence.
- ORCA's documented Bohr conversion is 0.5291772083 Angstrom, versus PySCF's
  0.52917721092. The observed nuclear-repulsion difference is about 3.225e-6 Eh,
  much smaller than the total discrepancy. This is excluded from the XC test.

No concrete translation defect was found. Coulomb fitting, quadrature and
functional constants remain separate possible contributors; their combined
effect has not been measured.

Sources: archived native output at
`workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/ecp_remaining_v1/tasks/ecp_pbe0__aquo8_central_La/execution/endpoint.out`,
especially lines 678–681, 723–740, 823–840, 987–991 and 1016–1018;
[ORCA functional defaults](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/DensityFunctionalTheory.html#customization-of-basic-functionals),
[ORCA unit conversion](https://www.faccts.de/docs/orca/6.1/manual/contents/orcaarchitecture/conversionfactors.html),
and the pinned PySCF/LibXC runtime and input assets.

## CP2K water benchmark: source inspection while jobs run

The developing native water monomer energy differs substantially from PySCF.
Terminal component comparisons, not that intermediate absolute offset alone,
will determine whether this benchmark passes. Both original jobs are preserved.

Literal coordinates, fragment membership, electron counts, spherical AO counts,
basis coefficients, ALL nuclear potentials and PBE selection agree. CP2K uses
primitive plus contraction normalization. Its isolated-fragment calculation
retains the parent's DFT settings and original coordinates. The offset therefore
appears before constrained ALMO optimization.

Coordinates are near a cell face and automatic centering is disabled. Source
inspection found modulo grid deposition and a nearest-image MT kernel, so this
alone does not prove a boundary error. The documented MT requirement concerns
cell size relative to the entire density distribution.

The decks use the standard GAPW one-center XC path, with
GAPW_ACCURATE_XCINT defaulting to false, GAPW_1C_BASIS=ORB, EPSFIT=1e-6 and
fixed atomic grids. The 400/600 Ry comparison tests plane-wave resolution only;
it cannot establish convergence of every GAPW approximation. No new parameter
scan or changed settings is authorized by this source observation.

CP2K's native PBE constants also differ slightly from LibXC7: beta=0.066725,
mu=beta*pi^2/3 and PW-LDA A=0.031091 in `xc/xc_pbe.F`, versus LibXC's
higher-precision beta/mu and modified-PW A. This is nominal PBE matching, not
bit-identical functional parameterization. Its numerical effect is unmeasured.

Relevant frozen CP2K source: `mscfg_methods.F`, `cp_subsys_methods.F`,
`input_cp2k_qs.F`, `qs_vxc_atom.F`, `grid/grid_api.F`, and `pw/mt_util.F` in
the pinned CP2K 2026.2 source tree recorded in tests/cp2k/RUNTIME.json.
[MT documentation](https://manual.cp2k.org/cp2k-2026_2-branch/CP2K_INPUT/FORCE_EVAL/DFT/POISSON.html).

## Measured fixed-density XC result (22429)

The explicit LibXC PBE0 expansion reproduces native PBE0 energy/Fock to
4.55e-13 Eh / 1.78e-13. Substituting only ORCA documented mu/kappa/beta
changes the same-density energy by -8.503404706e-5 Eh. The archive residual
remains +0.000740546245 Eh: most of the discrepancy is not explained.
No SCF, baseline admission, fitted constants or gate change. Exact results:
workspaces/lady_almo_implementation_20261023/xc_fixed_density_diagnosis_v3/result/RESULT.json.
Two previous metadata/dispatch failures and their saved evaluations remain.

## Native saved-orbital J/K result (22436/22437; collection22441)

Both native one-Fock calculations completed without orbital optimization. The
RIJCOSX replay differs from the archived electronic energy by -0.000045388230
Eh, within the original1e-4Eh gate. Direct J/K raises energy by0.000735320167Eh
(0.461420371kcal/mol) relative to that replay; D4 is identical. This is a
substantial numerical-backend contribution to the cross-code discrepancy.
The replay differs from the archived post-final-COSX energy, and the PySCF
density is separately optimized; these differences must not be silently
combined into an exact closure or used to admit the failed baselines.

A parser incorrectly treated startup auxiliary V-metric preparation as proof
of active RI-J. Saved-only correction22441 requires actual native J timing in
the direct run and SPLITRIJ/COSX timing in the control, unchanged explicit
inputs/Hamiltonian flags and every other original check. Both original parser
failures and two earlier MPI startup failures are retained. No native rerun
was needed for the parser correction. See
workspaces/lady_almo_implementation_20261023/orca_backend_saved_recovery_v1/result/REPORT.md.
