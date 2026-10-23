# Fixed explicit polarizable panel: complete, transfer unqualified

The native Tinker panel completes **28 geometries plus four isolated ligand fragments**, with all four declared water/bidentate-acetate derivative checks passing. This establishes an executable candidate and narrow numerical consistency. The asymmetric La/Dy parameter provenance remains a material physical limitation.

## Actual execution and numerical checks

The serial local controls reused the original La-water output and made six new native program invocations for three missing structures: 3.44 seconds elapsed. Maximum analytic-versus-numeric gradient component errors were 1e-8 (La water), 2e-8 (Dy water), 7e-8 (La bidentate acetate), and 1.1e-7 (Dy bidentate acetate) kcal/mol/Å. All pass the predeclared 1e-3 threshold at 1e-5 Å. Values are derived from the native eight-decimal gradient tables. No tolerance was changed.

Root-submitted job **17272**, one CPU, reused those four controls and completed the other 28 entries in 56 native program invocations; adapter elapsed time was 25.03 seconds. Every entry has finite total/component energies and full analytical gradients. The finite-difference executable internally evaluates displaced energies, so a program-invocation count is not an energy-evaluation count. No GPU, ALPB, geometry optimization, dynamics or parameter fitting occurred.

## Central fixed-geometry interaction energies

These are **complex minus isolated ligand** energies at the same ligand geometry, in kcal/mol. The spherical isolated ion has zero self energy in this classical model. Do not substitute Tinker's printed intermolecular component for this subtraction: induced polarization also changes intramolecular energy.

| Fixture | La | Dy |
| --- | ---: | ---: |
| Water | −91.27177843 | −98.76432953 |
| Acetamide | −182.51828142 | −158.47018439 |
| Bidentate acetate | −507.01945884 | −545.65259699 |
| Monodentate acetate | −417.57024714 | −448.27534039 |

The La results reproduce the published SI values to its quoted precision or about 0.01 kcal/mol. The SI's known monodentate-acetate weakness remains: its QM value is −466.67 kcal/mol versus AMOEBA22 −417.58. Numerical consistency does not repair that transfer limitation. The eight- and nine-water entries are fixed clusters, not sampled coordination populations or hydration free energies.

## Acetamide mapping review

The unexpected weaker Dy acetamide interaction is present in the declared model; the saved inputs show no atom-type or coordinate mismatch. La/Dy inputs differ only in metal element/type 500/501. Ligand types remain `105,106,109,109,107,108,109,107,104`; their coordinates and bonds are unchanged. Dy501 exactly preserves the author-contributed Dy10 charge +3, zero permanent higher multipoles, polarizability 0.7280 Å³, Thole 0.1390, and buffered-14-7 parameters 3.2300 Å/15.0000 kcal/mol.

La's published amide-specific `POLPAIR 105 500` and `VDWPR` terms correctly apply only to La. There is no invented Dy donor correction. At the central acetamide geometry, the permanent-multipole energy is identical (−122.82208715 kcal/mol); the La/Dy polarization terms are −102.64772661/−82.07695172 and vdW terms are 29.81939197/33.29671411. This identifies the implemented terms responsible for the difference without claiming a unique physical energy decomposition or validated preference.

The Dy ion originated with different GEM-water multipoles. Its common-water and donor transfer needs independent reference qualification. Preserve all donor classes and matched metals when assessing that gate; a water-only pass cannot select this Hamiltonian for protein sampling.

## Protein implementation boundary

The La SI parameter file is the small-molecule set. The native `amoebabio09.prm` supplies protein and amidated-terminal biotypes and compatible water parameters, but uses different atom numbers and donor multipoles/polarizabilities. Copying the SI key into a protein model would misassign corrections. Geometry/typing preparation may proceed, with explicit source joins, while donor correction transfer and a complete protein Hamiltonian remain unqualified.

## Evidence

- Native panel: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/panel_v1/RESULT.json`, SHA256 `b256fdca923c0f476651e73250b93672eb49f7634fbceca66d2e6ea1d398a052`.
- Controls: `controls_v1/RESULT.json`, SHA256 `f2d70fc06fad5ae19566c9fd779b14cbcf48f170614ccae24320a482c4c0deaf`.
- Frozen adapter `RUN_PANEL.py`: `4189602145853bc39d8dfaf30915b2a6b1d4632984e0c50fc0bf321644fa306c`; panel config `33cfc47dff736afcae68bd17ff4f19757993527e5bf384b4c732289bb480d82b`.
- Every row pins the exact prepared structure/key/parameters, raw native outputs and actual subprocess receipts. Original build and preparation failures remain preserved.
- Primary [La publication and SI](https://doi.org/10.1021/acs.jpcb.2c07237); [pinned author-contributed Dy parameter file](https://github.com/TinkerTools/tinker/blob/87050685eff8840d312e2a332cc82c33f63c7c3d/params/amoebail22.prm).
