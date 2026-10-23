# Hydrated panel handoff

Physical preparation:18 new configurations and2 original aquo8 controls.
Workspace prepared_v2/PANEL.json hash1191e3a00e3c02aa861faf83aae22eb250c94064e7ac780bf2fbebb00eb6c21a.
All geometric/stoichiometric gates pass; no energies were used in construction.
The first prepared/ output is preserved; v2 adds the correctly matched free
water bookkeeping geometry. It does not change the18 complex coordinates.

Root-ready binding is bound_v1/ under the matching workspace:

- NATIVE_PANEL.json a551d1187b73851c349f8bffc19e936a0bac561ca8538d7a06adf2a1b088c7c9:
  exactly6 central ECP PBE0 cells with unchanged REFERENCE.py and basis assets.
- MACE_CONFIG.json7e236ff044ea89cf1367b4c5dea9a782051a1d7f33c5f2fbb397c7a0e6909815:
  exact18 fixtures with root's unchanged RUN.py and original native model config.
- NATIVE_VALIDATION.json verifies all6 paired geometries and rendered decks.
- READY.json records source/input pins, command templates and preparation cost.

Actual preparation commands (one Python process; OMP, MKL, OpenBLAS, NUMEXPR
and BLIS thread controls all1):

```sh
python3 diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/PREPARE.py --out workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/prepared_v2
python3 diagnostics/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/BIND.py --out workspaces/lanm_selectivity_free_energy_20261021/reference_chemistry_v1/hydrated_mixed_v1/bound_v1
```

Internal elapsed times0.256s/0.508s, zero molecular calls. Root owns submission.
RUN_SPEC.md defines scoring and collection; the initial56-cell collector is
inapplicable. Reuse original aquo8 outputs only after their admission. Report
summed ligand generalized force for the actual rigid-ligand distortion; the
frozen native metal-projection field is a different quantity. No AE extension,
parameter fit, sampling or optimization is included.
