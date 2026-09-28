# Frozen-f embedded preparer delivered

A four-cell Hans8DQ2 manifest passes the existing dry-run workflow. Six real-artifact tests pass in1.77s, including both remaining prepared sources, sulfur coverage, Mex charge0, physical versus effective spin, and rejection of corrupted real-state/basis inputs. No ORCA call, scheduler submission, new coordinates or state preparation occurred.

The new `scripts/metal_environment_lady_embedded_prepare.py` supplies only `prepare` and `dry-run`. It reuses existing source validation, snapshots the existing native executor/renderer, and pins the literal target basis/ECP/AuxJ assets. No shared script, frozen-f adapter, production path or prior input was modified. Tests: `tests/test_metal_environment_lady_embedded_prepare.py`.

Manifest: `workspaces/metal_environment_response_20260926/lady_frozen_embedded_prepared_v1/Hans8DQ2/manifest.json`, SHA256 `c8fc40f9c46e47235597d570c64212dcb60b5b9c923d4a568a3f309cdceca98e`. Exactly La_A,Dy_A,La_B,Dy_B;195 atoms/cell; all same external field. Charge−1, backend effective multiplicity1; physical La1/Dy6. Both target representations have806 explicit valence electrons; prepared Hans8FNR/Mex8FNS have816/832 respectively under their actual charges−1/0. No remaining-source endpoint manifests were emitted.

The energy is the embedded electronic component only. `DoEQ false` excludes external–external Coulomb energy (the existing project convention), not electron–external or nuclear–external interactions. It is not a complete hybrid energy or total physical force. The prior core/cap mapping and finite-field exclusions remain exact. Root must integrate a new embedded parser that verifies actual ECP/electron/state/components and records `.pcgrad`; the old isolated parser must not be silently reused as a claim of embedded force qualification. Root also owns the execution gate after active derivative checks.

Four56-rank worker metadata is a proposal, not an allocation or cost measurement. The experiment was prepared before results, with no biological sign target. No classification or affinity is returned. Numerical and physical qualification remain unavailable until actual execution and appropriate checks.

Reproduction from repository root uses an unused output directory:

```bash
python scripts/metal_environment_lady_embedded_prepare.py prepare --inputs workspaces/metal_environment_response_20260926/lanm_ef3_preparation_v1/Hans8DQ2/INPUTS.json --plan diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/EMBEDDED_PREPARATION_PLAN.md --output workspaces/metal_environment_response_20260926/lady_frozen_embedded_reproduction/Hans8DQ2 --la-basis legacy/qmmm_lc/ecp_lib/orca_La.ecp_basis --dy-basis benchmarks/hans_lanm_dy_qmmm_correction_v1/qmmm_calibration/basis/lcecp1_tzvp_Dy.orca.inc --la-aux diagnostics/metal_environment_response_20260926/electronic_diagnosis_20260928/lcecp1_tzvp_La_autoauxj_orca611.inc --dy-aux benchmarks/hans_lanm_dy_qmmm_correction_v1/qmmm_calibration/basis/lcecp1_tzvp_Dy_autoauxj_orca611.inc --mpi-ranks 56 --workers 4
python scripts/metal_environment_lady_embedded_prepare.py dry-run --manifest workspaces/metal_environment_response_20260926/lady_frozen_embedded_prepared_v1/Hans8DQ2/manifest.json
python -m pytest -q tests/test_metal_environment_lady_embedded_prepare.py
```

Changing source, method, assets, fields, charge/state, mappings or implementation changes manifest/cache identity; valid baseline outputs cannot fill these new tasks.
