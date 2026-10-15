# Canonical geometry rebase: actual result

Job **5016** completed with **25/25 prepared sources**, zero energy evaluations and zero optimizations. Every source retained its original four angular modes and bounds, atom/state/cap inventory and complete source-contact mapping. Terminal amide χ2 remains outside the selected modes where originally absent.

The exact A0AC La-conditioned sample 3 regression passed. The new output aliases its already executed mapping, contacts, task origins and search results after comparing the actual rebased geometry, Jacobians, cap checks and full-source contact arrays. Its energies and searches must be reused. The remaining 24 sources have geometry-qualified q0 only; energy and projected-force qualification remain unrun.

The application receipt records **35.64470355899539 wall seconds** and **34.38098565200001 process CPU seconds**, on one allocated CPU. Root's terminal scheduler record is `COMPLETED 0:0`, 36 seconds. No GPU, SCC, DFT or model call occurred.

Actual evidence under `workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_relaxation_v1/`:

| Record | SHA256 |
|---|---|
| `manifest.json` | `4d7a83d44c0cf3b7631adfa334d273bc0db17389062b5908936e3a12b829415c` |
| `execution/PREPARATION.json` | `663ad699c2fae8baf6a42ee65078481a26c872d0fa6ad0d6cf17d3e8eeb0fbca` |
| `execution/EXECUTION.json` | `1f24f7a30d4f74c2a938e88a9871fcff1a5dd09b10d4ec3a5542093b08929330` |
| A0AC `A0AC_REGRESSION.json` | `041015ec60292f5ad27b6438c1ed2b4311bb0f4e2f60a5e7db590e700fb86b09` |

Admission establishes the declared local geometry, covalent, cap and severe-contact checks. It does **not** establish an equilibrium or clash-free conformer, validate every trial in the search domain, or imply an energy or classification gain.
