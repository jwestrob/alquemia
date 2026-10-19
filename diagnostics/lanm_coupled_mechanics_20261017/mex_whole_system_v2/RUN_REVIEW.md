# Native panel review: 14790 and 14792

**PASS for saved collection.** All eight declared whole2049 native MACE cells
are complete. Job14790 made one scout forward; job14792 reused its exact raw
and assembled records and made the remaining seven forwards. No failed native
request, additional state, solvent/MM/GFN/DFT call or optimization is recorded.

## Immutable results

- Scout: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_scout_v1/execution/RESULT.json`,
  SHA256 `5eb4ca1738a19098de9d4788d812aa31ac610b4b32aedcbcc022f62475a6dd8c`.
- Remaining: `workspaces/lanm_coupled_mechanics_20261017/mex_whole_remaining_v1/execution/RESULT.json`,
  SHA256 `4df04ed07607c7340a3cac1b7e6c17db87f74e2de8bc943e926f50996fa43fe5`.
- Both replay the admitted preflight manifest
  `acdaa4623a85a6e36fa7edc23d43a94fd94bc51ff9465aea69834a41f1b8d8d0`
  with unchanged source `dd611ca0` and checkpoint `9b64b4fd`.

The eight ordered tasks match the frozen two geometries × two metals × two
spin hypotheses. Actual batch receipts report2,049 atoms, one graph, charge+2,
target index2045 and the correct multiplicity10/15 or4/9. Before-state target
edge counts are206 and after-state counts198. No charge prediction is claimed.
Both jobs used NVIDIA H200, float64, the unchanged native forward, no execution
adapter or extra dispersion, and torch2.8.0/mace-torch0.3.16/ASE3.29.0.

Every raw-force and assembled-artifact pin was verified. All arrays retain
2,049 physical rows; native energy metadata, START/task bindings, source IDs,
four reused ALPB receipts and their saved component-gradient arrays match.
The first remaining row equals the original scout row plus its explicit reuse
pin. There is no repeated scout calculation.

## Recorded resources and continuation gate

| Allocation | New forwards | Model load, s | Application elapsed, s | Process CPU, s | Peak parent RSS, KiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 14790 | 1 | 6.613859 | 9.045296 | 7.080429 | 1,550,432 |
| 14792 | 7 | 6.720614 | 16.140243 | 14.838763 | 1,871,224 |

Each allocation reports12 CPUs. The scout forward took1.076935s, with peak
allocated device memory81,684,463,616 of150,109,880,320 bytes. Its recorded
600-second/90%-memory continuation gate passed. The remaining job reports the
same peak allocated memory. These are actual application measurements;
scheduler allocation costs are recorded separately by root.

This review inspected saved metadata and artifacts without recomputing energies,
force transformations, work or norms. The saved collector must still perform
its exact raw/component/total replay and declared comparisons. Execution
success establishes neither whole-protein force accuracy nor a physical local
spin arrangement. Original preflight14788 and the closed hybrid/reference
failures remain preserved.
