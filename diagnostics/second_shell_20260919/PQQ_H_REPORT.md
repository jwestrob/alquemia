# PQQ hydrogen preparation: lower energies, stable discrimination

Four native MACE proposals and four DFT single points completed. DFT still orders
XoxF above MxaF, with its gap changing **30.571889 → 30.071899 kcal/mol**. This does
not demonstrate improved accuracy or a larger DFT gap. Native MACE predicts a
slightly larger gap, 102.222180 → 104.027304, on its own vacuum energy scale.

| Case | DFT ΔR from H preparation | Native MACE ΔR |
|---|---:|---:|
| 1H4I / MxaF | −0.116776 | +0.841248 |
| 4MAE / XoxF | −0.616767 | +2.646372 |

`R = E_Ca − E_La`; units are kcal/mol. Individual DFT endpoint energies fall by
169.88–262.68 kcal/mol, almost entirely canceling between metals. Lower geometry
energies therefore supply no corresponding gain in discrimination here. They do
not invalidate the baseline's already observed PQQ signal.

## What changed

Only physical hydrogen coordinates move. Native OMOL proposes them in the
154/202-atom donor contexts, then transfers only hydrogens belonging to the
original 73/80-atom scoring cores. All original heavy atoms, artificial caps,
atom order, charge, cofactor state, protonation and water inventory remain exact.
Scoring retains the original r2SCAN-3c/CPCM(Water)/DefGrid3 single-point recipe.
No environmental energy is added to the original-core score.

All four endpoints meet the frozen 0.03 eV/Å criterion for constrained
stationarity and retain their covalent parents. Respectively 16, 17, 35 and 34
hydrogens reach the 0.35 Å displacement boundary. These are constrained
preparations, not unconstrained minima. 1H4I/La reaches SLSQP's 200-iteration limit
despite a projected maximum gradient of 2.34e−6 eV/Å. Its final iterate is retained
under the same declared proposal rule as the other endpoints. No criterion was
changed after inspecting results.

## Actual cost

- MACE job **1202091**: 264 GPU-seconds; 4,224 allocated core-seconds.
- Four searches: 2,113 recorded objective calls plus four original-core calls.
  ASE may cache repeated coordinates; this is not a count of independent forwards.
- Peak CUDA allocation: 6,122,767,872 bytes; peak process RSS: 1,847,408 KiB.
- DFT job **1202101**: 327 seconds on 64 CPUs, or 20,928 allocated core-seconds.
  All four endpoints terminate normally; no failures or retries.
- Combined: **25,152 allocated core-seconds and 264 GPU-seconds**. Preparation and
  local tests are additional, unmetered work. Earlier reference and static-context
  calculations are separate.

1H4I/La accounts for 184 of 257 worker-seconds. A future implementation could stop
at the already declared force criterion and avoid iterations spent satisfying a
much tighter energy criterion. This frozen pilot did not use that shortcut.

## Status

New protocols are `PQQ_native_OMOL_context_physical_H_preparation_v1` and
`native_r2scan3c_PQQ_context_prepared_physical_H_v1`. The baseline, default scorer
and canonical 25-case PQQ results remain unchanged. No bands are inherited and no
threshold is fitted. Both crystal controls were already consumed in development.

**Recommendation:** retain the current PQQ default. This pilot demonstrates lower
electronic energies and retained ordering; it supplies no accuracy-based reason
to promote the hydrogen preparation.

The exact result is
`workspaces/second_shell_20260919/pqq_h_dft_v1/collection_1202101.json`.
MACE proposals and trajectories are pinned in
`workspaces/second_shell_20260919/pqq_h_v1/collection_1202091.json`.
See the recorded agreement, submission files and actual-fixture tests alongside
this report.
