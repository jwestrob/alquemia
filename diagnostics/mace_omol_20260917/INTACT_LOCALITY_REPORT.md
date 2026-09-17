# Intact OMOL: local context survives, distant contributions cancel

All20 archived primary endpoint readouts were rechecked against their actual
execution receipts, coordinates and charges. All18 component/score/contrast
closure checks pass; maximum discrepancy2.575e-8kcal/mol. No new model forward,
geometry, force or score was calculated. Local execution cost7.99wall seconds,
10.31CPU seconds,462348KiB peakRSS.

The geometry-independent embedding and atomic reference terms cancel exactly
within each bound-minus-detached pair. For atoms beyond18A, the sum of absolute
node changes in any individual endpoint transfer is at most3.66e-11kcal/mol.
This is not merely cancellation of large positive and negative distant terms.
It is consistent with the pinned6A graph cutoff and three local interaction
blocks. Model-wide total-charge conditioning still exists in the nonlinear
features; this audit has not removed or validated that dependence.

Native learned score contributions, kcal/mol:

| Case | Metal node | Other atoms within6A | 6–12A | 12–18A |
|---|---:|---:|---:|---:|
| GGR1GLG | 63.046 | 2.143 | -2.572 | 0.099 |
| alpha1F6S | 73.461 | 20.509 | 5.640 | 0.259 |
| alpha6IP9 | 71.999 | 13.582 | 5.050 | 0.173 |
| MxaF1H4I | 78.367 | 16.468 | 7.134 | -0.356 |
| XoxF4MAE | 81.061 | 31.585 | -2.183 | 0.185 |

These are the network's native energy bookkeeping terms, not uniquely defined
physical atomic energies. In particular, the metal node already contains
information from its neighbors; its term is not a bare-ion potential.
Source-resolved water terms are retained, not omitted to improve a direction.
All sums reproduce the published five-case results and three contrasts.

**What transfers:** intact nearby chemistry and matched subtraction support a
useful finite-range candidate without chemical caps. The result does not
demonstrate full long-range electrostatics or isolate the cause of its accuracy
gain over cores. Total charge changed between core and whole representations,
as did the local chemical environment. A separately declared disconnected
charged-spectator test now examines that conditioning ambiguity. No new
classification bands, truncation policy or production change follows here.

Full source-mapped terms: `workspaces/mace_omol_20260917/intact_locality_v1/result.json`.
Frozen implementation: `intact_locality_source_v1/implementation/` in the same
workspace tree. Repeat into a fresh output with:

```bash
workspaces/mace_hybrid_20260916/software_v1/venv/bin/python workspaces/mace_omol_20260917/intact_locality_source_v1/implementation/mace_omol_locality.py --source-report workspaces/mace_omol_20260917/intact_report_v1/result.json --agreement diagnostics/mace_omol_20260917/INTACT_LOCALITY_PLAN.md --output workspaces/mace_omol_20260917/intact_locality_replay_v1
```
