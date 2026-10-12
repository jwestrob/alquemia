# Independent corrected expanded canonical review

12 October 2026. **The archived expanded-context candidate retains a valid
canonical calibration after the La parameter correction. Its accommodation
narrows the training separation relative to its own origins.** The frozen
operational calibration gate passes, permitting the already-declared structural
transfer evaluation. This is not a predictive accuracy result or a promotion.

## Independently recovered result

The review reads original xTB stdout and charge files, tblite scientific output,
execution receipts, archived MACE records, and source state signatures. It does
not import the scoring runner. All 25 designated canonical members are present:
14 Ca-class and 11 La-class biological groups. Both metals receive each source's
same archived origin/Ca-proposal/La-proposal pool.

| Score policy | Canonical gap, kcal/mol | Ca boundary | La boundary |
|---|---:|---|---|
| Operational accommodation | 3.472171047 | Q9Z4J7 | C5AXV8 |
| Mathematical row minima | 3.472171047 | Q9Z4J7 | C5AXV8 |
| Origin only | 10.297678912 | Q9Z4J7 | A0A3F2YLY8 |

The rule is unchanged: all 25 members must be available, the gap must exceed
0.02 kcal/mol, and bands are max(Ca) and min(La). Origin retention still requires
more than 0.1 kcal/mol improvement before selecting another candidate.
The origin full-panel gap differs from the earlier two-anchor origin ordering
because a different La reference defines the boundary. A larger training gap
does not establish better transfer; that comparison remains the next test.

Every independently recomputed selected geometry and score agrees with the
collection (score comparison tolerance 1e-7 kcal/mol). Exact bands and individual
works are retained in `RESULT.json`. Scores are the MACE-vacuum plus matched
GFN2 ALPB-minus-vacuum Ca-minus-La contrasts on this protocol's scale, not
binding free energies. No old native-parameter bands are imported.

## Numerical qualification and preserved failures

All 300 input records, comprising 600 backend evaluations including reused
work, were checked against raw results and output hashes. All 150 solvent
components pass the previously declared v2 gate. Maximum backend component
disagreement is 1.298692496e-7 Eh against a 1e-6 Eh limit; maximum all-atom
population discrepancy is 1.093457009e-7 e against a 1e-5 e limit. Convergence,
finite values, electron/charge closure, source charge/multiplicity and one-thread
environment checks pass.

**All 300 absolute backend-energy comparisons still fail the original v1 gate.**
Their maximum residual is 1.938589135e-4 Eh. The already-documented differing
non-SCF dispersion conventions cancel within each matched solvent subtraction;
absolute total-energy parity is neither claimed nor relabeled as passed.

## Actual incremental cost

| Backend | New calls | Sum of child CPU seconds | Sum of call wall seconds |
|---|---:|---:|---:|
| xTB | 276 | 1979.057753 | 1988.662059 |
| tblite | 276 | 8106.917830 | 8120.915695 |
| Total | 552 | 10085.975583 | 10109.577754 |

The 14-worker numerical stage lasted 744.688672 seconds (12.41 minutes),
excluding preparation and subsequent collection. Maximum recorded single-child
RSS was 127052 KiB; this is not an aggregate process-tree peak. No GPU or Slurm
allocation was used. An additional 24 input pairs (48 backend calls, 12 solvent
components) were reused from the completed factorial, incurring no new molecular
calls; their original costs remain separately recorded. This audit ran no model.

## Reproduce

From the repository root:

```bash
python diagnostics/orca_parameter_audit_20261012/expanded_canonical_review/REVIEW.py
```

The audited collection SHA-256 is
`929dedab7ba46c7088de7a4a695a8f19fc5af60066ecd3599645f687bd4e6896`.
The finite transfer may proceed under the frozen calibration; production,
historical failures, and the separate common-five candidate sensitivity remain
unchanged.
