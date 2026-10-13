# PedE saved-energy diagnosis

Question: where do the model's source-dependent origin contrasts and accommodation
works arise—OMOL vacuum energy, the matched GFN2 solvent contribution, or both?
This decomposition does not determine whether a fold or the energy model is wrong.
Root's separate fixed-band analysis handles classification margins and calibration
scale changes; this branch does not attribute calls from raw shifts alone.

## Frozen finite scope

- Every recovered PedE/Q88JH5 two-Ca sample 0–4: all six metal/geometry cells.
- Exact corrected-method La-conditioned samples 0, 2, 3 and 4, including the
  already-designated canonical sample 2. Sample 1 is absent from the pinned
  corrected expanded manifest and is explicitly unavailable for this comparison.
- No new electronic evaluations, optimizations, DFT, threshold fitting, or voting.
  These are consumed correlated development/calibration sources, not blind tests.

The analysis requires identical checkpoint/model, scalar solver settings/runtime,
and the same local chemical state signature across all nine sources. It verifies
charge/multiplicity and paired coordinates from the saved finite pools. Spectator
ions are retained in full sources but absent from these local energies; this is
not a global occupancy-response calculation.

## Algebra declared before execution

Use the exact existing conversion constants and component-v2 energy records:

`E_M(q) = V_M(q) + U_M(q)`

`V = E_OMOL,eV × EV_TO_KCAL`

`U = (E_GFN2,ALPB − E_GFN2,vacuum) × HA_TO_KCAL`

For each metal and every candidate, retain `W(q) = E(q) − E(origin)` and its
separate V/U differences. Use origin-relative differences before conversion to
avoid subtracting large converted totals. Retain both own and cross candidates.
Keep saved operational and mathematical selections; verify their algebra under
the unchanged 0.1 kcal/mol origin-retention policy. A vacuum-minimum proposal need
not be the composite energy minimum, and neither is an unconstrained composite
stationary point.

For each policy, report `R = E_Ca − E_La = R_V + R_U`. Compare origin components
relative to frozen canonical PedE sample 2. Within each declared source cohort,
decompose the origin contrast range using the *same sources* defining the minimum
and maximum composite R, so component differences add to the actual range.
No source is selected or excluded from the scoring population by these extrema.

Serialization/algebra check tolerance: 1e-6 kcal/mol, independent of classification.
Paired coordinate check: 1e-12 Å. Preserve the existing component-v2 residual and
population gates; report saved residuals rather than rerunning their raw parser
audits. The historical absolute-v1 backend parity failure remains failed.

## Execution and output

Only Slurm executes `ANALYZE.py`; it requires an actual job ID and scheduler record
before scientific algebra. No molecular software is imported. A one-CPU job reads
the pinned records and writes `RESULT.json`, `CELLS.tsv`, `WORKS.tsv`, `REPORT.md`
and its actual execution receipt. Existing output directories are never reused.

All 54 cells are retained. Missing or inconsistent required artifacts fail the
analysis visibly; no energies are replaced by zero and no source is omitted.
Root owns submission and interpretation. Syntax/JSON/hash inspection only occurs
on login. The final report and vault note follow the actual Slurm result.
