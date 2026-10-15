# Native geometry–charge factorial for the two Asn sample4 sources

Declared before new evaluations. Root owns Slurm execution under the existing
contained discriminator authorization. This is a diagnostic of the completed
Asn result, with no new scorer, class bands, fitted correction or production
change.

## Question and fixed inputs

Both Ca-conditioned4 and La-conditioned4 gain about9.5kcal/mol in the La solvent
term when their old `adaptive_La` geometry is replaced by `amide_flip`. Each
geometry had its own GFN2 vacuum charge vector. The observed change therefore
mixes geometric reaction-field response and geometry-induced vacuum charge
redistribution, despite being frozen against solvent-induced SCC response.

Use exactly these two sources, both exchanged metals, and both complete saved
204-atom geometries. Preserve the mapped atom order, physical identities,14 caps,
dry oxidized PQQ3−, neutral Asn, Ca−4/La−3 singlets and recorded636 GFN valence
electrons. No coordinates, protons, charges, radii or state assignments are
constructed or optimized. Short represented H contacts remain unchanged.

For each source/metal define:

- x0: actual `adaptive_La` coordinates; q0: that geometry's actual vacuum charges.
- x1: actual `amide_flip` coordinates; q1: that geometry's actual vacuum charges.
- Gij: unchanged native ALPB evaluated at geometry xi and charge vector qj.

The four source/metal groups each have four cells:16 total. **Reuse G00 and G11
in every group:8 diagonal cells. Evaluate only G01 and G10:8 new native calls.**
There are zero new diagonal calls, SCF/SCC calls, MACE calls, DFT calls, optimizers
or installations. The eight off-diagonal vectors are explicit counterfactual
geometry/charge combinations; they are not newly self-consistent electronic
states or proposed classifier endpoints.

## Qualification before any new call

Pin the actual completed4623 configuration, collection, execution and selected
eight input results, plus the existing native manifest, build, binary, runtime,
helpers and corrected-La parameter-lineage record. Reuse the maintained native
binary/factory with internal GFN2-xTB/ALPB water, P16 kernel, gsolv convention,
298.15K,230 angular points and unchanged radii. Do not recompile the engine.

The eight reused diagonal results are actual previously executed native results.
Before dispatch, verify their exact executable command/binary, input XYZ/charge
hashes, successful execution and raw stdout. Reparse all four printed native
terms and both rounding bounds and require exact equality to the saved records.
Verify component closure within the existing1e−9Eh allowance and charge sums
within the existing204×5e−9e rounding allowance plus1e−12e. The existing
component-v2, water-diagonal and vacuum-precision gates must still pass. This
qualifies genuine saved replays; it does not manufacture a zero-call fresh replay.

For every crossed pair, require the same physical maps, source-preparation pin,
state signature, ordered element identities and formal state. q remains at its
original mapped atomic indices. No interpolation, renormalization, atom exchange,
permutation or cap-charge omission occurs. The source of the geometry and source
of the charge vector are separately recorded in every off-diagonal receipt.

New evaluations call the unchanged `NATIVE.evaluate` interface directly, with
the saved geometry and the other geometry's vacuum-charge file. The new cavity
and Born-matrix hashes must equal the reused diagonal for that same geometry;
Gsasa and Gshift must match within the existing1e−9Eh allowance. All four terms
and the rounded charge sum are retained. Each term keeps the unchanged native
precision gate: rounding bound plus1e−9Eh, compared to the existing1e−6Eh
informative ceiling. No relaxed gate or automatic retry is introduced.

## Both paths and interaction

Report these quantities separately for Gelec, Ghb, Gsasa, Gshift and their sum,
in Eh and kcal/mol using the pinned627.509474 conversion:

| Quantity | Expression |
|---|---|
| Actual diagonal change | G11 − G00 |
| Geometry first, fixed q0 | G10 − G00 |
| Charge change at x1 | G11 − G10 |
| Charge first, fixed x0 | G01 − G00 |
| Geometry change at q1 | G11 − G01 |
| Mixed interaction | G11 − G10 − G01 + G00 |

Both paths must close to the diagonal change, and both equivalent interaction
identities must agree, within the existing1e−9Eh native allowance. Propagate each
stored-charge bound plus the existing allowance conservatively with absolute
coefficients, then sum term bounds for total energies. Preserve component and
aggregate informative flags against the same ceiling. An aggregate false flag
remains visible and is not silently promoted or used to discard an unfavorable
result. This propagation adds no invented charge precision or stochastic error
estimate.

“Geometry” includes changed interatomic distances, Born/reaction geometry and
geometric HB/SASA weights. It is not a pure cavity-volume intervention. A mixed
term means attribution depends on the path; neither path is privileged by its
sign. The adaptive_La geometry is compared under both metals to examine the same
two complete arrangements; it need not be the Ca winner of the original pool.
No full common-pool score or classification is inferred from this factorial.

## Execution and failure visibility

One shared standard-partition CPU; all numerical thread variables1. Serialize
the eight native calls. Existing GNU-time/native receipts record each actual
call and its elapsed/CPU/RSS; the parent writes actual Slurm placement and stage
duration. Root arms the existing completion wake on the submitted job.

Metadata preparation and source/syntax review may occur on the login host.
Scientific parsing, precision checks, evaluation and algebra run only in Slurm.
Use a fresh immutable workspace and lock; no automatic restart. Any missing or
failed required cell leaves its four-cell factorial unavailable, retaining the
diagonals and other results. Preserve partial results and execution receipts.
Historical absolute-v1 parity failures remain explicit beside accepted component
checks; the completed4623 result and every prior rejection stay unchanged.
