# Design review: whole GFN base plus a regional MACE correction

19 October 2026. **Scientifically coherent and worth a bounded fixed-fixture test while the 7CCN reference runs.** It is a new mechanically embedded candidate, not a demonstrated repair or a reason to restart optimization. Source implementation and actual runtime gates remain to be reviewed. No numerical work, imports or submissions were performed.

## What the model does

The proposed energy is

`GFN2_vac(full) + MACE_vac(region) − GFN2_vac(region) + fixed_FF_ALPB(full)`.

This has the standard subtractive two-layer structure documented by the [xTB developers](https://xtb-docs.readthedocs.io/en/latest/oniom.html). The full GFN calculation supplies the scaffold and interactions across the regional boundary, replacing the challenged MM ledger. Its self-consistent electrostatics and dispersion are explicit parts of GFN2, as described in the [native API documentation](https://xtb-python.readthedocs.io/en/latest/general-api.html). The observed absence of these explicit terms from the executed whole-MACE path motivates this hypothesis but does not prove their absence caused the measured error.

**Self-consistency remains at the GFN level.** The vacuum regional correction does not see the environment field or feed a higher-level density back into full GFN. Fixed-charge ALPB also remains additive. This is not mutually self-consistent MACE/GFN polarization. Subtracting an isolated capped region is an approximation to a local replacement, not an exact decomposition of the full GFN density. Identical caps still leave method-dependent boundary error. No old embedded-minus-isolated field term, MM term or extra dispersion belongs in this energy.

## Required state and cap identities

- Preserve the exact full 256-atom fixtures and terminal chemistry. Full charge is −2. Use the existing baseline region: the inspected map has 172 atoms including two hydrogen caps, also charge −2. Do not substitute the older 48-atom, charge −1 local-reference convention.
- Both regional evaluations must share identical symbols, ordering, coordinates, charge, cap lengths and connectivity. Preserve physical La multiplicity 1/Dy 6 for MACE and the published GFN effective multiplicity 1/UHF 0 convention separately. GFN frozen-f valence treatment is not the ORCA ECP46/55 implementation. Derive and check actual full/regional native electron counts; correction-profile qualification on earlier cores does not establish whole-peptide state accuracy.
- Pin the corrected La/Dy profile `635e111b…`, native executable and imported parameter evidence. Both GFN calculations must be vacuum, without external point charges or native solvent. Whole and regional GFN must use the same declared scalar-energy/gradient convention, electronic temperature and solver settings. In particular, do not mix an SCC iteration energy with a different final energy that includes additional terms.
- Rebuild cap coordinates and Jacobians at each fixed endpoint using the existing mapping rule. Both retained and omitted anchors receive chain-rule contributions. With `J` mapping full coordinates to the capped region, the gradient is `g_GFN_full + Jᵀ(g_MACE_region − g_GFN_region) + g_ALPB`. Frozen connectivity does not mean frozen cap positions. Source and replay checks must preserve the chosen arithmetic and raw-artifact authority; prior cap-scatter roundoff failures must not be erased by an improvised tolerance.

## Smallest discriminating comparison

Declare the full eight-state ceiling now: four consumed 7CCO states, then at most four already selected 7CCN challenge states. Freeze the 7CCN geometry, partition, state, settings and comparison rules before its new DFT outcomes; do not select or retune the candidate on that challenge. Maximum eight full-GFN, eight regional-GFN and eight regional-MACE evaluations are ceilings, reduced by exact compatible reuse. No new DFT, solvent, geometry search or optimization is needed.

**Whole GFN plus the same ALPB alone is a mandatory comparator and needs no additional calls.** Report it beside the regional-corrected candidate, existing whole-MACE, the available hybrid and whole-PBE0 reference. This distinguishes a useful base from a useful regional correction. Retain each component's before/after work, both per-metal errors and their differential. Report full 256-atom force errors and error changes, existing chord projections and both terminal regions; a selected projection or cancellation in the metal differential cannot replace those comparisons.

## Finite gates and decision

Root's proposed requirement that both 7CCO per-metal work discrepancies remain within 1 kcal/mol of whole PBE0 is defensible as a **new, predeclared development triage rule** for spending on the second fixture. It is stringent; it is not a DFT error bar or an accuracy qualification. Do not relax it after observing either metal. Report the force evidence even when this work gate passes; it supplies no automatic optimizer approval or new numerical force-accuracy threshold.

Stop on any missing component, failed original native SCC/parameter/state admission, nonfinite force, broken source/cap/component closure or exceeded declared call/resource budget. Preserve the failed slot; do not form a survivor comparison, change temperature/precision/solver limits, retry chemistry or alter the partition inside this test. If either 7CCO work gate fails, stop the candidate before the remaining four calls. A failed or incomplete 7CCN reference leaves its comparison unavailable rather than licensing another reference or favorable-subset result.

A pass warrants only the prebound challenge. Even agreement on both fixtures would remain conditional peptide evidence; numerical derivative qualification of the new composite and any broader accuracy requirement would precede a separately declared optimizer. The pending 7CCN reference therefore need not idle this small consumed-fixture screen, subject to measured scout cost, but no alternate-region sweep or sequence of candidates chosen from its eventual outcomes is justified.

## Evidence inspected

The saved [subtractive-design discussion](../whole_system_intervention_v1/INDEPENDENT_DESIGN_REVIEW.md), SHA `74a637ad363369d350a1ec40a0133a826857638c8cb2b0f08b7c5abd6a36e2af`, and [native physics review](../whole_model_physics_review_v1/SOURCE_PHYSICS_REVIEW.md), SHA `2849ace719a930c4804023f392a4fae98baeddb2d49e640c8fffba09c168d7a0`, were read. Existing baseline cap map `312a4e84…`, native isolated-GFN task `957bd524…`, the corrected-parameter report and preserved embedded-GFN failure were also inspected. An independent source-based cross-check by `/root/mechanics_mex_adapter` agreed on the limited self-consistency and state/cap conditions above.
