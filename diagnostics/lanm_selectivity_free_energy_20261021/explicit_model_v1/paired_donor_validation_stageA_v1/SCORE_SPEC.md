# Stage A orientation validation scoring

Source-only declaration, before candidate binding or any reserved numerical readout. Score exactly the existing 16 orientation cells: aquo8, hydrated mono-acetate, bi-acetate and acetamide, each minus/plus and La/Dy. No geometry, parameter, reference, weight or threshold change; no fit objective or candidate selection.

## Interface and provenance boundary

`evaluate(data, model, body, anchors)` is pure arithmetic. `data['rows']` has exactly16 metadata rows, each with `geometry`, `input.bonds_zero_based`, `bodies` and `reference`. `model` has exactly the same16 IDs. Reference/model records contain energy in kcal/mol and forces in kcal/mol/Å, optionally `status='complete'`; `None` or an explicit `unavailable`, `not_attempted` or `failed` status denotes missing data. An incomplete record cannot supply a partial energy or force result. Invalid finite values, geometry, chemical state, graph or paired atom mapping raise an admission error.

`anchors` must explicitly contain `model` and `reference`, each with La/Dy records giving the fixed **aquo8 central** energy in kcal/mol. Missing anchors remain explicit unavailable slots; they do not suppress orientation-work or force comparisons that do not depend on them. The integration adapter owns exact source/candidate/reference/raw-output admission and unit conversion. It must pin the original reserved manifest, the frozen passing candidate, both actual central-aquo model records and both unchanged central-aquo reference records. The scorer reads no files and chooses no default or historical anchor.

Use unchanged rigid-body helper `matched_dry_v1/OBJECTIVE.py`, SHA256 `bf65a20a379811f8d35e2cb5ae78365a5d34022691b332bae5df223196078819`. All16 geometry slots remain present even when numerical records are unavailable. Missing slots have `passed=null`, and overall validation requires complete16 cells, all four anchors and all44 declared comparisons. No missing-as-zero, dropping a slot, averaging away a failure or automatic retry.

## Twenty energy comparisons

Define `D_X(g)=E_X,Dy(g)−E_X,La(g)` for model/reference X. For each of eight configurations, compare

`C_X(g)=D_X(g)−[E_X,Dy(aquo8 central)−E_X,La(aquo8 central)]`.

Report `C_model−C_reference`; every absolute error must be ≤1 kcal/mol. This balances the metal substitution offset without introducing free fragments.

For each of four families and each metal, compare the orientation work `W_X,m=E_X,m(plus)−E_X,m(minus)`. Report eight individual work errors and four paired errors `[W_model,Dy−W_model,La]−[W_reference,Dy−W_reference,La]`. Each absolute error must also be ≤1 kcal/mol. The sign is always plus minus minus, then Dy minus La, then model minus reference. No radial or central-to-side work is invented for this orientation panel.

## Twenty-four rigid-response comparisons

For each of eight configurations, compare La, Dy and paired Dy−La forces using the original covalent-body mapping and helper. Individual metals remain separately required even if a shared force error cancels in the paired difference. The ion contributes net force; each intact ligand/water contributes net force and torque about its geometric centroid, divided by1 Å. Categories have equal weight, bodies within a category equal weight, and translation/rotation equal weight as in the helper.

For helper residual vector r (already divided by2), rigid RMS is `2*sqrt(sum(r²))`. Require RMS≤2 and maximum absolute unweighted body component≤5 kcal/mol/Å for **each** of24 fixtures. Report model/reference body responses and atomic residual diagnostics; internal atomic residuals are not an extra fitted or pass criterion. No aggregate loss or block mean determines acceptance.

## Interpretation and bounded source checks

Outputs retain16 cell availability records,8 balanced energies,12 orientation works and24 force fixtures, plus fixed declared/available/passing counts. A calculated failure is distinguished from unavailable evidence. `validation_passed` requires every declared comparison and input; it remains false for a partial panel. Passing is conditional on the admitted reference surface and these configurations, not protein/bulk/affinity validation. Preserve frozen-f ECP/state and reference derivative/method caveats from the parent protocol.

Validate source arithmetic using reserved geometry **metadata only** and explicitly synthetic energies/forces: zero self-error; signed single-cell energy perturbation; metal-offset cancellation; plus/minus work; shared-force paired cancellation with individual failures; missing cells and anchors; exact denominator; invalid/nonfinite records and geometry. Do not read reserved native/model outcomes or run engines for these controls. Candidate binding and execution remain root-owned and separate.
