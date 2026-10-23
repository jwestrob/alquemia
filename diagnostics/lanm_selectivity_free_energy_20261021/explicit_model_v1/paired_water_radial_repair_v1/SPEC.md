# Prospective radial-inclusive paired-water development

## Question and preserved outcomes

**Design declaration only: no implementation, fitting, native calls or validation opening yet.** Can the same four radius/damping parameters satisfy both individual metals and their difference when the consumed radial water configurations enter development explicitly?

Job17576 passes all four-center development controls. Frozen transfer17578 passes all10 Cartesian force-consistency checks and all water energy checks, but the compressed water fixture fails rigid response: La RMS/max2.39566/4.759, Dy6.43945/12.7679 and paired4.04577/8.0088kcal/mol/Å. Dy is worse than its original3.739/7.418 baseline there. Expanded water passes. Preserve both results: center success does not imply radial transfer, and this finite failure does not establish impossibility of the four-parameter family. No donor fitting around a failed water control is included.

## Exact eight consumed cells and fixed model

Use only `water_{minus,central,plus}_{La,Dy}` and `aquo8_central_{La,Dy}`. The four water-side geometries already exist at−0.15/+0.15Å; in those original fixtures the **metal** moves along the stored water-to-metal axis, while water atoms remain fixed. Do not replace these with the later rigid-water Water48 path or regenerate any coordinates. Four geometry pairs, eight actual cells, remain the full denominator.

Start exactly at17576's selected vector `[3.710148544089586,0.2001607792689767,3.389554168967775,0.16712500000000005]`, ordered La Rmin/damping, Dy Rmin/damping. Reuse all eight actual candidate E/force records admitted by17578. Preserve the original model baselines separately. Fixed epsilons remain La0.24429636467571708 and Dy0.37332330078584813kcal/mol. Original radius bounds remain La[3.1437505772349654,4.253309604494365] and Dy[2.82379559590138,3.820429335631278]Å; damping bounds remain[0.05,0.50]. No donor, water-H, bonded, intrinsic-polarization, permanent-multipole, geometry, polarization-tolerance or reference-setting change.

Pinned starting evidence:

- 17578 transfer RESULT: `workspaces/lanm_selectivity_free_energy_20261021/explicit_model_v1/paired_water_transfer_v1/run_v1/execution/RESULT.json`, SHA256 `702d97ababc232336e3b5405a1326cec4cafeccc9570ce9f514076cf7908e986`.
- Its BOUND inventory: `.../paired_water_transfer_v1/bound_v2/BOUND.json`, SHA256 `d7292365886c131f1c81b67abc6348b04f084531c4ceae4a3006926d4f17144d`; this binds exact original/candidate inputs, consumed radial reference E/forces, geometry and states.
- 17576 selected RESULT SHA256 `fab7016c4335a82dea983368abc6f23201dd83df0324f7794e5e0c978fd8e190`.

All targets stay on the same **uncorrected PBE0-D4** energy/force surface. CC/CP and Water48 method sensitivity remain separate diagnostics, never mixed into this force objective. Frozen-f/electronic-state and physical-reference uncertainty persist.

## Frozen symmetric objective and controls

Keep four equally weighted blocks:0.25 energy,0.25 individual La rigid response,0.25 individual Dy rigid response,0.25 matched Dy−La rigid response. Energy scale1kcal/mol; force/torque scale2kcal/mol/Å. For each force block, give the four geometry pairs equal weight. Within a fixture use unchanged net force and geometric-centroid torque/1Å, equal ion/water category weights, equal waters within category and equal translation/rotation for each water. Do not fit shared internal-water force distortion.

The energy block is **0.50 mean-square balanced exchange residual plus0.50 mean-square individual radial-work residual**:

- Three exchanges: C(q)=D(water_q)−D(aquo8), q=minus/central/plus, D=E_Dy−E_La. Each reaction contains one ion of each metal with matching total waters.
- Four individual works: E_M(water_side)−E_M(water_central), M=La/Dy and side=minus/plus, each equally weighted.

Only model-minus-reference residuals enter these means. The energy blocks share endpoints and are correlated engineering weights, not independent statistical observations. No raw unlike-element energy or offset is fitted. Report the two paired radial-work residuals as separate controls as well.

Stop only when **all** existing per-fixture controls pass: every individual and paired geometry RMS≤2 and max body-component error≤5kcal/mol/Å; all three exchanges, four individual works and two paired works have absolute error≤1kcal/mol. Preserve every existing center control and per-fixture result; no average may conceal a failing compressed cell. Also retain the original four-center objective/controls as diagnostics, distinct from the new declared loss. No tolerance or force/energy weight changes after execution.

## Finite joint refinement proposal

Reuse the successful bounded joint method, starting at the exact admitted eight-cell baseline. At most12 outer steps. Each step evaluates eight clipped central parameter displacements±0.005 normalized interval widths, dividing by actual clipped spans. Build the weighted-residual Jacobian and use SVD relative cutoff1e−6. Produce the same uniform-scaled Gauss–Newton and bounded linear-least-squares directions within original bounds and component trust≤0.025. Evaluate scaledGN half/full then boundedLS half/full; finish every derivative and trial batch before selection. Current wins ties; strict lowest original-for-this-branch objective wins. Record any passing-but-unselected row.

Stop when the selected point passes all controls, no actual improvement occurs in a complete trial batch, or12 steps finish. No fallback, family expansion, extra sensitivity selection or automatic continuation. A controls-pass stop is not a global optimum; a local/budget stop is not proof of family failure. Report Jacobian rank/condition at the actual iteration origins and bounds; do not label an earlier Jacobian as final-point qualification.

Ceilings: **144 new parameter vectors and1,536 new native program invocations**. Each derivative changes only one metal:8 derivative vectors×4 represented-metal cells×2 programs=64 programs/step. Four joint trials×8 cells×2programs=64 more;12×128=1,536. The initial eight cells are reused with zero calls. Exact physical-cell cache reuse can lower this ceiling. Full-vector request identity remains separate from represented-metal cell identity; changing La must reuse Dy and vice versa. Every new cell uses one analyze and one analytic testgrad, one thread. An8CPU allocation supports at most8 native processes; root owns submission. First native failure stops new starts and preserves siblings/denominator; no retry. Aggregate child CPU is measured once, not summed over overlapping receipts.

## Minimal source reuse and next gate

The original evaluator cannot be used byte-for-byte: it hardcodes four cells, the original baseline vector and its source/config identities. Make a contained new version with eight explicit cells and an admitted17578 baseline, preserving the renderer, parser, native invocation, per-cell locking and failure logic. Reuse the exact bound input graphs/types/reference arrays. The full Cartesian center artifacts and analytic side artifacts have different native-record shapes; admit each with its declared parse mode, without rerunning them. Import prior records only after exact scientific input/raw replay and retain their original authority pins; never relabel an old cache key as new or alter old ledgers.

Reuse the old rigid-body helper and joint direction routine. A small new objective must implement the explicit seven energy residuals and four-fixture force blocks; a thin new driver changes only finite inventory/baseline/loss/caps, preserving the numerical search order. Fresh outputs and ledger; immutable17576/17578 and all earlier caches. Zero-call preparation must replay the starting eight rows and new score plus self-target/force-sign/paired-cancellation controls before root release. No generic fitter framework is needed.

If a new candidate passes these consumed controls, separately declare its derivative and frozen transfer/validation checks. Reserved24 orientations/protein outcomes remain unopened. No protein, periodic-water, bulk hydration, affinity or physical-accuracy qualification is implied.
