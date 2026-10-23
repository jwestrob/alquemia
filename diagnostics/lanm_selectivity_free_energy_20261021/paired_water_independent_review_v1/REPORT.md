# Independent paired-water pilot review

**No material unresolved issue in the reviewed source bundle.** This accepts the implementation of the finite representability experiment; it neither releases execution nor qualifies the model physically. Root owns the execution release.

## Scientific arithmetic

The four parameters are the separate La and Dy water-O pair radii and damping values, with both resolved pair well depths fixed. The scalar energy contrast is
`[E(DyW1) - E(LaW1)] - [E(DyW8) - E(LaW8)]`, corresponding to `DyW8 + LaW1 → DyW1 + LaW8`. Its stoichiometry is balanced, so no unlike-element or free-fragment offset is fitted. The exact coordinates, covalent bodies and atom order match across metals within each geometry.

The objective uses four equal 0.25 blocks: balanced energy, La rigid response, Dy rigid response and paired rigid response. Its force blocks retain equal geometry/category weights, equal translation/rotation weights per water, and a 2 kcal/mol/Å scale. Water torque is computed about the geometric centroid and divided by 1 Å. Individual and paired blocks are correlated development weights, not independent likelihood observations. The separate pass conditions include both metals and the pair.

I independently recomputed the saved baseline using explicit sums of atomic-force residuals, centroid torques, category means and geometry means. The three force blocks agree with the adapter within **3.6e-15**; the complete scalar objective agrees within **2.4e-11**, from energy summation order. The baseline objective is **5.238728922966** and balanced energy error **−1.621251017 kcal/mol**. La/Dy/paired rigid RMS errors are **6.163346 / 1.462036 / 5.760335 kcal/mol/Å**; all three fail their separate controls, including the Dy maximum-component condition.

Saved Water48 sensitivity replay uses the outward metal-to-water-O unit vector and **−sum(water forces)·u**. The baseline La/Dy/paired derivatives reproduce **12.138589780 / 32.808390970 / 20.669801189 kcal/mol/Å**. The older training geometry's inward radial vector is not used. Uncorrected PBE0-D4 energies and forces remain the objective; CP/CC differences are separately reported sensitivity checks.

## Execution, cache and selection

- The renderer changes only the represented metal's water-O VDWPR/POLPAIR records. Water-H pairs, non-water pairs, atom parameters, charges, polarizabilities, geometry and covalent terms are unchanged. Native resolution proves baseline equivalence and fixed epsilons before baseline reuse.
- Full-vector requests retain all four values. Native cache identity uses the active metal's exact rendered input and complete source/cell lineage. A change to La alone reuses unchanged Dy cells, and vice versa. Per-cell locks prevent concurrent duplicate native calls.
- Complete polls finish before selection. Fixed poll order controls ties, independent of completion order. The search stops before an incomplete budget-limited poll; eight final sensitivity vectors are reserved and cannot replace the selected point. Clipped sensitivities divide by the actual normalized separation.
- The 241-vector and 1,920-native-program limits remain enforced. The driver requires one eight-CPU Slurm task, with at most eight workers and one native process per worker at a time. Numerical thread controls remain one.
- A material initial sequencing issue was corrected before this review was frozen: native exit-zero failure markers, missing/nonfinite energies and missing analyze components are now checked after each program, before releasing its native slot or scheduling the next program. Native parsing/failure sets the shared stop flag; no fallback or retry is added.
- A second correction records one aggregate start/end child user/system CPU difference. Overlapping per-worker child-CPU deltas are not summed. Parent CPU, wall time and per-program receipts remain separate.
- Actual zero-call preparation, pinned configuration and an explicit root release gate execution. A budget-limited result does not establish family failure or optimizer convergence. Rank, condition and bound flags do not become parameter uncertainty or physical validation. Aquo8 symmetry makes many force components correlated; identifiability is tested, not inferred from array length.

## Evidence and limits

Author controls cover deterministic complete polls, reordered callback completion, local stop, refusal to begin a partial poll, failure stopping, and sensitivity non-reselection. Evaluator controls cover four actual native baseline replays, effective-input caching and active-metal-only rendering, invalid inputs, exit-zero native failure markers, absence of release/allocation, and persistent stopping. All passed with zero native calls. I reviewed their source and saved results rather than rerunning them.

Exact reviewed source and control pins are recorded in the matching workspace `REVIEW.json`. Frozen source acceptance is distinct from actual preparation/configuration admission and root submission. No reference call, fitting call, validation outcome, donor model change or dynamics was performed by this reviewer.
