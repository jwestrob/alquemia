# After the matched pair result

**Decision: test the reference remainder in a common full-cluster basis before changing the collective model response.** Prepare one finite, energy-only ghost-fragment diagnostic: 68 new electronic single points, all four original complexes and all eight completed CP fragments reused, zero new model calls. This is a prospective design only; root has not released it in this event turn. No parameter fit, donor extension, target change or held-out test is proposed.

## What the actual decomposition establishes

The complete 36-reference/34-model collection gives, in kcal/mol:

| Paired W9−W8 contrast | Reference | Fixed model | Model minus reference |
|---|---:|---:|---:|
| Full C | −0.534207 | −4.872491 | −4.338284 |
| Isolated-pair P | −8.933375 | −7.670058 | +1.263316 |
| Cooperative remainder B=C−P | +8.399168 | +2.797568 | −5.601600 |

The larger remainder discrepancy is opposed by the pair discrepancy. A pair-only correction that removes +1.263316 while leaving B unchanged would expose the larger −5.601600 discrepancy; additive compensation would not repair the inferred collective response.

Saved native component arithmetic confirms what **the model's** B contains:

| Model component | C | P | B |
|---|---:|---:|---:|
| Buffered 14–7 vdW | −5.29009248 | −5.29009248 | approximately 0 |
| Mutual polarization | +0.41760198 | −2.37996595 | +2.79756793 |
| Permanent multipoles and all reported bonded terms | 0 | 0 | 0 |

The full-energy B is +2.79756794; its approximately 1e−8 difference from the component sum is printed precision. The pair-additive vdW term cancels from B to approximately 1e−14. This is the expected source-level identity, now supported by the 34 literal dimer energies. It does **not** make the reference B a uniquely defined polarization energy. The reference splits into electronic B=+8.47600205 and D4 B=−0.07683414; dispersion nonadditivity is small in this decomposition, whereas electronic finite-basis and collective-response effects remain combined.

All fragments used literal source coordinates. W9 also rearranges the first eight waters, so these are coordination/arrangement contrasts, not a pure ninth-water binding energy. The errors concern these consumed development configurations and the provisional frozen-f PBE0 surface.

## Why the remaining reference test is different from completed CP8

CP8 corrected metal-versus-whole-shell basis sharing, changing full C by +0.32472251 to −0.20948415. The unchanged model then misses that alternative by −4.66300635. It did not qualify the basis dependence of the isolated-pair sum. Raw dimers had 90 orbital basis functions; the original W8/W9 complexes had 391/434. Consequently the apparent +8.399168 cooperative remainder can include the benefit of basis functions present only in larger fragments.

The proposed test places every fragment within its parent cluster's **same complete orbital/auxiliary basis**. It therefore addresses whether the present allocation of error to B is a fragment-basis artifact. It cannot rescue the already computed full C merely by repartitioning it, and does not settle basis-limit, functional, frozen-f/state or spin–orbit uncertainty.

## Exact common-basis algebra

Let `Delta X = X_Dy,9 − X_La,9 − X_Dy,8 + X_La,8`.
For each parent `(M,N)`, all superscript `N` energies below use all of that parent's basis centers at their exact original coordinates, with inactive atoms represented as ghosts:

- `E_MN`: original full complex electronic energy, reused.
- `a_MN^N`: real ion with the entire water shell ghosted, reused from CP8.
- `s_MN^N`: real complete water shell with the ion ghosted, reused from CP8.
- `d_MNi^N`: real metal plus intact water i; every other water ghosted, new.
- `w_MNi^N`: only intact water i real; metal and all other waters ghosted, new.

Then define electronic quantities consistently within each parent basis:

```
C_common = Delta [ E_MN − a_MN^N − s_MN^N ]
P_common = Delta [ sum_i=1..N (d_MNi^N − a_MN^N − w_MNi^N) ]
B_common = C_common − P_common
         = Delta [ E_MN − sum_i d_MNi^N
                   + (N−1) a_MN^N − s_MN^N + sum_i w_MNi^N ]
```

The ion term in P is subtracted **N times**, once per pair. Do not substitute the newly calculated own-basis ions or assume that ghost-basis ion/water/shell terms cancel between different parents. The full-basis shell and each water also carry the appropriate metal ghost basis, so their La/Dy values need not agree.

This C_common is exactly the previously admitted CP8 electronic full-shell contrast. Retain the original D4 contributions separately and unchanged: C_D4=+0.10663833, P_D4=+0.18347248, B_D4=−0.07683414 kcal/mol. All 68 new calls are electronic PBE0 **without D4**, as in CP8; any displayed electronic-plus-original-D4 sensitivity adds the corresponding saved D4 quantity only once. Verify closure independently for electronic, D4 and their sum. Do not label a mixture of raw P and CP-corrected C as this common-basis decomposition.

## Finite run inventory and numerical admission

| New role | Count | Real electrons / charge | Basis retained |
|---|---:|---|---|
| M + water i; remaining waters ghosted | 34 = 2×(8+9) | 18 / +3, effective singlet | Entire original W8 or W9 basis |
| Water i only; ion and remaining waters ghosted | 34 | 10 / 0, singlet | Entire original W8 or W9 basis |
| New full complexes, ion/shell fragments, model calls, gradients | 0 | Reuse admitted results | No repeated chemistry |

Use the already tested literal maps from the pair panel, with **no rotation, averaging, new distance or geometry optimization**. Every W8 cell has 25 basis centers/391 AOs and every W9 cell 28 centers/434 AOs. A real La/Dy retains its exact 46/55-electron ECP and metal orbital/auxiliary basis. A ghost metal has zero nuclear charge, no ECP potential, and retains its exact orbital/auxiliary basis. All ghost waters likewise retain basis functions but contribute no electrons, nuclei or dispersion. Expected atom tables, basis counts, charges, ECP ownership and electron inventories must be checked against native output.

Preserve the original PBE0, RIJCOSX, DefGrid3, SCF thresholds/iteration ceiling and frozen-f effective state. Use the qualified CP8 native admission, retain actual stopping and final-integration metadata, and require all 68 new fragments plus the 12 reused complex/CP records before a complete decomposition. A failure remains explicit; no automatic retry, replacement state or partial sum.

Resource proposal: **68 native invocations maximum, four single-thread MPI ranks per cell, standard_expand, default memory/time and allocation-derived MaxCore**, scheduled only by root within current shared limits. Independent cells may run concurrently within those limits. The 391/434-AO ghost bases make these more expensive than the 90-AO isolated dimers; previous small-dimer timings are not a valid cost prediction. Record actual per-cell and total costs. There are no new fitting vectors, force calculations or model programs.

## How the outcome changes the next mechanics decision

Compare the reference B shift with the current 5.601600 kcal/mol B discrepancy and the existing 1 kcal/mol development energy resolution; neither number is an uncertainty estimate. If common-basis B moves substantially, the previous component attribution was basis-sensitive. Recompute both pair and remainder errors without changing targets; the full-C discrepancy remains the completed CP8 value. If B remains near its present value, the observed model deficit in metal-dependent collective response survives this particular basis-sharing test and warrants a separately declared response/functional-form experiment. It still does not uniquely identify polarization, charge transfer, exchange penetration or nonlinear susceptibility as the missing physics.

The completed metal-alpha test does not justify adding alpha to a fitter: its actual ±1/±5% probes gave no aggregate improvement, and favorable extrapolation required invalid/out-of-bounds parameters. H-Thole probing remains low priority. The consumed donor fit17621 also stopped without actual local improvement and still failed 21/33 energy checks and 27/36 force fixtures. Those independent development failures strengthen the case for testing the interaction response rather than more optimizer budget; they do not prove the same cause as water B. Donor qualification, fresh water validation and bulk sampling remain held.

## Evidence and scope

- Matched collection RESULT: `reference_chemistry_v1/matched_pair_decomposition_v1/collection/run_v1/result/RESULT.json`, SHA `9ec1b00961fb7c0309062d187f4c79f298aedb489c23ae89c7848bd3c21d9dac`.
- Model raw-row aggregate: `matched_pair_decomposition_v1/model/run_v1/execution/RESULT.json`, SHA `a223bfbb7ed7abd934ec74a89b9529dba16280f7a63676bc903348fb593b0ab6`. The table above is the signed sum of its saved per-cell component energies, with the four saved cluster rows reused.
- CP8 RESULT: SHA `3647298ae095a95ff1d2149eea2d273432c933cf33bc41e9928abc4a0c811e79`; original selected17676 model/result and prior alpha evidence stay pinned in the preceding notes.
- Donor final report: `explicit_model_v1/paired_donor_fit_v1/FINAL_REPORT.md`, SHA `69f269ee4a9bd1cf507d0c9feeebf12aad7e4bdfe2f8ebc75458a0f4dea716a3`.

This note uses saved arithmetic/source evidence at one local thread. No native call, submission, source implementation, target change or held-out outcome read was performed. Root/independent reviewers own actual collection admission; this is a concrete next-step proposal, not execution authorization or physical qualification.
