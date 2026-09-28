# Carbon-bound H repair: geometry improves on all three actual sources

**The bounded native bonded+steric preparation removes most verified C-bound H distortion at low cost, with unchanged heavy atoms, waters, exchangeable H and proton inventory.** Actual saved force/geometry checks support a separately named repaired-input electronic comparison. Original automatic admission remains unavailable because optimizer metadata was lost during result serialization; a generic signed-volume check also incorrectly identifies two nonstereogenic methylene centers. These limitations are preserved rather than silently changed to passes.

All three optimizations finished under the frozen1,000-iteration/5,000-evaluation limits. Complete coordinates, before/after forces and actual objective XML were written before JSON serialization failed on `numpy.bool_` in a check. Partial JSON preserved actual evaluation counts, energies and before/after angle audits. `RECOVER.py` reconstructs geometry admission from those immutable files without any additional molecular call. It does not invent optimizer flags, iteration counts or per-source wall times.

| Source | Free C-bound H | Actual force/energy calls | H-angle energy before→after, kcal/mol | Maximum free-H gradient, kcal/mol/Å | Max repaired angle error |
|---|---:|---:|---:|---:|---:|
| Hans8DQ2 |630|114|4846.851→674.714|.04315|6.388°|
| Hans8FNR |630|87|3143.874→586.648|.04051|6.593°|
| Mex8FNS |597|105|1904.494→546.201|.04665|6.381°|

H-angle totals include unchanged exchangeable H and water angles; residual values are not evidence that optimized C-bound H remained badly oriented. All free-H gradients meet both the proposed.1 admission threshold and the optimizer's.05 target directly in saved analytic forces. All frozen coordinates and source IDs match exactly; maximum C–H length error is.0226/.0334/.0210Å and displacement1.281/1.448/1.194Å. No new nonbonded free-H contact below.65Å appears. This is substantial motion of generated hydrogen coordinates, not protein heavy-atom accommodation.

## Stereochemistry and admission accounting

The original check indiscriminately compared signed volumes at every four-neighbor heavy center. It flags Hans8DQ2 Asp38CB and Hans8FNR Met132CG; each has **two chemically equivalent H neighbors**. Such CH2 centers are not stereocenters. Their named-H volume reversal is a label/geometry issue, not inversion of molecular chirality. All actual four-distinct-substituent centers pass; Mex has no volume flags. The original `no_inversions=false` remains in the recovered Hans check records. No favorable score motivated this diagnosis, and no electronic energies were evaluated.

The original requirement `optimizer_success` cannot be recovered because serialization failed before that metadata was written. `RECOVERED.json` therefore uses `geometry_recovered_original_admission_unavailable`, with that field null. The objective force residual and all other geometry evidence are independently assessable. Root owns whether to use these saved geometries as an explicitly reviewed candidate. No extra minimization was performed just to recover metadata or change a failed criterion.

## What this means scientifically

The objective is exactly native bonds/angles/torsions/CMAP plus native LJ, with all Coulomb terms zero and all ions absent. It is a cheap geometry-preparation surrogate, not a total force model, binding energy, metal-dependent correction or solvent-consistent protein relaxation. It must not be added to QM energies. Heavy atoms, the real water inventory, exchangeable H and source topology remain fixed; no protonation/occupancy state changed. No synthetic cap was optimized.

This resolves a concrete input weakness that affected both exterior and full QM region. It does not prove improved La/Dy discrimination or eliminate source-dependent chemical states. Any next electronic comparison must use a new identity, source-mapped repaired A, the same declared intact Asp rotation to derive B, regenerated external fields/maps, and both target metals on identical repaired coordinates. Never mix original and repaired endpoint energies. Compact transfer requires actual source-ID mapping, not positional matching or silently reusing the larger region composition.

## Artifacts, costs, technical recovery

`workspaces/metal_environment_response_20260926/scaffold_H_repair_v1/{Hans8DQ2,Hans8FNR,Mex8FNS}/` contains `coordinates_A.json`, `forces.json`, `system_objective.xml`, incomplete original `RESULT.json`, and complete `RECOVERED.json`. `COLLECTION.json` retains all three original serialization errors; `RECOVERY_COLLECTION.json` retains direct checks. Source inputs remain unchanged.

Job1220315 completed in32wall seconds on1sharedCPU:32allocatedCPU-seconds, noGPU, mem0. Actual molecular calls306 total, two recorded states per source plus bounded optimizer evaluations. Process peak RSS141,484KiB; Slurm sampled83,372KiB. Serial recovery/angle arithmetic outside the allocation adds unmeasured small CPU time, no molecular calls. Actual scheduler receipt: `ACCOUNTING.txt`.

Original `RUN.py` and plan are immutable execution evidence; **do not rerun this known serialization-bug version**. A future implementation must convert NumPy scalars to JSON-native types and persist optimizer receipts before aggregate checks. No repeats are needed to inspect the saved results. Vault updated; production, old inputs and active root electronic jobs untouched.

## Subsequent explicit chemical review

Root reviewed the saved evidence and authorized a separate chemically reviewed admission. `REVIEW_ADMISSION.py` admits all three geometries from their actual saved force/geometry checks. It excludes only flagged carbon centers with two chemically equivalent H neighbors from a stereocenter interpretation; every other center retains the original volume test. Missing optimizer status remains null, and the original gate records are unchanged. The saved gradient is direct geometric-stationarity evidence, not a recovered optimizer receipt. No observed preference influenced this correction.

Each source now has `REVIEWED_ADMISSION.json` pinning coordinates, forces, full source export, original/recovered records and permitted C-bound H IDs; the aggregate of the same name lists all three. The independent compact-transfer agent has these pins and owns preparation-only mapping into new compact fixtures. Root owns any electronic execution. No additional optimization or molecular evaluation occurred during review.
