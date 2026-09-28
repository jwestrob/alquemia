# Where the hydrogen defects occur and a bounded repair proposal

## Actual source tracing

The defects occur inside and outside the electronic region. `H_ORIGIN_AUDIT.json` partitions native H-containing angle strain:563.917kcal/mol wholly QM,12.793cross,4270.141MM. The largest QM defects include Asp83 CA–CB–HB2 at166.711° (native109.5°), Asp83 HB3–CB–CG161.334°, and Met92 CA–CB–HB2158.409°. These are real source-mapped atoms, not synthetic caps.

The actual moving Asp85 donor is less distorted: all its local H-angle terms total5.241kcal/mol, with CB-related H angles approximately101–119°. Its +2° motion preserves those angles. The active electronic response therefore includes badly oriented neighboring QM hydrogen chemistry, but this audit does not demonstrate that the Asp85 motion itself creates a covalent defect or identify how much electronic response changes after repair.

Original and normalized H angles in these examples agree to numerical precision. The archived protonation receipt is `workspaces/benchmark_set_20260915/prepared/hans_lanm_v1/hans_protonation_manifest.json`: PDBFixer1.12.0/OpenMM8.5.1, pH5, no missing heavy/terminal atoms added, original8DQ2 source hash retained. Current `protonate_cif.py` calls PDBFixer `addMissingAtoms`, then `addMissingHydrogens(ph)` and writes stable IDs. The receipt does not pin that historical script, so its precise historical implementation cannot be inferred from today's source alone. The bad directions demonstrably predate the subsequent whole-chain preparation because they are present in the actual archived protonated coordinates.

The pinned whole-chain preparation calls `mace_global_prepare.protein`, which uses native equilibrium X–H lengths and `mace_hydrogen.repair`. That helper rescales the existing heavy→H vector; it neither reconstructs hydrogen angles nor minimizes them. Water OH lengths are separately normalized. Bond-length checks passed but did not test these directional defects. Do not rewrite either the archived protonated source or active electronic inputs.

## Proposed versioned preparation diagnostic — not executed

First repair only **carbon-bound hydrogens** under the actual source graph. This directly targets the largest verified tetrahedral defects without simultaneously reorienting exchangeable proton donors or coordinating waters. Keep all source heavy atoms, metals, water coordinates, exchangeable N/O/S-bound H, atom IDs and proton inventory exactly fixed. Synthetic caps continue to be derived from unchanged heavy source cuts; no cap degrees of freedom.

Use the exact native ff19SB/TIP3P parent System as a cheap H-coordinate preparation objective, with only selected C-bound H positions free. This is an explicit ion-free preparation surrogate, not a physical metal-site relaxation or a hybrid energy; missing target/spectator interactions stay visible. Use analytic native forces and a bounded H-only minimizer with a declared iteration ceiling and maximum accepted H displacement before any run. No arbitrary restraint constants or Hessian fit. Avoid deleting/readding hydrogens, which could change identities/protonation; reuse the current maps and unchanged native parameters. The prior source export/decomposition validates mechanics bookkeeping, not the repaired coordinates in advance.

Admit only if heavy/frozen coordinates remain exactly identical; each H keeps its source bond; native bond/angle strain decreases; no new heavy/H clash or stereochemical inversion appears; minimizer residual/status and any remaining H defects are reported. Preserve the full original and repaired arrays plus actual before/after parent forces. A nonconverged or clashing candidate is a failure, not permission to loosen geometry checks.

After inspecting this purely mechanical repair's effect in QM and exterior, root can decide whether a paired electronic test is worthwhile. Construct repaired A first, derive B by the same intact Asp85 rotation, and give both target metals the identical repaired A/B geometries and identically re-prepared fields. Changed MM H coordinates imply a new field identity. Do not compare old-A/new-B or silently combine old endpoint energies with repaired geometries. A later exchangeable-H/water search is a separate chemically consequential branch, not automatic scope expansion.

This proposal addresses a verified source preparation defect, not a fitted route to a desired preference. All active electronic results remain immutable; no new model or repaired structure was evaluated here.
