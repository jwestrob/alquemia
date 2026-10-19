# Independent review of the exact declaration and inventory

19 October 2026. **PASS: no specification or inventory blocker found.** This is readiness for adapter review and actual zero-call Slurm preflight, not a runtime result. No scientific computation, module import, model load or submission was performed.

Reviewed immutable artifacts:

- SPEC.md: `b18386f62d9dd04e92f6baeaba10eb1c3a9d58a1e303886f45cd8a246546db91`.
- FEASIBILITY_INPUTS.json: `d94e1cadcfa12e15013cac8e5fbc1679a10c4c566ab64e47ccae4dec09a4cc6a`.
- Existing REVIEW.md remains unchanged at `2fa818b623a5daf8f6c3ca129013ed025b7760521bb3ca36904febd6fb247359`.

All exact path/SHA artifact records in the inventory were independently hash-verified. Literal source and certificate metadata support the declared availability; numerical admission replay remains for Slurm.

## Scope and gates

The declaration freezes both four-cell deformations before new scoring. Its ten-call ceiling matches the inventory: eight full-256 vacuum GFN gradients, one regional-172 isolated GFN gradient and one regional MACE energy/gradient. Both regional calls are solely 7CCO after-Dy. Seven existing regional pairs and all eight ALPB receipts must be reused. No duplicate scout, solvent, embedded-GFN, MM, DFT or optimizer call is included.

Zero-call preflight must cover all eight states before the one-CPU 7CCO before-La whole-GFN scout. Native admission, finite full gradient and measured wall time no greater than 600 seconds gate the remaining consumed-fixture work. Only complete four-cell results with both per-metal work discrepancies no greater than 1 kcal/mol permit the four prebound 7CCN calls. That last rule is explicitly development triage, not a reference uncertainty or accuracy threshold. The mandatory whole-GFN-plus-ALPB comparator needs no extra calculation. Failure stops new scheduling and preserves the denominator.

## Reuse details that the adapter must enforce

The exact 7CCN after-La replacement is NATIVE_RESULT `f40dbdf2…`, joined through recovered POINT `ca3954fd…` and certificate `4707e888…`. The certificate states same Hamiltonian, different solver protocol, old failure still unavailable, and exception limited to this exact isolated component. Its recorded damping is 0.05. Replay that exception's own native provenance and qualification; never rerun the failed default solve or apply its damping to a new state.

For 7CCO after-Dy, use the exact metal-independent after-La map `c470ece5…`, including both existing cap coordinates and two-anchor Jacobians, changing only the declared target symbol/state. Its supplemental ALPB is `26f188d8…`; the historical original inventory's missing solvent entry remains intact. Whole target index 255 and regional target index 0 are distinct. Both charges are −2, with 256 physical versus 172 capped atoms.

Use raw regional MACE forces and isolated GFN gradients for the new correction. The old hybrid gradient or embedded-minus-isolated component cannot substitute. Retain exact map and raw-input authority, the full cap chain rule, original conversions and immutable prior reconstruction failures. Preserve native GFN 6.7.1, parameter text `f6f05c6c…`, corrected profile `635e111b…`, 300 K, accuracy 0.0001, 500 iterations, effective UHF 0 and no restart/solvent/point charges for new calls. Physical MACE La1/Dy6 remains separate.

Implementation must bind these identities, count actual attempts, retain failed slots and gate later stages on actual receipts. In particular, the existing GFN admission helper checks native state, populations, SCC/settings, scalar/gradient consistency and atom order; invoking a source label alone cannot replace those checks. Full force comparisons remain descriptive, and no optimizer follows this finite test.
