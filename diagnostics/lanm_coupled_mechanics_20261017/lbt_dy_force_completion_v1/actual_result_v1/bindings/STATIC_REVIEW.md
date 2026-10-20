# Independent static review

**PASS for the declared single GPU job. No remaining source blocker found.** Actual saved-input replay and native admission must still pass within that job before its results are available.

## Frozen bindings

| Artifact | SHA256 |
| --- | --- |
| FORCE_COMPLETION.py | `d77cc40e6cb76f30beb5408c0fafa3d8107c0979f526beb069ee7060ed3a177f` |
| CONFIG.json | `b1655fb8d13b23e036365fc068dfcde2906d83469af693bd94bf526456fb65a6` |
| READY.json | `b64569d37c5a620898f877111485ed8e587832566e684816417a26cfbbd0af79` |
| SPEC.md | `5438ca99d6166e7164fe980e5461614254d6c4e2185712ddcf6fa08db31071ed` |
| SCIENCE_DECISION.json | `7c58e45bab117f93d5182ab48ea2d7501671476dc2f6f05ea3a236f65ed4ef21` |
| INVENTORY.json | `e2c196d1b5961cc74d50647e0565780928c3c1c103b282d563833f4b7e003957` |
| RUN.sbatch | `e8a0f2b21922061d8c07f3333c3e57c741b469245ca2be9bb065c7fa2abb2bfd` |

I read the complete source and native helper, inspected the actual inventory/parent schemas and original force-producing code, parsed source syntax without importing it, and verified the five direct CONFIG artifact pins plus launcher hash. This was metadata/source review, not execution of the scientific gates.

## Admission and call scope

The source freezes eight original cells and exactly three missing Dy force requests. Its saved checks precede lazy model construction: complete collector/config/source joins, native row identity and output/engrad pins, exact physical IDs/order and positions, original task/POINT/model/state authority, each metal's ALPB receipt/raw gradient, frozen directions, and five available MACE force comparisons.

The final schema handling correctly distinguishes the inventory protocol, compact recipe wrapper and task-style state declaration from native `input_state_check`. POINT-backed cells resolve their original full recipe; completed 7CCO supplements retain the pinned whole-model authority.

The unchanged native helper `0e48b8c6` and force configuration `d3c28300` enforce the original checkpoint, float64/head/versions, physical multiplicity, state parity and finite full forces. Runtime checks bind the qualified host, interpreter, 12 CPUs and one allocated GPU. The launcher has no memory/time override; its inherited `standard` header requires root's declared `--qos=standard_expand` submission override. READY's launcher arguments correctly omit the interpreter already executed by the launcher.

Only three serial `compute_force=True` calls are reachable, with one warm model. The existing `1e-7 eV` energy-replay gate, exact native input-state replay and force-sign/shape checks remain unchanged. The first failure stops further requests. Attempted native calls remain counted even if evaluation raises. All eight endpoint slots and four matched comparison slots exist before admission; unavailable entries remain explicit. Each force report links directly to its actual MACE raw receipt and task, including new supplements, while preserving the original energy-only provenance.

## Force arithmetic and interpretation

The reference conversion retains `(raw_ALPB * HA) / BOHR` and reference composite force `-(native_gradient + solvent_gradient)`. Model forces preserve the original producer order separately:

- 7CCO: `raw_MACE_force * EV + (-raw_ALPB * (HA / BOHR))`, matching MODEL_COMPLETION.py `f35e2225`.
- 7CCN: `-(-raw_MACE_force * EV + raw_ALPB * (HA / BOHR))`, matching REFERENCE.py `c0cca6dd` and its saved composite gradient.

No cross-order bitwise equivalence or relaxed tolerance is introduced. Matched vectors are `F_Dy - F_La`; their reported error is the MACE contrast minus the PBE0 contrast. Full vectors, atom-vector norms and only frozen directions are retained. Original energy comparisons are copied without a new energy/work calculation.

Design limitations remain as recorded in DESIGN_REVIEW.md `7f0ff3bb`: physical Dy6 MACE and frozen-4f9/effective-singlet DFT are distinct approximations; this is conditional descriptive mechanics evidence, not spin-matched validation. Shared solvent is not independently validated. Native SCF/state/gradient caveats and false force-accuracy qualifications remain. No optimizer, reference, solvent, other-model or geometry calls are introduced.

No scientific imports, numerical replay, jobs, source/config edits or production changes were performed by this review.
