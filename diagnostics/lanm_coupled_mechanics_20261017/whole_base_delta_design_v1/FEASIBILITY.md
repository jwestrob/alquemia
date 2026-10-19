# Whole-GFN base plus local MACE correction: fixed LBT inventory

19 October 2026. Source and metadata review only: no scientific imports,
energies, force arithmetic, preparation, jobs or changes to existing sources.
The running 7CCN DFT reference14795 is untouched.

## Finding

The eight specified states can use the **same existing baseline172 region rule**
and four exact archived maps. Seven local MACE/isolated-GFN component pairs and
all eight fixed-FF ALPB receipts exist. No compatible full256 vacuum-GFN result
was found for these coordinates. The minimum new fixed-geometry work, if every
call admits, is **eight full256 GFN2 gradients, one local172 GFN2 gradient and
one local172 MACE force evaluation**. The latter two are solely 7CCO after-Dy.
No MM, embedded-GFN, ALPB, DFT or optimization call is needed for this inventory.

[FEASIBILITY_INPUTS.json](FEASIBILITY_INPUTS.json), SHA
`d94e1cadcfa12e15013cac8e5fbc1679a10c4c566ab64e47ccae4dec09a4cc6a`,
contains the exact absolute paths and full hashes for every parent, cell,
component, map, native helper and parameter lineage. It is copied metadata,
not an execution manifest or runtime admission.

## Frozen coordinates and reusable components

| Source / geometry | Metal | Local POINT SHA prefix | Local MACE / isolated GFN | Full256 GFN | Fixed-FF ALPB |
| --- | --- | --- | --- | --- | --- |
| 7CCO before: cartesian_01 | La | `3fac44a7` | Both cached | New | Cached |
| 7CCO before: cartesian_01 | Dy | `6995daa6` | Both cached | New | Cached |
| 7CCO after: hybrid13017 La terminal | La | `2db7d738` | Both cached | New | Cached |
| Same 7CCO after coordinates | Dy | Unavailable | Both new | New | Cached supplement13142 |
| 7CCN before: geometry_20 | La | `5894bdbe` | Both cached | New | Cached |
| 7CCN before: geometry_20 | Dy | `b1f0ea87` | Both cached | New | Cached |
| 7CCN after: cartesian_03 | La | `ca3954fd` | Both cached; certified isolated solver exception | New | Cached |
| 7CCN after: cartesian_03 | Dy | `32faa4ab` | Both cached | New | Cached |

The 7CCO inventory is
`lbt_whole_reference/INVENTORY.json` (`48d9f40f…`); the 7CCN inventory is
`lbt_whole_7ccn_reference_v1/INVENTORY.json` (`4203d1aa…`). These are two different
already-declared motions. They are not interchangeable after geometries.
7CCO reference evidence is consumed development evidence; the independent
7CCN reference remains pending, and 7CCN is Lu-conditioned.

The complete local matrix is
`workspaces/lanm_coupled_mechanics_20261017/lbt_hybrid_pool_recovery_v1/execution/RESULT.json`
(`d98aafb6…`). The original partial matrix remains unavailable and is not promoted.
The original whole-reference inventory still correctly records missing after-Dy
solvent. Its later receipt is separately bound from
`lbt_whole_model_completion/result/COLLECTION_RESULT.json` (`64f27e85…`),
`rows[cell_id=after_Dy].ALPB` (`26f188d8…`), not the old `solvent_admission` entry.

The bounded cache audit inspected 138 full256 vacuum-GFN TASK records representing
116 physical records in the current `lanm_coupled_response_20261016` and
`lanm_coupled_mechanics_20261017` workspaces. None matches a target's exact
ordered IDs, elements and coordinates. The earlier electronic-substitution
matrix covers geometry_00–17; it does not supply geometry_20 or either new
Cartesian/hybrid endpoint. Among 390 local hybrid POINT records, no Dy POINT
matches the exact 7CCO after coordinates. This is an audit of these owned
workspaces and known native routes, not a claim about every historical file.

## One partition, exact archived caps

Use `lbt_partition/REGION_POLICY.md` (`7589391b…`) and
`lbt_partition/PREPARE_V2.py` (`0d02c79f…`). The original source-union rule retains
complete residues3–14, preceding residue2 C/O and following residue15 N/H,
with two source-derived H caps. Baseline172 contains170 real atoms including
the target and two caps. Whole256 contains no synthetic caps, waters or field.
Both whole and local charges are −2.

All seven cached states have the identical ordered core identity/boundary
topology. Each archived map's physical atoms and POINT positions exactly match
its reference physical JSON. No distance analysis or map regeneration was used
to establish this metadata identity.

| Physical geometry | Exact archived map SHA |
| --- | --- |
| 7CCO before | `312a4e8403065199950b152f0110a4401b7330da26d6cecf1ca8ab99a38b3342` |
| 7CCO after | `c470ece5cf05f36736a5b84ee7093d3d05a9fd73e3c44e104ab83a624b3e211c` |
| 7CCN before | `5d11b9e5450f6440dbd102398441ae42187df30cc49ef150249d287045609e16` |
| 7CCN after | `bcc6387e68860bd17bddaaa971814b36c2d5712a011559d659d206dc02581b82` |

For missing 7CCO after-Dy, use the exact after-La map, whose target is the
metal-independent placeholder M: local index0 maps to physical255. Substitute
the native target symbol/state only. The cap anchors, cap coordinates and both
anchor Jacobians already exist; reconstruction tolerances are unnecessary.
The map retains field metadata, but the proposed isolated calculations use no
external charges. `PREPARE_V2.scatter` accepts a zero field-gradient array and
applies the two-anchor chain rule without a sign flip.

## Native route and state distinctions

Existing full256 native precedent:
`diagnostics/lanm_coupled_response_20261016/LBT3_SOLVENT.py`
SHA `e3a906103b2982aed5f99b4bc01f4fd02fdf0986f8611ba1863987a8df3a4373`.
Its `admit()` reads native engrad, populations, state, method and SCC evidence;
job12510 admitted8/8 original whole-peptide vacuum/water cells. This precedent
does not make its old coordinates reusable or authorize new water calculations.
Its all-in-one `run()` is not the proposed panel: it hardcodes original sources
and vacuum/water pairs. Reuse the isolated native staging/admission functions.

The dimension-general native helper is `lbt_hybrid/GFN_ADMISSION.py`, SHA
`3563fb16f5bb986420176e78020afd1dd593e5cb69ef36f15fa1fe3021b99dab`.
It can admit isolated256 or172 native gradients with explicit count/state/task
fields. Existing component producer `lbt_hybrid/EVALUATOR_V2.py`, SHA
`686eba84fd3dfcf28724d2d4e7957528c987445d7137b18c88a07dfec055f661`,
shows the local native call and arithmetic; its full hybrid evaluator should
not be invoked for the new formula.

- Pinned executable: `/home/jwestrob/miniconda3/envs/lanm_bench/bin/xtb`,
  SHA `481a1e7d8151fc08e0341b8e06ab60ac09cf1c2628738eeb83346271274bf8f`.
- Full native parameter text: `param_gfn2-xtb.txt`, SHA
  `f6f05c6c5264cb2f75c54f8ac6838044541c1e493b16f1a4d4ee4fa64de5f545`.
  Preserve the complete text, exact staged loading and runtime-library pins.
- Corrected La/Dy profile: `635e111bac31e7c5755bf612c8e148ce73b401c110ff31887fd5e31b590fa85f`;
  qualification12488 result `7b2d23d6…` admits all8 native/standalone controls.
  This is parameter/backend parity, not independent force or physical-state validation.
- Unchanged standalone settings: GFN2, charge−2, effective UHF0,300K,
  accuracy0.0001,500 iterations, no restart, native analytic `--grad`, no solvent
  or point-charge argument. Whole256 has750 effective valence electrons;
  local172 has538. Physical La multiplicity1/Dy6 remains separate metadata.
  Neither GFN model is an explicit open-f Dy calculation.
- Local MACE retains q−2 and physical multiplicity1/6, the exact native
  float64 model/runtime and force sign. Its pinned native helper is
  `LBT3_FD.py` (`0e48b8c6…`); the exact model/config evidence is in the inventory.

## The existing SCC exception must remain visible

7CCN after-La is exactly the earlier isolated point059 failure. Do not call that
failed default solve again. The already-admitted replacement is
`lbt_hybrid_scc_recovery_v1/execution/failed_005/NATIVE_RESULT.json`
(`f40dbdf2…`), joined through recovered POINT `ca3954fd…` and
`SOLVER_EXCEPTION_CERTIFICATE.json`
(`4707e888b429a2edbea8864ea1b785feee1dbf6a9fdda629529bf455a36280c7`).
It uses the previously qualified broydamp0.05 numerical exception with the same
Hamiltonian/state. Preserve original failure, qualification13009 and certificate;
this is not permission for damping fallback at any new state. The later closed
embedded-response failure is a different calculation and is absent from this
proposed formula.

## Assembly boundary and unresolved gates

The proposed fixed-state energy is
`GFN2_vac(full256) + MACE_vac(local172) - GFN2_vac(local172) + fixedFF_ALPB(full256)`.
No MM ledger, GFN embedded-minus-isolated term, external field, extra dispersion
or second solvent term belongs in that sum. The full base retains electronic
scaffold/environment interactions; the isolated local correction does not itself
provide high-level environment polarization. Caps remain an approximation.

Use native-energy before→after differences before conversion, retaining each
component and solver lineage. Raw MACE forces are eV/Å and require a minus sign
for gradients; GFN and ALPB raw derivatives are Hartree/bohr. The local correction
requires the archived cap chain rule on both terms. Saved full hybrid gradients
cannot substitute for this new total derivative. Saved MACE component gradients
exist; a separately mapped isolated-GFN gradient is not a POINT component—the
old POINT stores embedded-minus-isolated. Raw isolated gradients and exact maps
are available for a new, explicitly recorded assembly.

Prior cross-host cap scattering differed at last bits, so source identity alone
does not promise bitwise reproduction of old assembled arrays. This inventory
adds no reconstruction tolerance and does not claim new force qualification.
Keep archived native inputs/maps authoritative; record any new assembly and
numerical checks separately rather than weakening old failed replay gates.

The meaningful next decision is a fixed eight-state comparison against the
unchanged whole-DFT reference, with complete force/energy denominators. New SCC
failures remain unavailable and stop the declared route; no survivor selection,
recovered old optimizer or favorable-sign criterion follows. Full-gradient
qualification would still be required before any optimization. Conditional
terminal chemistry, Lu-source provenance, frozen-f/SCF limitations and absent
binding-free-energy interpretation remain. Existing hybrid after-Dy is still
missing even if a new whole-base-delta value is later obtained.
