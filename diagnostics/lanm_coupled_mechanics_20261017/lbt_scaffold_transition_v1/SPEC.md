# Fixed LBT scaffold transition: saved references and MM attribution

## Question and frozen cells

How do the modeled classical contributions and the discrepancy from whole-peptide PBE0 change across the exact deformation already tested? Separate the initial force offset from the change in force during that motion. Two endpoints do not determine a Hessian, unique equilibrium parameters, or stiffness.

Authority is `ROOT_SCOPE_v2.json` (SHA e7493639861f02a5309dda574f7780147d81f2fb8ff209b780eeb06e79ba6a40). The original scope's 197 aggregate was a transcription error: the unchanged XML counts 34 + 61 + 90 + 2 sum to **187**. No term or scientific selection changed.

Use all four exact 256-atom 7CCO cells already frozen by the whole-reference inventory: before_La, after_La, before_Dy, after_Dy. Both metals use the La-selected cartesian_01 before coordinates and the same saved 13017 La-proposed after coordinates. Full reference result SHA 9229f7be and completed whole-model result SHA 64f27e85 remain authoritative. CONFIG pins every physical parent, positions, existing hybrid POINT, ALPB and reference force table. Physical order, charge −2, target255, physical multiplicities La1/Dy6 and conditional frozen-f reference state remain distinct.

Before-Dy is point023, **not** the own-Dy point024 used in old MM audit13122. Its whole-MACE request was energy-only; forces remain unavailable. After-Dy has whole-reference/MACE/ALPB outputs but no complete hybrid POINT. Its new MM query is component-only and cannot supply hybrid work or hybrid force change.

## Permitted computation and finite limits

One CPU, actual Slurm, OpenMM8.5.1 Reference, native threads1, scheduler memory/time defaults, no requeue. All scientific imports and analysis occur inside that allocation. No DFT, MACE, GFN, ALPB, model loading, new geometry, optimization, finite differences, parameter changes or fitting.

For each of four states:

1. Query the four original MM groups and full total: 5 queries / 1 context.
2. Relabel the four retained bonded classes only; query each and unchanged full total: 5 queries / 1 context.
3. Query every retained bonded term intersecting **any** atom of Phe1 (indices0–21) or Ala17 (244–254). Unique union: 34 bonds,61 angles,90 periodic torsion terms,2 CMAP terms =187. Seven batches of at most31 individual groups, with all remaining terms and all nonbonded forces retained exactly once in group31. Query each term and each batch full total:194 queries /7 contexts.

Thus **9 contexts /204 energy-and-full256-force queries per state; at most36 contexts /816 queries total.** Every attempted context/query is recorded. Partial failures remain unavailable; no retry or survivor-only transition. All four inputs must pass before any query. After-Dy requires an admitted after-La common-group control before its queries.

## Unchanged MM model and closure

Reuse the frozen old audit's native parameter inspection, Bond/Angle/Torsion subset copy, geometry descriptors and closure helper. Extend only the CMAP subset: retain all16 native map tables, exact map indices/energies, atom indices and periodic flags in every copy. Preserve both group3 containers (`QM_MM_LJ_generic` and its native exception replacement), groups1/2, masses, all force parameters and zero constraints. Native-to-retained XML term values and ledger decisions must match exactly.

Replay three archived MM results/configurations/executions and POINT MM components. After-Dy has no historical MM receipt: validate its full group sum and its groups0–2 against the same-coordinate after-La calculation, after the existing pinned partition helper proves common physical MM definitions. That proof only normalizes the already-declared inactive target sigma in comparison copies (q=epsilon=0); actual systems/group3 are unchanged.

Inherited closure gates: absolute energy residual ≤1e−8 kcal/mol and maximum Cartesian force-component residual ≤1e−5 kcal/mol/Angstrom. Record measured residuals for archived group/total replay, full group sum, bonded-class sum, relabeled full total, every split-system full total, and term-force sum at **all33 terminal atoms** for each bonded class. These are implementation closure checks, not physical accuracy thresholds.

## Saved reference and force comparison

Use admitted stored POINT component/full gradients and admitted reference force tables as authority. Enforce exact internal gradient sums, exact physical coordinates/IDs/state, shared fixed-FF ALPB and source pins. Do not repeat the failed cross-host cap-scatter replay or relax its original exact gates. Existing failures13446/13455/13460 remain unchanged and separate.

For available comparators report all256 endpoint forces, after−before forces, model−reference errors at both endpoints and their change. Use both original frozen full-array directions from DIRECTIONS SHA6c306212: the initial Phe1-CA force direction and normalized whole-deformation chord. Never regenerate or reselect a direction. Include full256 vector RMS/max-atom errors, per-atom arrays, named Phe1CA3/Ala17C246 vectors, component projection contributions and full region/boundary classification. Small projections do not establish force agreement.

Report native saved PBE0/whole-MACE/hybrid work separately; no absolute energy comparison across Hamiltonians. MM group/class and whole-interaction works are exact same-model after−before differences. Terminal sets are overlapping interaction memberships, not residue energies. Union sums count shared interactions once; each terminus and their overlap are explicitly reported. REPORT.md includes at most15 leading interactions ranked descriptively by the larger absolute La/Dy work with deterministic class/native-index tie breaks; all187 are retained in RESULT.json. This ranking selects no new calculation.

## Outputs, interpretation and execution

Fresh output directory only. Preserve each query's serialized System, group mask, exact positions pin, energy and full256 forces; each state's provenance, all187 interactions, native values/units and actual bond/angle geometry; and RESULT.json, FORCE_TRANSITIONS.json and REPORT.md. A failure report retains all four slots and partial checks. No separate reporting job is needed.

Conditional free termini, vacuum dianion, frozen-f representation, SCF evidence, nonstationary endpoints and absent independent whole-reference derivative qualification remain explicit. No quantum bond-energy attribution, accuracy/preference claim, blanket geometry repair, parameter refit or production promotion follows. Root owns submission after static review; executable command and exact pins are in READY.json.
