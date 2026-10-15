# Parallel Ca4 selected-energy reference — staging v2

Declared before binding, numerical checks or native evaluation. This is a separately versioned **scheduling-only revision** of `../reference_candidate_v1/`; that unbound/unrun source package stays immutable. Root now permits this independent local-energy check to proceed alongside uniform-transfer job5661 and reporter5666. No final transfer outcomes existed or were inspected when this decision was made.

The four selected-energy comparison independently tests the already observed incremental MACE vacuum response (+17.988344487 kcal/mol) and is informative under either eventual transfer classification outcome. It neither depends on nor becomes a qualification prerequisite for the uniform transfer. Endpoint choice, method, coordinates, source/state gates, numerical criteria, interpretation and four-cell scope are unchanged. Only the former terminal-full-transfer execution prerequisite changes.

## Design decision and four cells

The earlier reference-feasibility note proposed origin → a common La-generated geometry. That tests total accommodation, but can agree while the alternative basin remains above the old accommodated state. This separately declared design instead tests the vacuum contribution to the **actual selected-energy change** for A0ACD6B9F2 Ca-conditioned sample4:

| Metal | Old selected geometry | New selected geometry | Charge / multiplicity |
|---|---|---|---|
| Ca | adaptive_Ca | relaxed_flip_Ca | −4 / 1 |
| La | adaptive_La | relaxed_flip_La | −3 / 1 |

Exactly four fixed singlepoints; source `a0acd6b9f2-pqq-la_model__conditioned_Ca__seed-1_sample-4`, envelope `2eda778295ffed58`. All retain the current 204 atoms / 14 caps, dry oxidized PQQ3−, original protons and saved coordinates. No La4 extension, origin DFT, cross-metal extra endpoints, gradients, optimization, state change or GFN2 parameter file. Each metal traverses its own old/new selected geometries; **this is not a common perturbation, DFT pool minimum, isolated Asn effect, affinity or classifier test**.

Actual source: `asn_relaxation_remaining_v1/execution/FINAL_COLLECTION.json`, SHA256 `2ed05f5a241f4f1eccbb2a680428a4b032a51449df56ad323e138b83cf32fd82`; actual preparation v2 SHA256 `6405de39b9d3ffa88aea7f99f895a685ddcd8ec5b1b37529766b72c365d9d57d`. The full archived common6 paired geometries remain admission evidence even though only the four selected endpoints receive DFT. Source maps, both original and rebased proposal lineages, exact MACE receipts and accepted full-source gates remain pinned separately.

The selected old XYZ hashes `696cbb38224149c2034679da39fe269a9f423d7571bf26314cd7eef64d261a9e` (Ca) and `1aa9804090aa28469303003cc577193fcd7f724991f632456009dbacd06943f4` (La) match their actual proposal/MACE-request bytes. They were additionally checked against the same ten historical DFT manifests named in the earlier feasibility note, with no exact match. This closes that specific omission from the earlier 16-hash cache audit; it does not claim an exhaustive archive search. New selected hashes are `975422eeaf5d909efc81608046a498faca3cccb0851a0d798cd84a22783b7590` and `bbf2aa9a03d1ca00f8e4d0bf7f3934763f76a167103d638adca506e3e29a5aa4`. Both saved vacuum and composite common3/common6 selection records support these choices; no selection is made by this adapter.

## Method, checks and outputs

Use the unchanged pinned ORCA 6.1.1 native r2SCAN-3c / DefGrid3 / TightSCF / NoAutostart vacuum protocol. Copy the already executed MMOL native task runner, allocation renderer, original scalar/population parser, corrected printed-TRAH collector and count-aware collector into any future immutable workspace. No mutable script import or inherited MMOL193/13/charge assumption supplies this source.

Metadata binding only copies/hashes exact inputs and code. All geometry, electron counting, parser fixture replay, execution and collection require Slurm. Preflight checks the actual 204/14 graph/state, unchanged physical inventory, both original proposal and relaxed proposal coordinate/receipt joins, fixed caps, full six-candidate Ca/La coordinate pairing, current selected choices, singlet all-electron/46-electron La-ECP inventory, and native input bytes. Replay the actual inherited parser fixtures, all eight ExaA and all six MMOL completed native collections. Preserve original energy-based convergence admission, actual printed solver criteria, density residuals and electronic-state warnings; acceptance is not electronic stability or force qualification.

One endpoint per explicit execution command, 12 MPI ranks on one standard Slurm node, one thread each; memory comes from the unchanged allocation-aware renderer with its reserve. No memory/time/priority override, retry, automated expansion or job submission by this adapter. Existing task locks and native runner refusal of incomplete/stale attempts remain intact. Per-attempt preflight, native execution/output/GBW, collection and failure receipts are preserved. A final collection always retains all four declared cells; a missing endpoint leaves its work and the paired change unavailable.

After the parallel-stage root decision and actual preflight, the two old-selected endpoints serve as independent scouts. The two new-selected endpoints additionally require root's immutable review of both actual admitted scout collections, including convergence/state/runtime. This is the same four-cell denominator, not a six-cell origin/proposal protocol. A failed scout calls for diagnosis; no automatic iteration escalation or retry is permitted.

Report ΔE_M = E_M(new selected) − E_M(old selected), then D = ΔE_Ca − ΔE_La, for native DFT and saved MACE, plus individual DFT−MACE discrepancies. Use the existing 1 kcal/mol interpretive resolution, neither a fitted acceptance threshold nor a DFT error bar. Do not turn positive D into a pass: the tiny saved Ca work may be unresolved while the strong La lowering is supported. The saved MACE increments are Ca −0.286035604, La −18.274380091, D +17.988344487 kcal/mol; this is most of the saved composite +20.883598044 response, not independent validation of its solvent part. No bands, row-minimum selection or experimental labels enter the native comparison.

## Internal execution gate and costs

Every molecular execution requires an immutable root decision file binding the **new v2 manifest** and the actual completed canonical25 and ten-source A0AC solvent collection/execution pins below. The adapter checks the canonical25 availability, ten-source 40/40 completion, exact producer config joins and A0AC MACE lineage. It must explicitly request `execute_four_selected_Ca4_native_endpoints`, name execution stage `parallel_with_pending_uniform_transfer_v2`, record `transfer_status_at_decision: "pending"` and `transfer_outcomes_inspected: false`, and give the scientific interpretation. Root may additionally record the actual pending/running5661/5666 state and its observation time. Those statements describe the immutable decision time; later scout execution does not invent or require a still-pending scheduler state.

| Completed evidence | SHA256 |
|---|---|
| `asn_canonical_solvent_v2/execution/FINAL_COLLECTION.json` | `1506f040eace36620a53a1efbda94638f5cec444910097d0016b359af462b438` |
| Its `EXECUTION.json` | `df010cf7b2922e1dc8bc4782d7e0d42228946ca1bf24b0043eb52945c124cd41` |
| `asn_relaxation_solvent_v1/execution/FINAL_COLLECTION.json` | `64668c7ec758da42d34b039855bdc370b39d06815db47a00dfaf69c41a69368d` |
| Its `EXECUTION.json` | `338de13151763a0a110664b477d841abea63665f7fb169036f07f07c3f0efc83` |

All paths are under `workspaces/pqq_failure_causal_synthesis_20261015/`. No terminal full-transfer collection/execution is supplied, generated, inferred or required. No desired class, favorable response sign or transfer improvement gates admission. Missing/negative root authorization blocks native calls. The two-new-cell scout review remains unchanged and requires both actual admitted old-selected scout outputs. This is internal scientific staging, not fresh user permission; this source package creates no decision file.

The scientific protocol identifier remains v1 because no scientific protocol changes. Manifest schema/staging identify the v2 adapter. `PREPARE.py` records exact hashes of all five original v1 files as parent evidence, adds the three completed-evidence pins (the A0AC collection was already pinned), and copies the same native helpers. `REFERENCE.py` replaces only the execution-decision gate and adds its evidence check to preflight; sources, parser fixtures, physical checks, scout gate, native execution and numerical collection are byte-identical function bodies. `RUN.sbatch` is unchanged.

Readiness commands are proposed only. This delivery is source preparation only: **binding, scientific preflight and all four molecular cells are unrun**. The source package does not submit jobs. Agreement would independently support this finite selected-energy improvement; disagreement would motivate a local surface/state diagnosis. Neither validates solvent, basin populations, biological preference or universal physics.

Real planning precedents: the six MMOL193 cells used 12 MPI ranks and 3683–5497 wall seconds each (jobs4603/4604/4637/4638/4639/4640: 47736/45468/49116/55020/65964/44196 allocated CPU-seconds). The ExaA192 eight-cell reference used 522696 allocated CPU-seconds, 3670–6994 seconds per endpoint. These are actual allocation costs, not A0AC runtime forecasts. See `../../mmol_reference/FINAL_REPORT.md` and `../../../exaa_donor_response_20261014/FINAL_REPORT.md`.
