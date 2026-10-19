# Independent static review: fixed 7CCN whole reference

Reviewed 19 October 2026. **PASS for the frozen source and configuration; actual Slurm preflight and native execution remain pending.** No scientific module was imported, numerical replay performed, or job submitted in this review. Only this review file was written.

## Frozen evidence

| Artifact | SHA-256 |
| --- | --- |
| REFERENCE.py | `c0cca6dd4f44619fc2e3201f024b7913f6b7e9d8e8c4fca694fe17bd942d0ecf` |
| SPEC.md | `a8cb688c8d5f8f1f109a6c1436da24e68410cb92528f503d123cc83650307447` |
| CONFIG_PREFLIGHT.json | `876b5a5665a52c8a5ff4345a7318a743245ea61f797e172f2e15fd3ab2fa675e` |
| INVENTORY.json | `4203d1aabfc4f312b8627469e0b94450927e06b73d38294b3db6a5a6739cb0e2` |
| READY.json | `d804ac9686945724418dd0e4ec32f44a86c0b005c5a56a184936aaaecb30dc43` |

All artifact pins in CONFIG_PREFLIGHT and READY were independently hash-checked. The old whole-reference helper is `c22c989674e280806abace63c67148b0181caee0b7b20bcd85b8e3d75d6c2a86`; the compact native manifest is `18ed7c0f4141b95cba20191bf36867bc839ae8722945b81dfb74fe3eeeda7f25`.

## Findings

- **Selection and interfaces:** Literal saved schemas match the source joins: job 12968's exact 7CCN-La start and terminal, common-pool geometry_20/cartesian_03, four original POINT children, preparation and physical inventory. Both metals use each same coordinate array. The adapter uses the old whole helper's low-level state, deck and solvent functions; it avoids the old hardcoded 7CCO preparation and comparison paths. The fixed charge is −2, the 256th atom is the target, and no caps, waters or external field are introduced.
- **Electronic state and native admission:** Whole state construction preserves physical La multiplicity 1 and Dy multiplicity 6 separately from effective singlets and ECP cores 46/55. The compact parser accepts the supplied whole state, rather than its own charge −1 preparation recipe. It retains version, charge, electron count, RHF, ECP, analytic-component, finite-gradient and coordinate/energy checks. Native receipt verification binds executable, manifest, template, output, gradient and runtime assets. `expected_engrad_complete` is an existing runner field. Canonical JSON conversion preserves exact saved-row replay across tuple/list serialization.
- **SCF meaning:** The original PBE0-D4/RIJCOSX/DefGrid3/VeryTightSCF/PModel deck and native stopping semantics are unchanged. Printed residual failures remain observations under native admission. Neither normal termination nor a full analytic gradient establishes numerical derivative consistency, electron binding or a qualified physical state. The closed Mex external-point-charge derivative result is not transferred to this nuclear-gradient panel.
- **Solvent and force reuse:** All four ALPB receipts exist and are checked through their native command, input, charge, component and gradient evidence. Saved whole gradients retain `-f * EV + g_ALPB * (HA / BOHR)`; reference solvent gradients retain `(g_ALPB * HA) / BOHR`. The source does not impose false bitwise equality after reassociation. Force signs are negative gradients. Both La force arrays are retained; both archived Dy MACE forces remain unavailable. A full reference force table and the predeclared full-deformation chord prevent a small selected projection from replacing the complete evidence.
- **Finite execution and collection:** All stages require Slurm. Preflight makes no molecular call. A run selects one task in a fresh directory and invokes the pinned runner once. Remaining cells require replay of the actual admitted before-La scout and a pinned measured scheduler-cost decision. Missing or failed slots remain unavailable, and all four references are necessary for work comparisons. The final complete-parent check requires the same preflight and membership in its phase allowlist.
- **Resources:** Molecular execution calls the original allocation helper, which requires 12 actual MPI slots on one node and derives memory from observed AllocTRES; the runner receives one worker and 12 ranks. Parent RSS and application timing are explicitly distinguished from MPI peak memory and scheduler cost. No launcher is included in this bundle: root must bind the existing 12-MPI molecular and one-CPU saved-data wrappers and retain their actual allocation receipts.

## Interpretation and remaining checks

The comparison forms per-metal endpoint differences before conversion, adds the same solvent contribution to both methods, and reports signed work errors and the Dy-minus-La work difference. The 1 kcal/mol resolution is descriptive, not an accuracy or native-admission gate. This is a new reference outcome for a consumed, Lu-conditioned peptide source and an existing nonconverged deformation, not an independent-protein validation, affinity result or fitted model.

No static blocker was found. Actual zero-call Slurm replay must establish array arithmetic, prepared states and task signatures before the scout. Native admission, resource behavior and actual cost remain unmeasured for these new cells; subsequent continuation remains gated on their own evidence.
