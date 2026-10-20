# Independent review of continuation scout 15254

## Verdict

The saved result supports completion of the declared 7CCN–La continuation: 12 new distinct composite evaluations and 12 accepted steps, reaching cumulative counts of 92 accepted steps and 92 evaluations within the frozen 92/120 limits. All requests completed; every saved trial passed the original physical and state guards. The search stopped at its accepted-step limit and remains nonconverged.

This review inspected saved fields, native receipts, source bindings and file hashes only. It performed no numerical replay, scientific imports, molecular calls or scheduler actions. In total, 215 referenced file/hash pairs were verified.

## Exact evidence

| Artifact | SHA-256 |
| --- | --- |
| Actual run `execution/RESULT.json` | `c445cf4b622fa2dbc7d7c47456e4b20b7c8db1a0b7718365fab2310c3c90a8c6` |
| Frozen `CONTINUE.py` | `0c14a03e0c3d13884739d8be33a6c2f2c39f734c5374a36de0f4d511d7573b8e` |
| `CONFIG_RUN.json` | `7127a7fcba909e71a0ec3675a89159eecf627e23ce61def59e8aebc6935a14f7` |
| Actual preflight15192 result | `2dd9db660a7907077872b35a50b71e571fa2782ab5a12fa50200bdc3deadd245` |
| Root run decision | `a21d2c53f7dc2f5c6598f12462ce0b468370d3a9a1b1e9e1d6edfecf6a6fae32` |
| Saved `SEARCH.json` | `4ac43f90cb4827115915bc0b11c25b9b845c0d9d7cce3e1471e23da1c778931d` |
| New `CHECKPOINT_FINAL.json` | `317c2b24eff75a2af4a6fb7066f1754d1df1b9a284405a8aad2ea3b20d5c74cd` |
| New terminal `evaluation_012/POINT.json` | `eb96df292085e05f3eed1788e506c54babab646c372698a41443705945fa2a1b` |
| New terminal cached ALPB result | `37e905285b5dc549039824ff2ab45bec5c71949232bd37ddf67e3390b476fffd` |

The actual workspace is `workspaces/lanm_coupled_mechanics_20261017/lbt_7ccn_continuation_scout_v1/run/`. The scheduler receipt records COMPLETED, exit 0:0, 26 seconds and 12 CPUs on the qualified H200 host. The producer records 24.575708642 seconds elapsed.

## Scope and admission

- The original checkpoint, initial point, prepared source, history, trust radius and 80/80 consumed counts were preserved. The geometry source was not rebased. Original maximum/RMS force stopping thresholds remain 0.10/0.03 kcal/mol/Å.
- Each of the 12 new requests joins exactly to its POINT and trace entry, with one MACE call and one fixed-charge ALPB call. There are no new DFT, GFN, MM, opposite-metal or other-source calls and no cache hits counted as new evaluations.
- All MACE forces have 256 rows, with the exact charge −2, La multiplicity 1, target index 255, physical atom inventory, model head and recipe. Each ALPB execution returned zero and retains the expected native binary, physical coordinates, charges and raw-gradient artifacts.
- The 99 original archived records remain unchanged. Two additional saved-only aliases of the original origin/start have different replay-directory paths but byte-identical positions, gradient and MACE children; they are not molecular calls. The subsequent 12 records are the actual new requests.
- All 12 steps accepted their first trial with the original Armijo and physical guards. Final geometry/state admission passes, with no reported bond, chirality, peptide or overlap violation. Old and new physical metadata are identical apart from coordinates.

## Outcome and limitations

The stored incremental work from the old checkpoint is −1.705988079 kcal/mol MACE vacuum, −2.520989587 kcal/mol fixed ALPB and −4.226977666 kcal/mol composite. These are producer values, not a new calculation in this review.

Forces did **not** improve monotonically. For example, the printed maximum falls from 9.739698 to 3.379517 and later rises to 9.153297 kcal/mol/Å. The final maximum/RMS are 5.674648/1.086248 kcal/mol/Å, both above the unchanged convergence criteria. Physical-state qualification, force-accuracy qualification and scientific-accuracy claims remain false.

No blocker was found to the separately declared single PBE0 energy/gradient at the exact new terminal, using the unchanged whole-256 La recipe and reusing the old terminal14880 native result plus both cached ALPB endpoints. That test should compare new-minus-old native energies before conversion and assess the incremental model lowering. The existing 1 kcal/mol reporting resolution is neither a native admission threshold nor an uncertainty estimate. Agreement on this segment would not erase the earlier 7CCN–La work discrepancy, qualify the total gradient numerically, select an electronic state, establish preference or authorize another optimizer extension. Native SCF residual, electronic-state, vacuum-charge and fixed-solvent limitations remain explicit.
