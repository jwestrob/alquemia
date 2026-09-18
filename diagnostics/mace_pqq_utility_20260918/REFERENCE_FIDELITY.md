# Frozen PQQ reference fidelity

| Set | DFT correct / total | Masked MACE correct / total | MACE unavailable |
|---|---:|---:|---:|
| Canonical calibration | 25/25 | 25/25 | 0 |
| Consumed crystal transfers | 3/3 | 2/3 | 1 |

All 56 archived DFT outputs/receipts were checked and final energies reparsed.
The frozen masked-MACE two-call values, algebra and original bands were replayed.
No new energies or threshold fitting were needed for this reference inventory.
The new full-panel timing/reproduction experiment is reported separately.

1KB0 remains unsupported by whole-chain MACE because its available preparation
cannot represent the deposited connectivity/cofactor chemistry faithfully.
The fixed-core DFT result remains valid within its recorded scope. This is a
coverage limitation, not a hidden failed classification or a MACE success.

1H4I and 4MAE repeat sequences in the canonical panel. Calibration accuracy is
not prospective accuracy; sequence groups are not independent protein families.
Motif and composition already separate the canonical panel. The task here is
retaining practical functional-class discrimination, not proving general affinity.
The two methods' numerical scores use distinct scales and references.

Full joined rows, exact source hashes and replay result:
`workspaces/mace_pqq_utility_20260918/accuracy_v1/{scores.tsv,result.json}`.
Replay: `accuracy.py --root REPOSITORY --output NEW_OUTPUT_DIRECTORY`, with explicit
executable absolute-path examples in COMMANDS.md. No scientific executor is used.
