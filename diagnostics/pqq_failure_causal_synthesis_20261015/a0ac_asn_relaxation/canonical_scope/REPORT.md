# Canonical terminal-amide rule has complete graph support

**All 25 existing canonical inputs passed the saved identity/graph inventory in Slurm4847.** Each already declares a neutral ASN at the fixed-core `anchor_asparagine` role, with complete mapped sidechain/H identities and an existing terminal χ2. No new proximity cutoff, remote amide replacement, source selection or canonical membership change was used. This supports preparing a uniform discrete-amide test; it does not admit any flipped geometry or demonstrate a score improvement.

## What passed

The actual checks establish paired physical maps, neutral amide microstates, complete named sidechains, terminal bond connectivity, both physical HD21/HD22 parents on ND2, exact physical-to-QM links, and the existing CB–CG χ2 moving group. No terminal-mode moving group intersects a retained or omitted cap anchor. Exact source maps, source preparations, protonated-source pins, states and old search manifests remain in every result row.

The proposed common chemical alternative remains an intact terminal-amide orientation change: OD1/ND2/HD21/HD22 move together around CB–CG with CG on the axis. No O/N renaming, independent H motion or chemistry change follows from this inventory. All actual role residues are ASN; the declared analogous GLN rule was not exercised.

## Existing mode membership, without reselection

All rows below are the actual seed1, La-conditioned canonical sources. Membership is copied from their original paired four-mode tasks; it is not recomputed or selected from a favorable energy. The full source and envelope IDs are retained in the [saved table](../../../../workspaces/pqq_failure_causal_synthesis_20261015/asn_canonical_scope_v1/execution/SOURCES.tsv) and [INPUTS.json](INPUTS.json).

| Protein | Sample | Existing anchor | Terminal mode | Selected for Ca | Selected for La |
|---|---:|---|---|---|---|
| a0a3f2yly8 | 1 | A:ASN285 | A/285/chi2 | No | No |
| a0acd6b9f2 | 3 | A:ASN290 | A/290/chi2 | Yes | Yes |
| a8r3s4 | 0 | A:ASN300 | A/300/chi2 | No | No |
| atq70401.1 | 0 | A:ASN288 | A/288/chi2 | Yes | Yes |
| bbl57595.1 | 2 | A:ASN290 | A/290/chi2 | Yes | Yes |
| c5atj3 | 2 | A:ASN275 | A/275/chi2 | Yes | Yes |
| c5axv8 | 1 | A:ASN275 | A/275/chi2 | Yes | Yes |
| c5b120 | 0 | A:ASN276 | A/276/chi2 | Yes | Yes |
| i0jwn7 | 3 | A:ASN290 | A/290/chi2 | Yes | Yes |
| mmol_1770 | 3 | A:ASN318 | A/318/chi2 | Yes | Yes |
| mmol_2048 | 4 | A:ASN313 | A/313/chi2 | No | No |
| o24759 | 3 | A:ASN295 | A/295/chi2 | No | No |
| p12293 | 4 | A:ASN293 | A/293/chi2 | No | No |
| p15279 | 2 | A:ASN288 | A/288/chi2 | No | No |
| p16027 | 2 | A:ASN288 | A/288/chi2 | No | No |
| p38539 | 0 | A:ASN257 | A/257/chi2 | Yes | Yes |
| q4w6g0 | 2 | A:ASN284 | A/284/chi2 | No | No |
| q60ar6 | 1 | A:ASN289 | A/289/chi2 | No | No |
| q88jh0 | 0 | A:ASN281 | A/281/chi2 | Yes | Yes |
| q88jh5 | 2 | A:ASN308 | A/308/chi2 | Yes | Yes |
| q89gy2 | 3 | A:ASN278 | A/278/chi2 | Yes | Yes |
| q8gr64 | 2 | A:ASN272 | A/272/chi2 | No | No |
| q92wy9 | 4 | A:ASN277 | A/277/chi2 | No | No |
| q9l935 | 1 | A:ASN283 | A/283/chi2 | Yes | Yes |
| q9z4j7 | 0 | A:ASN300 | A/300/chi2 | No | No |

A mapped χ2 is not necessarily an active optimization coordinate. For a row marked **No**, preserving its original four-mode selector would leave terminal χ2 fixed during any later rebased local search, while the originally selected modes could relax. The discrete flip and selection of local optimization coordinates are distinct. No modes were added or reselected.

## Still unrun

**Flip construction, full-source contact admission, rebased q0/Jacobian replay, optimization and all new energies are UNRUN.** The inventory did not measure donor orientation or impose a contact cutoff. Graph support does not establish that a 180° alternative avoids clashes, that every original active mode remains admissible after rebasing, or that the alternative is an equilibrium conformer. No class decision, calibration or production result changed. The original 25-reference calibration and all prior results remain immutable.

## Actual cost and pins

Job4847 completed0:0 on `node-64-768g-1`: one scheduler second on one CPU, TotalCPU00:00.759, batch MaxRSS117824KiB. The application receipt measured 1.3032269770046696s across its own boundaries. There were zero model calls and no coordinate construction. This is inventory cost only, not the cost of any future geometry or energy test.

| Artifact | SHA256 |
|---|---|
| Actual result | `e0f7c3d11e7e92f618dd2f10b83eca3616b5e4a4b9b4840077124d701d1c2ccf` |
| Execution receipt | `873800a268d14c21d77290dca1c6ca44ca1cd4a10f9785cb6abfec95ed7011e9` |
| Source table | `c5634dbea3533b909eb7d8bf9a879c9fe80267507aadb33f40c0bb0d8da4ba69` |
| Frozen input manifest | `8ba06c39543e6601bbbb73d59e4cf6b9eec68a317801a2ec9b4272005e90238f` |
| Existing canonical collection | `ee00fa03aba64c23e6a737fb2b3b4a9f40dc139e58bdb9225e8ed10c9af2b58d` |

[RESULT_POINTER.json](RESULT_POINTER.json) records exact paths, existing counts, submission, costs and qualification scope. This report copies actual saved Slurm results and metadata; no new scientific analysis was executed on login.
