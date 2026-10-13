# PedE PQQ audit: no demonstrated preparation defect

Slurm3356 completed on 2026-10-13. The exact scored cofactors retain their intended atoms, hydrogens, coordinates and formal state. The unusual quinone-region geometry is shared by successful and failed predictions; this audit does **not** identify a justified cofactor repair. Continue the independent fixed-geometry solvent/density comparison before changing protonation, redox state or PQQ coordinates.

## What passed

All nine admitted origins passed the declared checks: all24 PQQ heavy atoms plus three invariant hydrogens, `pqq_ox_3minus_v1` formal−3, unchanged source-heavy coordinates, identical Ca/La coordinates, correct H parents/construction lengths, and no covalent cut involving PQQ. Source-to-local and paired-coordinate errors were exactly zero in the stored inputs. The complete generic-to-CCD heavy-atom graph correspondence passed. These checks cover PQQ preparation and its boundary; they do not validate every surrounding fragment's energetics.

All ten raw sources remain represented. The previously unsupported La-conditioned sample1 has raw metrics but no admitted local geometry; it was not repaired or scored. The analysis reused the completed nineteen-local-geometry result and evaluated all six PQQ instances in three consumed crystal sources. No energies, new coordinates, protonation or optimization were computed.

## The geometry concern is real but not failure-specific

Distances below are from the actual saved structures, in Å. Generic C5–O1 maps to CCD C5–O5, the metal-facing quinone-region bond; generic C6–O2 maps to CCD C4–O4. The graph establishes atom correspondence, not an electronic oxidation-state assignment.

| PedE source | C5–O5 | C4–O4 |
|---|---:|---:|
| Two-Ca sample0 | 1.3392 | 1.3082 |
| Two-Ca sample1, wrong accommodated call | 1.3372 | 1.3104 |
| Two-Ca sample2, wrong accommodated call | 1.3396 | 1.3143 |
| Two-Ca sample3 | 1.3425 | 1.3132 |
| Two-Ca sample4 | 1.3386 | 1.3125 |
| La-conditioned sample0 | 1.3375 | 1.3086 |
| La-conditioned sample1, local unavailable | 1.3396 | 1.3127 |
| La-conditioned sample2, canonical | 1.3444 | 1.3112 |
| La-conditioned sample3 | 1.3368 | 1.3071 |
| La-conditioned sample4 | 1.3394 | 1.3064 |

The available origin/Ca-proposal/La-proposal PQQ bond metrics are unchanged within each source because these proposals do not move PQQ. Therefore the accommodated classification change does not arise from a newly distorted PQQ geometry. A shared cofactor distortion could still affect both endpoints or their sensitivity to neighboring geometry; this comparison does not rule that mechanism out.

The archived crystal comparators differ substantially:

| Source and instance | C5–O5 | C4–O4 | Ring-plane RMS, Å |
|---|---:|---:|---:|
| 1H4I A601 | 1.2363 | 1.1417 | 0.0130 |
| 1H4I C601 | 1.2366 | 1.1418 | 0.0129 |
| 4MAE A602 | 1.2487 | 1.2492 | 0.0054 |
| 4MAE B602 | 1.2626 | 1.2476 | 0.0041 |
| 1FLG A701 | 1.3900 | 1.4164 | 0.1235 |
| 1FLG B702 | 1.4043 | 1.3858 | 0.1309 |

These are coordinate measurements, not ideal bond targets or independent redox observations. The [1H4I primary paper](https://doi.org/10.1016/S0969-2126(01)00148-4) interprets its PQQ as semiquinone. [4MAE](https://www.rcsb.org/structure/4MAE) contains Ce/PQQ/PEG; the historic dry scoring preparation omitted the nearby PEG species. [1FLG](https://www.rcsb.org/structure/1FLG) is an ExaA comparator with no redox assignment made here. The disagreement among stored geometries makes a bond-length-only “oxidized PQQ” repair indefensible.

## Chemical-state qualification remains open

The original Ca, new two-Ca and canonical La folding inputs have exactly the same PQQ SMILES: two quinone carbonyls and neutral carboxylic acids. Scoring deliberately uses the deprotonated PQQ3− hypothesis while preserving heavy coordinates. This is a fold-time/score-time protonation difference, not evidence of a reduced-PQQ template substitution. The official [CCD](https://files.rcsb.org/ligands/download/PQQ.cif) confirms the chemical graph, but provides no ideal coordinates; its model coordinates come from1CRU.

Formal fragment charge does not impose that integer charge or oxidation state on the self-consistent density. Accordingly, the graph/state/mapping audit passes while electronic compatibility and equilibrium geometry remain unvalidated. No new fine bond-length threshold was introduced, and broad connectivity admission is not such validation.

## Decision and limitations

Keep the present cofactor preparation unchanged for the current independent solvent test. Do not hand-shorten a bond, transplant crystal coordinates, add protons or change redox state to repair a class call. The immediate useful test remains the already-planned independent electronic solvent response at saved origins. If that test establishes cofactor-localized strain or density instability, a separate bounded, explicitly state-preserving PQQ internal-coordinate response test could assess it; no such calculation was run or launched by this audit.

This result narrows the diagnosis: no lost atom, wrong H parent, PQQ cap or source-to-core distortion was found. It does not establish that oxidized PQQ3− is the experimentally occupied state, that predicted internal geometry is an electronic minimum, or that the ALPB response is correct.

## Reproducibility and cost

- Exact pins and command: [READY_v2.json](READY_v2.json); compact evidence: [EVIDENCE.json](EVIDENCE.json).
- Actual result: `workspaces/pede_solvent_response_20261013/cofactor_review_v2/execution/RESULT.json`, SHA256 `2401adc5fbfb7aec443e743709547278f55d4725cd6f441b032f041eeb546116`.
- Slurm3356: COMPLETED0:0, one CPU,3s scheduler elapsed. Analysis receipt:1.674528s wall,1.529600s process CPU,219232KiB process peak RSS. Scheduler sampled MaxRSS182068K is a separate measurement.
- Slurm3353: FAILED1:0, one CPU,2s scheduler elapsed. Its reader encountered the already-excluded La1 metadata before writing a result. [Versioned recovery](RECOVERY_v2.md) preserves that failure and changes only unavailable-schema handling. Analysis receipt:0.358119s wall,0.244904s process CPU,78616KiB process peak RSS. No scientific conclusion is taken from that attempt.
- Both attempts together used5 allocated CPU-seconds by scheduler elapsed. No GPU or molecular evaluations. All scientific imports/analysis ran on Slurm; login work was source/hash/JSON inspection and documentation.

The full saved result retains unrounded bond measurements, exact atom correspondence, source occupancy/B factors, H/cap records and every unavailable entry. Production protocols, sources, bands and scores are unchanged.
