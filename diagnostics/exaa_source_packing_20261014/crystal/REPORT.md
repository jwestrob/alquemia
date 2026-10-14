# 1FLG preparation: chain A admitted; chain B unavailable

Job4093 completed the declared two-chain preparation. **A passes; B fails the
unchanged catalytic-Asp partner gate (3.511 Å versus a3.5 Å maximum).** Both
chains remain in the denominator. No threshold was relaxed and no new energy,
optimization or folding calculation belongs to this branch.

## Source and chemical-state audit

The [1FLG deposition](https://www.rcsb.org/structure/1FLG) provides a2.60 Å ExaA
dimer. The previously completed alignment maps each582-residue crystal chain to
Q9Z4J7 positions35–616 with100% identity. This supports the same resolved sequence,
not identical full-length constructs or strain provenance: the primary article
names ATCC17933 whereas the current accession mapping names PAO1. Q9 positions
1–34 and617–623 are absent from the crystal construct.

Actual saved-source audit found no alternate conformers or PQQ external covalent
connection. Each chain has one24-heavy-atom PQQ and two monatomic Ca, each with
deposited occupancy1.0. The raw active/spectator sites are A751/A761 and B752/B762;
PQQ instances are A701/B702. All4563 retained heavy atoms per selected chain have
explicit source→mapped identities and unchanged coordinates. The exact previously
audited CCD→AF3 PQQ name map passes element, graph and donor-set equivalence.

The source redox state remains unqualified. Existing saved cofactor metrics show
long carbonyl-labelled bonds/nonplanarity; these were neither repaired nor used
to invent a new threshold. The primary full-text PDF could not be retrieved
(HTTP403), so no detailed redox discussion is attributed to it. This preparation
tests the existing **conditional oxidized-PQQ3− model**, not a verified matched
experimental oxidation state. See [source evidence](SOURCES.md).

## Actual completion and local-state checks

Job4088 generated both standard-protonated structures, then the original AF3
wrapper rejected one missing terminal heavy atom in each. The original manifests
truthfully retain internal-completion0/terminal-completion1. V3 reused those exact
outputs without repeating hydrogen placement. Slurm verified in both chains that
the sole addition was OXT on mapped A:TRP616 (crystal582), with no removed or moved
deposited heavy atom and no nonstandard H on either Ca or PQQ.

| Final chain A check | Actual result |
|---|---:|
| Expanded context |192 atoms,14 caps |
| Ca/La charges and multiplicities |−3/−2, both singlets |
| Local state signature against frozen Q9 |Exact match |
| Ca/La physical maps |Identical |
| OXT to nearest retained context heavy atom |41.494340 Å |
| Terminal residue to nearest context heavy atom |37.545263 Å |
| Spectator Ca to nearest context heavy atom |24.543982 Å |
| Active/spectator Ca to nearest PQQ heavy atom |2.731787/25.282441 Å |

The terminal residue is absent from context atoms, retained bonds and caps.
The existing spectator boundary check passes independently. The complete protein
identity is precisely the mapped mature subset, while its full-protein hash is
correctly different from Q9. The exact local state key is
`dbdead29335579d38decb6e1b85ae542a7b16390e8f98e90fb01225309aca561`.

Chain B passed the sole-OXT/deposited-heavy checks but stopped before a context
was admitted because its mapped catalytic Asp–Arg partner is beyond the existing
3.5 Å gate. Its final context/state/terminal-distance admission is unavailable;
chain A's checks must not be attributed to B.

## Retained limitations and handoff

The dry policy excludes73 chain-A and33 chain-B crystal waters. For admitted A,
nine omitted waters lie within the existing descriptive4.3 Å radius of retained
context heavy atoms; no atom from the other crystal chain does. Source waters,
full dimer and spectator remain in immutable source evidence. The local energy
model includes neither those waters nor distant spectator/global dimer response.
This is an explicitly dry structural comparator, not the intact crystallographic
chemical environment.

Root owns later scoring of admitted A. B stays unavailable; no automatic retry,
gate relaxation, crystal-coordinate repair or production promotion follows.

## Reproducibility and cost

Actual result: `workspaces/exaa_source_packing_20261014/crystal_v3/execution/FINAL_PREPARATION.json`.
Manifest SHA256: `360455889a1063f1ea577f08ae0f71fbf74215220be95ccc1d9e5d652e96894a`.
Exact artifact pins and both failure/admission paths are in [EVIDENCE.json](EVIDENCE.json).

Both jobs used one standard-partition CPU on node-64-768g-5. Scheduler elapsed:
4088=20s;4093=4s, totaling24 allocated CPU-seconds. Actual analysis receipts:
4088 wall19.474620s/processCPU18.990454s/peakRSS226488KiB;
4093 wall3.825200s/processCPU3.393746s/peakRSS217664KiB.
V3 made zero new protonation calls. All structural preparation/checks ran on
Slurm. Login work was code, metadata and saved-result reporting only.
