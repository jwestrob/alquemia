# Existing-evidence benchmark set for the East River paper

## Delivered

Built a source-linked ledger of **68 evidence records**, acquired/inventoried **28 published structures**, recovered the existing **B2/C5 native sequences and AF2 models**, and prepared **five new site states / ten La–Ca endpoint inputs**. No new energy calculations or cluster jobs were run in this construction phase. Baseline defaults and prior results remain unchanged.

These are not 68 independent labeled tests. Calibration cases, related proteins, mutations, structural replicates and unresolved candidates remain identifiable. No new accuracy estimate or claim of broad discrimination follows from constructing the set.

- Human-readable ledger: [BENCHMARK.tsv](BENCHMARK.tsv); compact citations, measurements and qualifications: [EVIDENCE.json](EVIDENCE.json).
- Pinned artifact index and counts: [RELEASE.json](RELEASE.json).
- Full evidence, source captures, coordinate inventories and prepared inputs: `workspaces/benchmark_set_20260915/`.
- Full machine ledger: `workspaces/benchmark_set_20260915/release_v1/benchmark_manifest.json`.
- Ready task manifest: `workspaces/benchmark_set_20260915/ready_tasks_v4/manifest.json`.
- Operations: [COMMANDS.md](COMMANDS.md); authorization/scope: [AGREEMENT.md](AGREEMENT.md).

## What the ledger contains

| Package | Records | Interpretation |
|---|---:|---|
| Historical PQQ | 29 | 25 calibration cases, three already-inspected crystal transfers, one secondary crystal not run |
| Additional PQQ evidence | 18 | Functional controls, direct PqqT comparisons, six-blade controls/B2/C5 and unresolved candidates |
| Non-PQQ | 15 | Direct, condition-limited and supporting comparisons; unresolved claims retained |
| Lanthanide-binding proteins | 6 | Lanmodulin/Lanpepsy comparisons and candidates lacking verified La/Ca labels |

Current dispositions: **32 records reuse existing scores; 2 have new prepared inputs; 9 have structures but preparation/interpretation gates; 8 have sequences but structure/identity gates; 2 have apo models requiring holo preparation; 1 has evidence but an input gate; 14 remain unresolved or unlabeled.** These are record counts, not independent-protein denominators. A prepared record can contain several ordered sites or source geometries.

The ledger retains source captures and hashes, conditions, construct/site mapping qualifications, biological groups, conservative fold/family holdout groups, prior score exposure, existing scores and failures. Exact sequence matches are inventoried; absence of an exact match does not establish independence by homology. Eligibility must use evidence type and limitations, never `direction` alone. None of these newly assembled cases is represented as a fresh blind test without an exposure audit.

## Immediately prepared additions

| Record | Frozen structural choice | Site states | QM atoms | Crystal waters |
|---|---|---:|---|---|
| Hans-LanM | 8DQ2 author chain A; EF1 → EF2 → EF3 | 3 | 50 each | 0 each |
| Bovine alpha-lactalbumin | 1F6S Ca-conditioned and 6IP9 La-conditioned strong site | 2 | 40 / 43 | 2 / 3 |

All five use the **existing** `generic_peptide_amide_vertical_native_r2scan3c_v3` policy. Every source heavy atom is preserved (observed maximum displacement 0 Å); paired La/Ca nonmetal coordinates match exactly; endpoint charges differ by one. No PQQ threshold or universal zero is inherited. All source-conditioned geometries and water inventories were retained before energies.

**Hans-LanM:** first author chain selected, all three primary sites retained in sequence order; the weak fourth Na site is not relabeled a high-affinity site. Biological dimers are A–C and B–D. Standard-residue protonation uses pH 5, matching the reported conformational assay; heavy coordinates remain from the pH 7 crystal. The protein-level La/Ca response includes folding/dimerization and does not supply three independent site affinities. Full deposited coordinates were used for preparation, but each QM core excludes the surrounding protein. [Primary study](https://doi.org/10.1038/s41586-023-05945-5).

**Alpha-lactalbumin:** the primary study supports stronger La than Ca binding at the strong Ca site, but its complete conditions and exact assay-to-structure construct correspondence remain unresolved. Preparation uses the existing generic pH 7 default, not an asserted reproduction of the 1996 assay. The 2019 La-only ITC experiment does not supply the missing matched Ca measurement. Both source geometries are supporting tests, not two independent labels. [Original comparison](https://doi.org/10.1016/0167-4838(95)00223-5), [La-bound structure/ITC](https://doi.org/10.1038/s41598-018-38024-1).

## Most valuable further computational targets

1. **B2/C5 and six-blade PQQ controls.** Exact primary-paper gene IDs map to local full native sequences (467/440 aa) and existing sequence-identical AF2 models. All five ranks were inventoried; unrelaxed rank 1 was selected deterministically, with its relaxed counterpart retained. The models are apo: PQQ/metal placement and a declared six-blade preparation protocol remain necessary. Signal-peptide trimming/tag state of the assayed construct is not fully recovered. The preprint reports same-assay La/Ca competitive retention at the protein level; this is not a site Kd. These targets are especially relevant to extending the PQQ discriminator beyond its canonical architecture. Ca-supported six-blade structures 1C9U and 6JWF provide an opposing functional class, not matched direct affinity measurements. [B2/C5 primary preprint](https://doi.org/10.64898/2026.07.10.737718).
2. **Seven additional PQQ functional controls.** Primary La/Ca functional directions and sequences are available. Exact experimental-construct reconciliation and holo input preparation remain gates. FAM1 XoxF5 is deliberately excluded from a strict La-positive class: the inspected purified-enzyme evidence is Ce-supported. Functional use must remain separate from affinity.
3. **PqqT WT/K142A/K142D.** Same-study La/Ca ITC is valuable, especially as a mutation series. All three belong to one family. 9B1U contains Na, not a La/Ca binding core; the sparse Gd sites in 9B1V are not automatically the assayed event. Metal-site mapping is the blocker. [Primary study](https://doi.org/10.1073/pnas.2405836121).
4. **Lanpepsy and Mex-LanM.** Lanpepsy now has a mapped La structure (9VY7), resolving the old missing-structure gap, but its five adjacent/shared-donor metals need a coupled-site preparation. Mex has mapped Nd/Y geometries, with the weak fourth crystal site separately identified; source-metal geometry and construct qualifications must remain explicit.
5. **Aqualysin.** The primary low-affinity-site comparison favors Ca over La, but mapping that event to either 4DZT crystal site is unresolved. Coordination number is not a substitute for experimental site identification. [Primary study](https://doi.org/10.1271/bbb.66.1281).

ConA remains conditional qualitative exclusion. Parvalbumin remains cross-study support. Aequorin remains an ordered vector associated with an experimentally unmapped site. Calbindin opposing-site labels are unverified; the 1KSM structural citation is [Bertini et al.](https://doi.org/10.1023/A:1012422402545), not evidence of opposite La/Ca affinities. Different-condition CD2 graft comparisons and candidates with only within-lanthanide evidence cannot inflate the direct La/Ca denominator.

## Validation and provenance

Ten tests passed (0.456 s, no skips): six new real-artifact integrity tests plus four existing baseline regressions. Coverage includes metal inventory, paired-coordinate/output-path invariants and grouping, rejection of an explicitly corrupted real capture, historical exposure/unresolved labels, exact coordinates on preparation replay, missing-result collection, archived energy/reference algebra and receipt validation. The existing runner dry-run passed for all ten inputs. These are software/preparation checks, not newly executed quantum validation. A collector bug for campaigns lacking aequorin/parvalbumin was fixed: empty groups are skipped, and new Hans/alpha observation groups retain ordered vectors. The unrun collection contains null scores with explicit incomplete status, never invented zeros. Historical baseline regressions still pass.

Original alpha preparation reports recorded a transient builder hash without preserving its source bytes. The final preparations therefore replayed carving/amide repair from the **exact archived protonated coordinates**, with preserved implementation copies. Both La/Ca XYZs are byte-identical to the first preparations; no new hydrogen placement or scientific selection occurred. Earlier staging bundles remain with explicit supersession notices. Final `ready_tasks_v4` preserves implementation bytes and takes family/evidence labels from the curated ledger. Hans retains its actual one-shot preparation script and pinned dependencies/results.

No new ORCA energies, folds, environmental corrections, threshold fits, comparator analyses or robustness campaigns were run. Preparation/acquisition CPU time was not instrumented; no cluster allocation was requested. Previous baseline job 1199299 cost remains 807 s × 64 allocated CPUs = 51,648 core-s for 26 endpoints, not a claimed cost for these new inputs.

## Recommendation for the paper

Retain the baseline and use this release as the expanded, stratified benchmark specification. Score the prepared additions as supporting grouped cases, then prioritize B2/C5/six-blade holo preparation and the additional functional PQQ controls. Report canonical class transfer separately from broader affinity/compatibility evidence. Scores on unlabeled East River proteins remain predictions; they do not establish in situ occupancy, physiological metal use or a causal explanation of transcription.

The complete strategy—including inexpensive comparators, family holdouts, preparation robustness and future experimental validation—is saved in the requested vault note: `agent-captures/2026-09-15_laca-paper-benchmark-existing-evidence.md`. Those later components do not block using the curated set now.
