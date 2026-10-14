# Inconclusive sources: identities and existing case-specific evidence

14 October 2026. Metadata and saved-report lookup only; no new geometry,
energy, classification, fit, biological label or Slurm submission. Exact source
records and document hashes are in [EVIDENCE.json](EVIDENCE.json).

**The nine requested sources belong to five known PQQ benchmark proteins.
They have admitted, complete current pools. The existing records do not identify
a common preparation failure that explains their inconclusive calls.** Several
well-known historical failures concern different samples and must not be carried
over to these sources.

## Defensible names and reference classes

| Identifier | Local protein annotation | Organism | Existing benchmark class |
|---|---|---|---|
| A0ACD6B9F2 | xoxF1, XoxF1 | *Methylacidimicrobium thermophilum* AP8 | Ln-associated; modeled as La |
| Mmol_1770 | XoxF4-2; NCBI ACT48674.1, PQQ-dependent methanol/ethanol-family dehydrogenase | *Methylotenera mobilis* JLW8 | Ln-associated; modeled as La |
| C5AXV8 | exaF, ExaF | *Methylorubrum extorquens* AM1 | Ln-associated; modeled as La |
| Q9Z4J7 | exaA, ExaA | *Pseudomonas aeruginosa* PAO1 | Ca-associated |
| Q88JH5 | pedE, PedE | *Pseudomonas putida* KT2440 | Ca-associated |

Authority: benchmark FASTA headers in
`/groups/banfield/projects/multienv/corkscrew/supplementary_structures/`:
`ln_binder_positive_controls.faa` lines59/118/25 and
`ca_binder_positive_controls.faa` lines107/119, respectively. The A0AC header
records an accession update from A0A7R8L7W4 after a demerge. Mmol_1770 and PedE
are recorded as tier2; the other three as tier1. Those are inherited curation
annotations, not a fresh literature assessment or measured La/Ca affinity.

The preparation/scoring records retain the same expected classes La, La,
La, Ca and Ca. `conditioned_Ca`, `conditioned_La` and `conditioned_2Ca_v1`
describe prediction inputs. They do not replace those benchmark labels.

## Leads for the exact requested sources

### A0ACD6B9F2, Ca-conditioned sample4

This is a recurring difficult XoxF1 source, not an off-site-PQQ exclusion. Its
current204-atom context is the same expanded probe reused by the Ca125 panel,
with charges Ca−4/La−3, both singlets. Its corrected baseline already remained
wrong after accommodation; the fixed-charge result is now an abstention.
The baseline selected separate Ca and La proposals.

The pre-existing donor screen reports Glu204 O distances2.532/3.388Å and
extra Asp347 O distances2.731/2.531Å; neither triggers its <2.2Å warning.
Thus the known short-carboxylate warning supplies no direct explanation here.
The other A0AC difficult source often mentioned in reports is **La sample4**;
its corrected accommodation recovered a correct call, so its result is not
interchangeable with this Ca-conditioned source. See
`orca_parameter_audit_20261012/expanded_transfer_review/REPORT.md` and the
exact existing geometry-screen row copied into EVIDENCE.

### Mmol_1770, Ca-conditioned sample0

This is XoxF4-2, with a193-atom Ca−3/La−2 singlet context. The corrected
baseline went from correct origin to inconclusive accommodation, selecting
separate Ca and La proposals. Its Glu200 O distances are2.322/3.399Å and
extra Asp362 O distances2.567/2.761Å; the existing screen is unflagged.

Two tempting historical diagnoses concern **other samples**:

- The collapsed-H preparation and later H-only repair concern **La sample4**.
- The unusually short extra Asp362 contact and three-point response experiment
  concern **Ca sample1**, not Ca0.

The separate historical canonical MMOL La search also failed on an oversized,
overlapping optimizer trial despite a valid origin; that is not a present
missing-cell failure. These distinctions are documented in
`accommodation_controls_20260920/FOLD_FAILURE_GEOMETRY.md`,
`accommodation_reference_profile_20260920/REPORT.md`, and
`adaptive_metal_20260922/SEARCH_FAILURE.md`. A sample0-specific structural
explanation remains unresolved by these older reports.

### C5AXV8/ExaF, two-Ca sample4

The five original single-Ca predictions placed Ca about25Å from intact PQQ.
The new sample4 is an actual distinct refold that passed unique-PQQ-site and
spectator-aware preparation: active ion **chain B**, nearest PQQ2.509520Å,
typed CN7, spectator-to-PQQ25.122580Å. Its192-atom context preserves the frozen
Ca−3/La−2 singlet state. The corrected baseline selected the La proposal for
both metals and was already inconclusive.

The old off-site exclusion therefore does not describe this new source.
Two-Ca occupancy in ExaF remains a declared hypothesis; the direct experimental
second-site support used in preparation was strongest for ExaA. The full source
retains the distant ion, whereas the local Hamiltonian omits it after the frozen
spectator-distance gate. Fine geometry or missing global spectator response is
not established as this sample's cause. Exact audit:
`ca_coverage_recovery_20261012/refold/audit/CONTINUATION_REPORT.md`.

### Q9Z4J7/ExaA, two-Ca sample3

This source also passed the complete recovery admission. It is the Q9 sample
whose active PQQ Ca is **chain D**, with chain B the spectator; an assumed
chain-B active selector would inspect the wrong site. Saved nearest PQQ distance
is2.403414Å, typed CN7, spectator-to-PQQ25.197612Å. The192-atom context has the
same state-signature hash as the admitted Q9 La controls. The corrected baseline
selected separate Ca/La proposals and was correct; the fixed-charge candidate
becomes inconclusive.

The existing 1FLG audit supports the N-terminal second Ca site in ExaA. That
motivates occupancy, but does not prove the predicted active-site geometry or
the adequacy of omitting the distant ion's global response. The new sources used
explicit AF3 3.0.4/ten-recycle settings; historical folding settings are not fully
matched, so one cannot attribute a source-family difference solely to occupancy.
See `refold/audit/REPORT.md` and `spectator_preparation/FINAL_REPORT.md` under
`ca_coverage_recovery_20261012/`.

### Q88JH5/PedE, La-conditioned samples3/4

These are admitted original controls, **not** the off-site La sample1 and not
the new two-Ca samples. Both already have measured typed CN7, the same192-atom
Ca−3/La−2 singlet state, and no <2.2Å Glu warning. Their saved Glu221 O pairs
are2.366/2.928Å and2.304/3.061Å; Asp358 nearest O is2.818/2.860Å.

The completed PedE geometry study gives active-ion pLDDT94.74/94.82 and PQQ
mean97.19/97.26. These do not suggest an obvious low-confidence source, but
prediction confidence is not experimental validation of fine donor geometry.
The D+2 homolog is Ser360, so an extra-Asp compression diagnosis does not apply.
The prior study measured these controls' origins; its detailed new-source
proposal analysis should not be misrepresented as a La3/4 motion validation.
Authority: `pede_failure_diagnosis_20261013/geometry/REPORT.md` and its completed
`geometry_local_v2/execution/RESULT.json`.

### Q9Z4J7/ExaA, La-conditioned samples2/3/4

All three are admitted192-atom, Ca−3/La−2 singlet sources with the same state
signature. Existing Glu213 O pairs are2.283/3.019Å, **2.167/3.114Å**, and
2.342/2.971Å. Only **La sample3** triggers the old <2.2Å warning; it is a
diagnostic flag, not a preparation rejection or a new exclusion rule. The
D+2 homolog is Ser352 rather than Asp. This is a concrete distinction to carry
into the source comparison, but does not explain the shared abstention of all
three by itself.

The detailed saved Glu-balancing/boundary/H-bond diagnosis in
`orca_parameter_audit_20261012/q9_response_diagnosis/REPORT.md` concerns
**La sample0**. Its0.8Å boundary, omitted Leu485 packing and lack of compatible
motion-DFT validation are leads, not measured findings for samples2/3/4.

## Common limits of the existing evidence

All nine exact state records retain empty explicit-water inventories, the
oxidized PQQ−3 microstate and matched singlet endpoints. Context preparation
records no new protons/waters or core displacement. This establishes consistent
saved state, not that the approximate protonation or local environment is exact.
The historical H diagnostic found long H-parent bonds even in successful controls;
that shared preparation limitation is not a demonstrated selective defect here.

The two-Ca sources retained every raw heavy atom and both ions in the full source;
their spectator passed the existing4.3Å retained-heavy boundary before omission
from the local model. Neither a passed connectivity/CN gate nor a high pLDDT
validates fine coordination energetics. Root's score analysis and Curie's separate
saved-geometry comparison can test case-specific differences without importing
another sample's failure, changing labels or discarding inconclusives.
