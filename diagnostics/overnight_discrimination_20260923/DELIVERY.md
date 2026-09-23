# Nikasha: useful accommodation, with a practical three-source route

**The four-angle MACE accommodation method now has demonstrated classifier
utility and a working experimental PLM command.** The strongest consistent-pocket
benchmark repairs both released wrong calls. A same-pocket ablation shows that
the permitted donor movements themselves add useful decisions. In two actual PLM
proteins, they also relieve compressed donor contacts and substantially reduce
source disagreement. The proteins' biological preferences remain unknown.

All owned jobs and collectors are complete. Jacob received the substantive email
with the physical-response figure; the local relay accepted it. No delivery
confirmation is claimed.

## What is supported

| Comparison | Actual result | Meaning |
|---|---|---|
| Strongest ten-fold-context candidate,225 structures |207 correct /0 wrong /1 inconclusive /17 unavailable|Both released errors repaired;25 consumed protein groups|
| Same pockets/strict solver, before→after accommodation |199/1/8→207/0/1 on208 supported sources|Real movement contributes beyond preparation/calibration|
| Practical envelope,100 declared triples, primary |91 correct /9 unavailable|One failed scalar blocks3triples;6old exclusions|
| Separate qualified restart sensitivity |94 correct /6 old exclusions|Restores all supported median-call fidelity; original primary unchanged|
| Two actual PLM triples |6/6 sources; both median predictionsCa-supported|Usability and physical response; unknown labels|

The practical envelope has four inconclusive individual source/context pairs,
versus one under the strict-tenfold comparison. A separate Ca-conditioned A0A3
probe also becomes inconclusive. Across the91 primary complete triples, source
ranges are not uniformly smaller:39 decrease/52 increase versus same-context
origins. Do not transfer the two PLM examples' large range reductions to every
protein, or call unchanged correct group summaries an accuracy gain over release.

The single failed nativeGFN cell oscillated through500cycles. Two prescribed
same-metal saved-state restarts converged to the same energy within
1.44×10⁻⁷kcal/mol. The origin-seeded value was chosen by policy before either
result. Its separate pooled sensitivity restores three missing summaries with
unchanged geometry and bands. This is not an automatic API retry policy.

## Concrete deliverables

- Existing production fastPQQ, explicitDFT, source graphs, MACE modes, shared-pool
  scoring and archived results were reused and preserved.
- The [experimental API](../pqq_three_source_envelope_api_20260923/COMMANDS.md)
  now accepts explicit complete three-source preparations, with finite plans,
  dry-run, execution, collection and report. v2 generalizes input scope only;
  the actual six-source molecular run is v1. v2 preflights are unexecuted.
- The [physical-response report](../plm_envelope_response_20260923/REPORT.md)
  retains all six sources, donor identities, distances, constraints, components
  and editable SVG/PDF figures. La-selected geometries relieve the short Asp
  contacts; two Ca contacts remain short and seven candidates hit the bound.
- The [protein overlay](../plm_candidate_overlay_20260923/REPORT.md) retains all
 176 original protein rows/34 fields and adds separate candidate fields for only
  the two computed groups.174 remain explicitly unscored; originalDFT and actual
  transcript/genome/phylogeny evidence are unchanged.
- A [concise manuscript draft](../nikasha_manuscript_20260922/PLM_ACCOMMODATION_DRAFT_20260923.md)
  separates the method, demonstrated benefit, PLM application and limitations.
  Vault notes and agent entry points reflect the completed state.

No new DFT, folding, full cohort rescore, default promotion, remote push or
goal-status change occurred. Historical reference IDs and score gauges remain
intact. No binding free energy, probability or physiological metal assignment is
inferred from these descriptors.

## Closed ideas and what they taught us

Adding two more donor angles lowered native energies but did not resolve C5AX;
extra Ca/La composite stabilization canceled and probe separation worsened.
The [six-angle branch](../six_angle_accommodation_20260923/REPORT.md) is closed.

Reusing all old DFT-tested donor geometries,32 new strict scalar calls plus four
exact reuses change differential response by at most0.081kcal/mol. The solvent
term's poorer DFT agreement on those tested paths remains. This closes the
[loose-stopping explanation](../strict_donor_response_20260923/REPORT.md), while
preserving the evidence for nativeMACE donor response.

Earlier [LanM La/Dy transfer](../lanm_series_followup_20260923/DY_TRANSFER_REPORT.md)
failed across actual source structures. It establishes no within-series affinity
prediction. It was not rerun to obtain a favorable result.

## Actual latest execution costs

| Work | New molecular work | Allocated CPU-s | Requested GPU-s |
|---|---|---:|---:|
| Full100 envelope transfer |3408 MACE;1212 GFN attempts,1211 pass|82,389|1,571|
| Two-PLM-triple integration |183 MACE;72 GFN,all pass|8,032|251|
| Strict donor-response check |32 GFN,4 exact reuses|1,216|0|
| Two-start failed-cell diagnostic |2 GFN,both pass|42|0|

These are distinct completed scopes, not the total cost of all historical
development. Full100 includes21 allocations, including technical failures and
empty dependent stages; no successful molecular calculation was rerun for the
metadata repairs. The practical six-source batch took251allocation seconds;
separate preparation functions took16.955s and reused prior protonation. Local
preparation/tests/report CPU and historical reuse costs are additional. This is
batch timing, not isolated single-protein latency or a matchedDFT speed ratio.

Focused actual-artifact tests pass: primary transfer9, practical integration7,
API6 plus its7 historical regressions, PLM join4, donor-response6, restart4,
pooled sensitivity3. No new scientific integration used invented outputs.

## Recommendation and next runnable operation

**Pursue the four-angle candidate as an opt-in three-source interpretation tool.**
It adds useful local response and preserves supported group fidelity. Retain its
individual-source and numerical qualifications; the current standard stays intact.

This read-only command inspects an actual prepared single-protein plan and
launches nothing:

```bash
/groups/banfield/users/jwestrob/conda_envs/lanm_qmmm/bin/python \
  scripts/pqq_three_source_envelope_execution.py dry-run \
  --plan workspaces/pqq_three_source_envelope_api_20260923/preflight_v1/07ab/plan.json
```

The linked API commands show explicit preparation/execution paths. Existing
results should be read or exported without rerunning their chemistry.

Primary transfer: [report](../motion_envelope_transfer_20260923/REPORT.md).
Separately named recovery: [pooled sensitivity](../native_failed_cell_recovery_20260923/POOLED_SENSITIVITY.md).
