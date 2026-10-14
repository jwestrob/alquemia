# Interior ion placement lowers energies but does not establish better discrimination

The fixed midpoint finds lower-energy configurations for the consumed ExaA sources,
but the ambiguous source's Ca-minus-La contrast changes by only **−0.177 kcal/mol**.
Its Ca and La energy gains almost cancel. The control instead shifts **+1.618
kcal/mol**, toward La on this descriptor's sign convention. Close this finite
probe; it has not earned routine deployment or an enlarged optimizer.

## Actual result

Job4137 completed all eight new cells. Both metals used the identical seven-candidate
pool for each source. All20 historical/full-transfer cells were reused exactly;
no source was omitted. Mathematical and operational selections agree.

| Source | Ca selection | La selection | Ca energy change | La energy change | Change in Ca−La |
|---|---|---|---:|---:|---:|
| Control, sample0 | Original adaptive_Ca | Midpoint adaptive_La | 0.000000 | −1.617772 | +1.617772 |
| Ambiguous source, sample3 | Midpoint adaptive_Ca | Midpoint adaptive_Ca | −1.442252 | −1.265282 | −0.176971 |

Energies are kcal/mol relative to each metal's selected historical-three candidate.
The previous five-candidate pool reproduced that historical selection exactly.
Origin-only contrasts are unchanged. The target-minus-control contrast changes
by −1.794742 kcal/mol, mostly from the control's La stabilization. This is not a
classification improvement or a newly independent validation observation. No
bands, aquo offset, threshold refit or affinity claim were applied.

The control's newly selected La configuration is **0.386260 kcal/mol higher** in
MACE vacuum energy, offset by **−2.004032 kcal/mol** in frozen-charge solvent energy.
Thus vacuum-force descent alone would miss this particular composite preference.
In the target, vacuum changes favor La more than Ca, while solvent changes oppose
that differential; their sum is small. Components are retained in RESULT.json.

## What was actually tested

The sole changed coordinate was the exchanged metal position: halfway along the
previous source-derived reciprocal displacement, about0.182Å. Each host retained
its exact saved adaptive_Ca and adaptive_La donor arrangements,192-atom membership,
14caps, empty explicit-water inventory, oxidized PQQ3− state and Ca−3/La−2
singlets. Preparation4136 admitted both sources/all8 cells at the unchanged gates.
This modeled geometry is not an observed structure, equilibrium population or
stationary minimum.

The endpoint descriptor is unchanged MACE-OMOL vacuum energy plus native ALPB
at that endpoint's published-parameter GFN2 vacuum charges. Actual execution used
8 warm MACE calls,32 xTB/tblite SCC calls and16 native component calls. The existing
component-v2 acceptance and matching-water-charge checks passed for all8 cells;
the historical failed absolute crossbackend gate remains failed. ORCA's erroneous
default La shell assignment was not used. Exact code, model, runtime, parameter,
coordinate and state pins are in the frozen CONFIG and RESULT.json provenance.

Actual saved-data replay executed156 checks, maximum residual2.743583e−11kcal/mol
against the unchanged1e−8 algebra tolerance. These are implementation checks, not
an accuracy estimate. No new DFT, GPU work, optimization, folding or labels.

## Cost and restart

Scoring4137:278 scheduler seconds on7CPUs (1946 allocated CPU-seconds), measured
TotalCPU1496.373seconds, batchMaxRSS6621956K. Internal MACE phase114.341seconds;
solvent phase155.214seconds; the total includes setup/collection. Internal elapsed
277.826seconds differs slightly from scheduler timing. Preparation4136:3seconds
on1CPU, measuredTotalCPU2.607seconds, batchMaxRSS52004K. No production speedup is
claimed from this contained diagnostic.

All jobs are terminal and the actual completion event was acknowledged. No rerun
is required. Read `workspaces/exaa_metal_midpoint_20261014/scoring_v1/execution/FINAL_COLLECTION.json`
for the unrounded matrix and `EXECUTION.json` for phase/resource receipts. Exact
commands are in scoring/READY.md and the immutable workspace SUBMISSION.json.

## Decision

Do not promote this probe or automatically refine the line/search. Preserve the
physical lesson: donor-only accommodation leaves some metal-position sensitivity,
but reducing individual endpoint energies does not ensure a useful paired change.
The DFT-confirmed donor response, weak packing effect, zero full-transfer gain and
small target midpoint gain favor examining the solvent/environment representation
before spending more on this coordinate. A separately scoped feasibility audit may
inspect protein shielding; it is not an implemented or scored correction.
