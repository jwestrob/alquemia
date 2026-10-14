# Independent interpretation: two effects, no common preparation diagnosis

Read-only review of completed job4037 and its frozen plan/source,14October2026.
No molecular calls, numerical reanalysis, submissions, recalibration or changes
to the scoring pipeline. Root owns score interpretation and Curie owns geometry.

**The fixed-distribution candidate has useful gains, but the remaining nine
inconclusives do not justify one universal correction.** They are nine source
contexts across five benchmark proteins, not nine independent biological tests.
All have complete admitted pools. Their reference classes describe existing
metal association; Ca/La conditioning describes prediction inputs.

## Keep the two comparisons separate

**Baseline versus frozen-distribution candidate:** eight targets retain the same
selected Ca and La geometry IDs. Q9Z4J7 La sample2 is the exception; its changed
selection lowers the Ca-minus-La contrast by1.0073kcal/mol and therefore helps
the expected Ca class. New abstentions cannot generally be attributed to a newly
different optimization. The energy-model change at fixed selections and movement
of the frozen class bands explain those changes. In particular, the candidate's
Ca boundary moves4.3986kcal/mol toward more negative contrasts; this contribution
to the lost Ca margins must remain visible.

**Origin versus accommodation within the candidate:** many shifts are dominated
by the MACE component. ExaA two-Ca sample3 shifts+16.6917kcal/mol, comprising
+15.4957 MACE and+1.1961 frozen-charge solvent. Mmol_1770 Ca sample0 shifts in the
opposite direction:−4.1861 total,−4.0049 MACE and−0.1812 solvent. These are real
model responses, but neither direction proves an erroneous motion. They do not
by themselves explain why changing the solvent model changed a classification.

The same-band origin comparison is a diagnostic only. PedE La samples3/4 and
ExaA La sample4 appear Ca-like if their origins are read on the accommodated
bands, but their origins are already inconclusive on the separately calibrated
static bands. Conversely, ExaF two-Ca sample4 has a La-supported static call even
though its origin lies in the gap on the accommodated bands. Selecting whichever
policy supplies the expected answer would invent a new classifier.

## What the individual patterns support

- **ExaA/Q9Z4J7:** four affected contexts, including the new two-Ca sample3.
  The latter has a large differential accommodation response and is a useful
  mechanical test case. La sample3 alone has the previously recorded short-Glu
  warning; it cannot explain all four. The historical missing-Leu485/boundary
  observation belongs to La sample0, not automatically these sources.
- **PedE/Q88JH5:** the two La-conditioned abstentions do not erase the genuine
  corrections of two-Ca samples1/2. They show conditioning/source dependence.
  High local confidence and complete inventories do not validate fine energetics;
  neither is there evidence here for a shared broken fold or missing cofactor.
- **XoxF4-2/Mmol_1770:** Ca sample0 is La-supported at its origin and becomes
  inconclusive after accommodation, only0.8119kcal/mol below its correct band.
  This makes it a complementary motion-response case, but the collapsed-H and
  compressed-extra-Asp diagnoses from other samples must not be transferred.
- **ExaF/C5AXV8:** two-Ca sample4 remains an existing abstention, now only
  0.1956kcal/mol short of its correct band. That distance is not an uncertainty
  estimate or evidence warranting a threshold adjustment. Accommodation moves
  it toward La on fixed bands rather than causing a wrong-direction response.
- **XoxF1/A0ACD6B9F2:** Ca sample4 remains a substantial structural/model outlier.
  Its origin is Ca-like even under the static calibration; accommodation moves
  toward La but no admitted common candidate gives La on the pooled bands.
  Merely selecting another currently scored pose cannot resolve this source.

The identities audit establishes no common preparation defect. All nine retain
empty explicit-water inventories and the declared oxidized PQQ state; consistency
does not establish that these states are physically complete. Peer medians are
descriptive, and the five abstaining existing triples remain actual regressions.

## Highest-value next physical test

Prioritize **omitted short-range scaffold resistance to the existing donor
motions**, conditional on the saved-geometry review identifying actual contacts
or an incompletely represented covalent/packing restraint. Do not infer missing
resistance solely from a large MACE energy gain.

A small useful test would keep the existing three geometries and both metals
for ExaA two-Ca sample3 plus the first prepared same-condition peer by source
index. Define one expanded, chemically complete local context from actual
source contacts to the moved groups, using existing contact/mapping rules before
scoring. Hold this membership and proton/water state fixed over the whole pool.
Evaluate those six cells per source with the existing warm MACE and qualified
cheap solvent machinery; no optimization, new fold or high-level sweep is needed.
Reuse the smaller-context values. Compare each metal's within-context deformation
work and its Ca-minus-La difference, never absolute totals across the two regions.

This can test whether explicitly represented packing changes the claimed
differential accommodation work. A positive result would motivate a general
boundary/context policy; it would not supply a calibrated new classifier or a
full scaffold free energy. Do not add arbitrary springs, copy one favorable
neighbor, or assign new bands from these targets. If no concrete omitted contact
is found, that test lacks its rationale: instead qualify one actual donor
displacement with an independent, finite electronic reference before spending
on more search freedom. Existing PedE DFT checked origin solvent response, not
these ExaA or Mmol motion energies. Neither route is launched by this review.

The current evidence does not justify another scalar solvent factor, routine
reprotonation/water expansion, a different policy per protein, or promoting
the new candidate solely because its wrong-call count is zero.

## Evidence

- Actual result: `workspaces/pqq_inconclusive_review_20261014/scores_v1/execution/RESULT.json`;
  SHA256 `2a6c3419703f9ff7aea0f54eee734530185ccd2e3647e78b26757ace08cb50ab`.
  All numbers above are quoted from its saved rows, without new calculation.
- Frozen scope: [PLAN.md](../PLAN.md), scorer [SCORES.py](../SCORES.py).
- [Identity and exact-source limits](../identities/REPORT.md).
- Upstream completed candidate collection SHA256
  `82981de62dfa2f2ff5eaedd273e38d8b20bf36034dc39110c104908599470d25`.
