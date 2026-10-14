# What remains inconclusive, and what to test next

**The remaining cases do not share a demonstrated preparation defect.** Three
of the four inconclusives in the recovered Ca-source benchmark sit within
1.2 kcal/mol of their expected-class band. The fourth, XoxF1, is a more substantial
source-geometry problem. Saved donor motions and omitted packing provide specific
next tests; another blanket refold or threshold adjustment is not supported.

This review adds no new classification improvement. It explains the remaining
abstentions of the completed fixed-charge candidate: recovered125 is still
121 correct / 0 wrong / 4 inconclusive. All previous predictions are preserved.

## The four recovered Ca-source cases

Distances below are to the **expected-class** edge of the frozen candidate bands,
in model kcal/mol. They are not uncertainties, affinities or probabilities.

| Protein and source | Distance short of expected band | What the actual saved results show |
|---|---:|---|
| ExaF C5AXV8, two-Ca sample4 | 0.196, La edge | Near-boundary result; both metals select the same La proposal. No optimizer-bound hit. The old corrected scorer was also inconclusive. All four same-conditioning peers and all five La-conditioned peers are La-called. |
| XoxF4-2 MMOL_1770, Ca sample0 | 0.812, La edge | Origin is correctly La-called under its own static calibration. Accommodation shifts the contrast −4.186, chiefly MACE (−4.005). The Ca proposal reaches the0.8Å bound and shortens the nearest Glu200 O contact2.322→2.151Å; the La proposal does not hit the bound. All other available peers are La-called. This identifies a specific differential deformation to qualify, not proven excessive relaxation. |
| ExaA Q9Z4J7, two-Ca sample3 | 1.191, Ca edge | Correct Ca origin becomes ambiguous after a+16.692 accommodation shift, of which+15.496 comes from MACE. It begins with Glu213 O at2.099Å; accommodation lengthens that contact. Active-metal pLDDT84.21 and metal→PQQ PAE7.283Å are worse than its same-conditioning peers. The active ion is chainD, not the spectator onB. All four other two-Ca samples are Ca-called. |
| XoxF1 A0ACD6B9F2, Ca sample4 | 5.815, La edge | The solvent change removes the old wrong-Ca call but does not establish La. Accommodation helps by+12.176. Asn290 OD1 lies4.595Å from the metal versus2.552–2.672Å in Ca0–3 and stays distant in both proposals. The local search cannot fully restore this contact inside its fixed-metal0.8Å displacement domain. However, correctly called La sample4 also has a distant Asn, so this is not a sufficient cause or an exclusion rule. |

Source conditioning is the ion supplied to structure prediction, distinct from
the experimentally curated protein class. The two-Ca sources are the admitted
recovery structures, not their obsolete off-site single-Ca predecessors.

## Five additional inconclusive control structures

La-conditioned **PedE samples3/4** fall1.835/1.425 kcal/mol short of the Ca edge.
Their origins are already inconclusive under the separately calibrated static
model. Both retain typed CN7 and high local confidence; the saved analysis does
not reveal a shared broken fold. One La proposal reaches the motion bound.

La-conditioned **ExaA samples2/3/4** fall1.182/4.005/1.583 kcal/mol short of the
Ca edge. Accommodation shifts them+10.051/+11.530/+5.159, mostly through MACE.
All six proposals reach the motion bound. The correctly called La sample0 also
hits that bound, so the flag alone is not diagnostic. Sample3 alone starts with
the old <2.2Å Glu warning; the other two do not. Sample4's origin is already
inconclusive under the static calibration; samples2/3 have correct static calls.

These cases account for the five inconclusive existing triples: PedE078 and
ExaA096–099. The original six unavailable triples remain unavailable. No new
median/voting policy is introduced, and related samples are not independent
biological observations. The review covers all48 available same-protein
source/context peers; exact identities, old exclusions and triple members are
retained in [EVIDENCE.json](EVIDENCE.json).

## Two different effects must stay separate

1. **Changing solvent models and recalibrating:** eight of the nine targets keep
   exactly their previous selected geometry IDs. The new abstentions arise from
   changed solvent energies and movement of the canonical Ca boundary. ExaA La2
   alone changes selection, and that change helps its Ca margin. This is not a
   newly introduced optimizer failure.
2. **Moving away from the origin within the new model:** the differential works
   above are mostly MACE-driven. Their magnitude asks whether local deformation
   is priced correctly. Origin scores on accommodated bands are only a diagnostic;
   they do not establish a better static classifier or justify per-case fallback.

Neither interpretation proves that the MACE surface, folding, solvation or
omitted scaffold is uniquely responsible. The prior PedE DFT comparison tested
origin solvent response, not these ExaA/MMOL donor movements.

## Concrete next experiment

**Test the cost of the existing ExaA donor motions with their missing local
packing present, before widening the optimizer.** The new full-source contact
lookup verifies an actual omitted neighbor in two-Ca sample3: Glu213 CG approaches
Leu485 CD2 from3.839Å to3.528Å in the Ca proposal and3.512Å in the La proposal.
Asp350 OD2 also changes its contact with omitted Phe442. These are source-derived
contacts, not a hypothetical outer spring or evidence that a clash occurred.
Correct two-Ca peers also approach Leu485, ending roughly3.57–3.67Å away. Thus
this is a shared omitted interaction to test, not a target-only geometry flag.

Proposed bounded comparison: this source plus the first admitted same-conditioning
peer, sample0; all three existing geometries cross-scored for both metals in a
chemically complete expanded context. Freeze the contact-based membership rule
before energies and use it on both sources, preserving each source's state,
waters and caps. That is12 composite cells before reuse; record every underlying
MACE/native solvent call. Compare **within-context deformation works**, never
absolute endpoint minima across different contexts. No optimization is needed.

The nearly equal Leu approach for the two proposals may make its effect largely
common to both metals. If the expanded environment changes both works similarly,
it has not explained discrimination and this branch should close. A useful
differential effect would justify a general context policy and broader frozen
validation, not immediate promotion or fitting these test cases. An independent
reference on the specific MMOL motion is the next distinct question if packing
does not explain the effect. Do not repeat the previously unsuccessful generic
alternative-start or joint-metal searches under a new name.

No new molecular test is launched by this diagnostic report. ExaF remains an
honest near-boundary abstention; the larger XoxF1 donor rearrangement is a separate
structural hypothesis, not an excuse to replace its source with a favorable fold.

## Execution and implementation limits

Three completed one-CPU Slurm jobs4037/4038/4041 used53 allocated CPU-seconds in
total, zero GPU time and zero new molecular evaluations. Score algebra closes
against archived values at the frozen1e−8 kcal/mol tolerance. All48 geometry rows
and their saved three-candidate mappings completed; confidence/optimizer evidence
is retained per source. Exact commands, result hashes, software/helper pins,
process timings and peak memory are linked from [COMMANDS.md](COMMANDS.md).

One reporting defect was caught and preserved:4038's nearest-"nonbonded" field
included fixed PQQ covalent bonds missing from the kinematic graph. It must not
be cited.4041 replaces only that descriptor with explicit inter-residue distances;
donor distances, scores and source coordinates did not change. Proximity alone
is not an energy measurement. All nine target states remain the recorded dry,
oxidized-PQQ−3 state; consistent bookkeeping does not establish optimal hydration.

See the [identity audit](identities/REPORT.md),
[independent interpretation](independent/REPORT.md), and the saved numerical
tables for the complete denominator. Production defaults and calibration remain
unchanged; this is development evidence, not blind validation.
